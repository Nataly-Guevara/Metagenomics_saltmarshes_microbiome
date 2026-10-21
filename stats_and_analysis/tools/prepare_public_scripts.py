#!/usr/bin/env python3
"""Create a separate, script-only publication copy from a local Chapter 2 collection.

The source collection is read only. No analyses, commits or pushes are run.
"""
import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ACCOUNT_PATH = re.compile(
    r"/(?:scratch/[^/\s\"'<>]+/[^/\s\"'<>]+|"
    r"home/(?:[0-9]+/)?[^/\s\"'<>]+|Users/[^/\s\"'<>]+|"
    r"g/data/[^/\s\"'<>]+/[^/\s\"'<>]+)"
)
EMAIL = re.compile(r"[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}")
SECRET = re.compile(
    r"(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]+|"
    r"-----BEGIN [^\n]*PRIVATE KEY-----)"
)

def relative(value):
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError("Unsafe relative filename in local collection")
    return path

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-project", required=True, type=Path)
    parser.add_argument("--expected-scripts", type=int, default=182)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    stats = repo / "stats_and_analysis"
    collection = (stats / "chapter2").resolve(strict=True)
    publication = stats / "scripts"
    if publication.exists() or publication.is_symlink():
        raise SystemExit("Public scripts folder already exists; nothing overwritten.")
    source_project = args.source_project.resolve(strict=True)
    parts = source_project.parts
    identifiers = []
    if len(parts) >= 5 and parts[1] == "scratch":
        identifiers = [parts[2], parts[3]]
    elif len(parts) >= 3 and parts[1] in ("home", "Users"):
        identifiers = [source_project.parent.name]

    def sanitise(text):
        text = text.replace(str(source_project), "/path/to/metagenomics_project")
        text = ACCOUNT_PATH.sub("/path/to/local", text)
        text = re.sub(r"(?m)^(\s*#PBS\s+-P\s+)\S+", r"\g<1>YOUR_PROJECT", text)
        text = re.sub(r"(?m)^(\s*#PBS\s+-M\s+)\S+", r"\g<1>YOUR_EMAIL", text)
        text = re.sub(r"((?:scratch|gdata)/)[A-Za-z0-9_.-]+",
                      r"\g<1>YOUR_PROJECT", text)
        for item in identifiers:
            replacement = "YOUR_PROJECT" if item == identifiers[0] else "USER_ID"
            text = re.sub(r"(?<![\w])" + re.escape(item) + r"(?![\w])",
                          replacement, text)
        return EMAIL.sub("YOUR_EMAIL", text)

    with (collection / "_collection/copied_files.tsv").open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["kind"] == "script"]
    if len(rows) != args.expected_scripts:
        raise SystemExit("Unexpected script count; publication was not created.")
    prepared = {}
    grouped = {}
    changed = 0
    python_count = 0
    shell_count = 0
    r_count = 0
    for row in rows:
        group = row["group"]
        if not re.fullmatch(r"[A-Za-z0-9_-]+", group):
            raise ValueError("Invalid analysis folder name")
        src = (collection / relative(row["collection_path"])).resolve(strict=True)
        if collection not in src.parents or not src.is_file():
            raise ValueError("Script is outside the local collection")
        raw = src.read_bytes()
        if hashlib.sha256(raw).hexdigest() != row["sha256"]:
            raise ValueError("A source script changed after collection: " + src.name)
        original = raw.decode("utf-8")
        if SECRET.search(original):
            raise ValueError("Possible credential requires manual review: " + src.name)
        cleaned = sanitise(original)
        if ACCOUNT_PATH.search(cleaned) or EMAIL.search(cleaned):
            raise ValueError("Personal path or email remains: " + src.name)
        for item in identifiers:
            if re.search(r"(?<![\w])" + re.escape(item) + r"(?![\w])", cleaned):
                raise ValueError("Account identifier remains: " + src.name)
        original_name = relative(row["source"])
        if original_name.parts[0] not in ("scripts", "betaNTI"):
            raise ValueError("Unexpected source script directory")
        short = Path(*original_name.parts[1:])
        target = Path(group) / short
        if target.as_posix() in prepared:
            raise ValueError("Duplicate public script filename")
        if src.suffix == ".py":
            ast.parse(cleaned, filename=src.name)
            python_count += 1
        if src.suffix in (".pbs", ".sh"):
            result = subprocess.run(["bash", "-n"], input=cleaned,
                                    universal_newlines=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode:
                raise ValueError("Shell syntax check failed: " + src.name)
            shell_count += 1
        if src.suffix == ".R":
            r_count += 1
        prepared[target.as_posix()] = cleaned
        grouped.setdefault(group, []).append(short.as_posix())
        changed += cleaned != original

    def put(base, name, text):
        path = base / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    temp = Path(tempfile.mkdtemp(prefix=".public_scripts_", dir=str(stats)))
    try:
        for name, content in prepared.items():
            put(temp, name, content)
        for group, names in sorted(grouped.items()):
            title = group.split("_", 1)[1].replace("_", " ")
            text = "# " + title + "\n\n"
            text += "Read the [setup instructions](../../README.md) before execution.\n\n"
            text += "Pilot, backup and versioned scripts retain their names; a filename alone does not establish the final version.\n\n"
            text += "\n".join("- [{}]({})".format(n, n) for n in sorted(names)) + "\n"
            put(temp, group + "/README.md", text)
        temp.rename(publication)
    except Exception:
        if temp.exists():
            shutil.rmtree(temp)
        raise

    readme = """# Chapter 2 analysis scripts

Public source code for the metagenomic and community assembly analyses, organised by analysis.

This publication contains scripts and setup documentation. Inputs, outputs, sample metadata, logs and private collection reports are not included.

Personal filesystem paths, account identifiers and email addresses are replaced by generic placeholders. Biological labels, analysis parameters and software versions are retained.

## Before running

1. Replace the /path/to/... placeholders with your own project, database, tool and library locations.
2. Replace YOUR_PROJECT in PBS directives with an allocation you can use; adapt modules and resources to your system.
3. Install the software and packages referenced by each script and supply its input files.
4. Review all input and output paths. Relative paths refer to a working data project, not to these publication folders.
5. Place scripts called by other scripts in their expected locations (usually scripts/ or betaNTI/ in the working project), or update those calls.

This is a source-code collection, not an automatically executable workflow. Scientific analyses were not rerun. The existing scripts may create or replace result files when run.

## Analysis groups

| Folder | Scripts |
| --- | ---: |
"""
    for group, names in sorted(grouped.items()):
        readme += "| [{}](scripts/{}/) | {} |\n".format(group, group, len(names))
    readme += """
## Coverage and validation

All supplied script versions are retained, including pilot and backup versions. Generating scripts for the soil metadata summaries and the matched 16S–MAG turnover comparison were not available in the reviewed collection.

Python and shell/PBS syntax was checked during publication preparation. R syntax was not independently checked by this helper. Input availability, package compatibility and scientific reproducibility still require verification before execution.

The full local collection remains separate and is excluded from Git.
"""
    put(stats, "README.md", readme)
    root_readme = """# Metagenomics Saltmarsh Microbiome

Source code and setup documentation for saltmarsh metagenomic and community assembly analyses.

| Folder | Purpose |
| --- | --- |
| [metagenomics_project/](metagenomics_project/) | Preprocessing workflow structure and initial quality-control script. |
| [stats_and_analysis/](stats_and_analysis/) | Chapter 2 publication scripts, grouped by analysis. |

The Chapter 2 publication contains scripts with personal paths and computing account identifiers replaced by generic placeholders. Read the [setup instructions](stats_and_analysis/README.md) before execution.

Inputs, outputs, sample metadata, logs and private collection reports are not included. Scientific analyses were not rerun during publication preparation. Pilot, backup and versioned scripts retain their names.
"""
    put(repo, "README.md", root_readme)
    print(json.dumps({
        "publication_scripts": len(prepared), "analysis_groups": len(grouped),
        "scripts_sanitised": changed, "python_syntax_checked": python_count,
        "shell_syntax_checked": shell_count, "R_scripts_not_parsed": r_count,
        "personal_identifier_scan": "passed",
        "inputs_or_outputs_copied": 0, "scientific_analyses_rerun": False,
        "next_step": "Review and stage only stats_and_analysis/scripts and both README files."
    }, indent=2))

if __name__ == "__main__":
    main()
