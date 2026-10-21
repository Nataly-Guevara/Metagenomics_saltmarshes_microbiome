#!/usr/bin/env python3
import ast
import csv
from pathlib import Path

base = Path("results/ecological_strategies")
source = Path("scripts/build_community_capability_matrix.py")
output = base / "Comment02_01_MAG_capability_prevalence.tsv"

# Read the definitions from the ORIGINAL script.
# Do not import it, because importing would execute its analysis.
tree = ast.parse(source.read_text())
definitions = None

for node in tree.body:
    if isinstance(node, ast.Assign):
        if any(isinstance(t, ast.Name) and t.id == "GROUPS"
               for t in node.targets):
            definitions = node.value
            break

assert isinstance(definitions, ast.Dict)

groups = {}
for key, value in zip(definitions.keys, definitions.values):
    group = ast.literal_eval(key)
    filename = ast.literal_eval(value.elts[0])
    fields_node = value.elts[1]
    fields = ast.literal_eval(fields_node.func.value).split()
    groups[group] = (filename, fields)

all_records = []
reference_ids = None
redox_counts = {}

for group, (filename, fields) in groups.items():
    path = base / filename

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    ids = [r["Final_MAG"] for r in rows]
    assert len(ids) == 373
    assert len(ids) == len(set(ids))

    if reference_ids is None:
        reference_ids = set(ids)
    else:
        assert set(ids) == reference_ids

    for capability in fields:
        assert capability in reader.fieldnames

        calls = {}
        for r in rows:
            value = r[capability].strip()
            assert value in ("0", "1", "0.0", "1.0")
            calls[r["Final_MAG"]] = int(float(value))

        count = sum(calls.values())

        all_records.append({
            "GRESP": group,
            "Capability": capability,
            "N_MAGs_positive": count,
            "N_MAGs_total": 373,
            "Prevalence_pct": round(100 * count / 373, 4),
            "Source": str(path)
        })

        if group == "Redox":
            for mag, call in calls.items():
                redox_counts[mag] = redox_counts.get(mag, 0) + call

assert len(all_records) == 37
assert len(redox_counts) == 373

with output.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(all_records[0]),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(all_records)

print("=== COMMENT 2.1: MAG DENOMINATOR AUDIT ===")
print("MAGs per source:", len(reference_ids))
print("Curated capabilities:", len(all_records))
print("Strategy groups:", len(groups))

print("\n=== ALL 37 CAPABILITIES ===")
for r in all_records:
    print("{:<14} {:<48} {:>3}/373  {:>6.2f}%".format(
        r["GRESP"],
        r["Capability"],
        r["N_MAGs_positive"],
        r["Prevalence_pct"]
    ))

n_low_redox = sum(v <= 2 for v in redox_counts.values())

print("\n=== REDOX FLEXIBILITY CHECK ===")
print(
    "MAGs encoding 0-2 of 10 redox capabilities:",
    n_low_redox,
    "/ 373 =",
    round(100 * n_low_redox / 373, 2),
    "%"
)

print("\nSUCCESS")
print("All six groups use the same 373 MAGs")
print("Output:", output)
