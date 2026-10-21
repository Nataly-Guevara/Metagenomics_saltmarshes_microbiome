#!/usr/bin/env python3
"""Collect scripts, small metadata files and a file inventory for local review.

Uses the Python standard library. Source files are never changed and nothing
is uploaded. An existing output ZIP is never replaced. This is an inspection
bundle, not a complete reproducible analysis package.
"""

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import stat
import sys
from datetime import datetime, timezone
from zipfile import ZIP_DEFLATED, ZipFile


CODE = {'.py', '.r', '.sh', '.bash', '.pbs', '.rmd', '.qmd', '.ipynb'}
SMALL_DATA = {'.csv', '.tsv', '.txt', '.md', '.json', '.yaml', '.yml',
              '.toml', '.ini', '.cfg', '.xlsx'}
SKIP_DIRS = {'.git', '__pycache__', '.ipynb_checkpoints', '.venv', 'node_modules'}
FILE_LIMIT = 10 * 1024 * 1024
TOTAL_LIMIT = 100 * 1024 * 1024


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scripts-dir', action='append', default=[])
    parser.add_argument('--metadata-dir', action='append', default=[])
    parser.add_argument('--inventory-dir', action='append', default=[])
    args = parser.parse_args()
    source = args.source.resolve(strict=True)
    output = args.output.absolute()
    if not source.is_dir():
        parser.error('--source must be a directory')
    if output.exists() or output.is_symlink():
        parser.error('Output already exists; choose a new ZIP filename.')

    roots = {}
    for mode, values in [('scripts', args.scripts_dir),
                         ('metadata', args.metadata_dir),
                         ('inventory', args.inventory_dir)]:
        for value in values:
            relative = Path(value)
            if relative.is_absolute() or '..' in relative.parts:
                parser.error('Scan directories must be paths relative to --source.')
            root = source / relative
            if not root.is_dir() or root.resolve() != root or root.is_symlink():
                parser.error('Missing or symbolic scan directory: ' + value)
            if root == output.resolve() or root in output.resolve().parents:
                parser.error('Place the output ZIP outside all scanned directories.')
            roots.setdefault(root, set()).add(mode)
    if not roots:
        parser.error('Provide at least one directory to collect or inventory.')

    files = {}
    issues = []
    for root, modes in roots.items():
        def scan_error(error):
            issues.append({'kind': 'scan_error', 'message': str(error)})

        for base, dirs, names in os.walk(root, followlinks=False, onerror=scan_error):
            retained = []
            for name in sorted(dirs):
                child = Path(base) / name
                if name in SKIP_DIRS:
                    continue
                if child.is_symlink():
                    issues.append({'path': child.relative_to(source).as_posix(),
                                   'kind': 'symbolic_directory_skipped'})
                else:
                    retained.append(name)
            dirs[:] = retained
            for name in sorted(names):
                path = Path(base) / name
                files.setdefault(path, set()).update(modes)

    report = {'created_utc': datetime.now(timezone.utc).isoformat(),
              'copied_files': 0, 'copied_bytes': 0,
              'per_file_limit_bytes': FILE_LIMIT, 'total_limit_bytes': TOTAL_LIMIT,
              'issues': issues, 'complete_analysis_package': False}
    table = io.StringIO(newline='')
    writer = csv.writer(table, delimiter='\t', lineterminator='\n')
    writer.writerow(['relative_path', 'bytes', 'modified_utc', 'collection_status', 'sha256'])
    output.parent.mkdir(parents=True, exist_ok=True)
    created = False
    try:
        with ZipFile(output, mode='x', compression=ZIP_DEFLATED) as archive:
            created = True
            for path, modes in sorted(files.items()):
                relative = path.relative_to(source).as_posix()
                size, modified, digest = '', '', ''
                status = 'inventory_only'
                try:
                    info = path.lstat()
                    size = info.st_size
                    modified = datetime.fromtimestamp(info.st_mtime, timezone.utc).isoformat()
                    eligible = path.suffix.lower() in CODE or (
                        bool(modes & {'scripts', 'metadata'}) and
                        path.suffix.lower() in SMALL_DATA)
                    if not stat.S_ISREG(info.st_mode):
                        status = 'non_regular_or_symlink_skipped'
                    elif eligible:
                        if size > FILE_LIMIT:
                            status = 'file_size_limit'
                        elif report['copied_bytes'] + size > TOTAL_LIMIT:
                            status = 'total_size_limit'
                        else:
                            with path.open('rb') as handle:
                                data = handle.read(FILE_LIMIT + 1)
                            if len(data) > FILE_LIMIT or len(data) != size:
                                status = 'file_changed_or_size_limit'
                            elif report['copied_bytes'] + len(data) > TOTAL_LIMIT:
                                status = 'total_size_limit'
                            else:
                                digest = hashlib.sha256(data).hexdigest()
                                archive.writestr('source/' + relative, data)
                                report['copied_files'] += 1
                                report['copied_bytes'] += len(data)
                                status = 'copied'
                except OSError as error:
                    status = 'read_error'
                    issues.append({'path': relative, 'kind': status, 'message': str(error)})
                if status not in {'copied', 'inventory_only', 'read_error'}:
                    issues.append({'path': relative, 'kind': status})
                writer.writerow([relative, size, modified, status, digest])
            report['inventoried_files'] = len(files)
            archive.writestr('file_inventory.tsv', table.getvalue())
            archive.writestr('collection_report.json', json.dumps(report, indent=2) + '\n')
            archive.writestr('README.txt',
                'Review bundle: scripts, small metadata, and an inventory of selected directories.\n'
                'Most result and input data files are inventoried, not copied.\n'
                'All script versions are candidates; no final version is selected automatically.\n'
                'See collection_report.json for omissions and file_inventory.tsv for provenance.\n'
                'Review the contents before sharing. This tool performs no upload or Git operation.\n')
    except BaseException:
        if created:
            output.unlink()
        raise
    print('Created: ' + str(output))
    print('Copied files: {copied_files}; inventoried files: {inventoried_files}'.format(**report))
    print('Omissions or scan issues: ' + str(len(issues)))
    print('ZIP size: {:.2f} MiB'.format(output.stat().st_size / (1024 * 1024)))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as error:
        print('Collection failed: ' + str(error), file=sys.stderr)
        sys.exit(1)
