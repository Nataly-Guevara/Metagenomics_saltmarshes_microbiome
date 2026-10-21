#!/usr/bin/env python3

from pathlib import Path
import re

tree_file = Path(
    "/path/to/metagenomics_project/"
    "results/phylogeny_bacteria_352/bacteria_352_bac120.tree"
)

output_file = Path(
    "/path/to/metagenomics_project/"
    "results/phylogeny_bacteria_352/tree_tips_352.txt"
)

newick = tree_file.read_text().strip()

# Extract labels that occur immediately after "(" or ","
pattern = re.compile(r"(?<=[(,])([^():,;\s]+)(?=[:),;])")
tips = pattern.findall(newick)

# Remove accidental duplicates while preserving order
seen = set()
unique_tips = []

for tip in tips:
    if tip not in seen:
        seen.add(tip)
        unique_tips.append(tip)

output_file.write_text("".join(f"{tip}\n" for tip in unique_tips))

print(f"Tree file: {tree_file}")
print(f"Tips detected: {len(unique_tips)}")
print(f"Output: {output_file}")

print("\nFirst 10 tips:")
for tip in unique_tips[:10]:
    print(tip)
