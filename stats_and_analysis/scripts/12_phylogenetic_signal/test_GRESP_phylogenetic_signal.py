#!/usr/bin/env python3

import csv
import math
import random
import re
from pathlib import Path

BASE = Path("results/ecological_strategies/phylogenetic_signal")
SCORES = BASE / "GRESP_six_strategies_352_MAGs.tsv"
TREE = Path("results/phylogeny_bacteria_352/bacteria_352_bac120.tree")
NPERM = 1999
SEED = 20261009

STRATEGIES = [
    "Osmotic_adaptation",
    "Redox_flexibility",
    "Nitrogen_acquisition_transformation",
    "Carbon_phosphorus_resource",
    "Iron_strategy_score",
    "Colonisation_host_strategy_score",
]

# Read Newick tree using Python standard library.
class Node:
    def __init__(self, name="", children=None, length=0.0):
        self.name = name.strip("'")
        self.children = children or []
        self.length = length


def read_tree(path):
    tokens = re.findall(
        r"'(?:[^']|'')*'|[(),:;]|[^(),:;\s]+",
        path.read_text()
    )
    pos = 0

    def parse():
        nonlocal pos
        children = []
        name = ""

        if tokens[pos] == "(":
            pos += 1
            while True:
                children.append(parse())
                if tokens[pos] == ",":
                    pos += 1
                elif tokens[pos] == ")":
                    pos += 1
                    break
                else:
                    raise ValueError("Invalid Newick structure")

            if tokens[pos] not in (":", ",", ")", ";"):
                name = tokens[pos]
                pos += 1
        else:
            name = tokens[pos]
            pos += 1

        length = 0.0
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 1
            length = float(tokens[pos])
            pos += 1

        return Node(name, children, length)

    root = parse()
    if pos >= len(tokens) or tokens[pos] != ";":
        raise ValueError("Tree must end with a semicolon")
    return root


root = read_tree(TREE)
postorder = []
tips = []

def collect(node):
    for child in node.children:
        collect(child)
    postorder.append(node)
    if not node.children:
        tips.append(node.name)

collect(root)

if len(tips) != 352 or len(set(tips)) != 352:
    raise ValueError("Tree tip count or uniqueness check failed")

# Check branch lengths; regularize exact zeros.
lengths = [n.length for n in postorder if n is not root]
if any(x < 0 or not math.isfinite(x) for x in lengths):
    raise ValueError("Invalid or negative branch length")

positive = [x for x in lengths if x > 0]
if not positive:
    raise ValueError("Tree has no positive branch lengths")

epsilon = max(min(positive) * 1e-6, 1e-12)
zero_count = 0

for node in postorder:
    if node is not root and node.length == 0:
        node.length = epsilon
        zero_count += 1

# Root-to-tip distances for Brownian covariance trace.
depths = {}

def get_depths(node, depth=0.0):
    if not node.children:
        depths[node.name] = depth
    for child in node.children:
        get_depths(child, depth + child.length)

get_depths(root)
n = len(tips)
trace_C = sum(depths.values())
tip_index = {name: i for i, name in enumerate(tips)}

# Independent contrasts calculate the GLS quantities
# needed for Blomberg's K without a matrix inversion.
def contrasts(values):
    states = {}
    for node in postorder:
        if not node.children:
            states[node] = (
                values[tip_index[node.name]],
                node.length,
                0.0
            )
            continue

        first = node.children[0]
        mean, var, ss = states[first]

        for child in node.children[1:]:
            m2, v2, ss2 = states[child]
            total = var + v2
            if total <= 0:
                raise ValueError("Invalid contrast variance")
            ss += ss2 + (mean - m2) ** 2 / total
            mean = (mean * v2 + m2 * var) / total
            var = var * v2 / total

        if node is not root:
            var += node.length

        states[node] = (mean, var, ss)

    return states[root]


_, root_var, _ = contrasts([0.0] * n)
expected_ratio = (trace_C - n * root_var) / (n - 1)

if expected_ratio <= 0:
    raise ValueError("Invalid Brownian normalization")

def blomberg_k(values):
    mean, _, ss = contrasts(values)
    if ss <= 0:
        return float("nan")
    ordinary_ss = sum((x - mean) ** 2 for x in values)
    return (ordinary_ss / ss) / expected_ratio


# Load and validate six GRESP scores.
with SCORES.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

data = {r["Final_MAG"]: r for r in rows}

if len(rows) != 352 or set(data) != set(tips):
    raise ValueError("MAG identifiers do not match the tree")

rng = random.Random(SEED)
results = []

for strategy in STRATEGIES:
    values = [float(data[t][strategy]) for t in tips]

    if not all(math.isfinite(x) and 0 <= x <= 1 for x in values):
        raise ValueError(f"Invalid scores: {strategy}")

    if len(set(values)) < 2:
        results.append([strategy, None, None, None])
        print(strategy, "NO VARIATION", flush=True)
        continue

    observed = blomberg_k(values)
    exceed = 0

    for _ in range(NPERM):
        shuffled = values.copy()
        rng.shuffle(shuffled)
        k_null = blomberg_k(shuffled)
        if k_null >= observed:
            exceed += 1

    p = (exceed + 1) / (NPERM + 1)
    results.append([strategy, observed, p, None])
    print(strategy, "K =", round(observed, 4),
          "p =", round(p, 5), flush=True)


# Benjamini-Hochberg FDR correction.
valid = [i for i, r in enumerate(results) if r[2] is not None]
ordered = sorted(valid, key=lambda i: results[i][2])
m = len(ordered)
previous = 1.0

for rank in range(m, 0, -1):
    i = ordered[rank - 1]
    previous = min(previous, results[i][2] * m / rank)
    results[i][3] = previous

output = BASE / "GRESP_phylogenetic_signal_K.tsv"

with output.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "Strategy", "Blombergs_K",
        "Permutation_p", "FDR_q"
    ])
    for row in results:
        writer.writerow(row)

print("\nTree tips:", n)
print("Root descendants:", len(root.children))
print("Zero-length branches regularized:", zero_count)
if len(root.children) != 2:
    print("WARNING: Check tree rooting before interpretation")
print("Permutations:", NPERM)
print("Results saved:", output)
