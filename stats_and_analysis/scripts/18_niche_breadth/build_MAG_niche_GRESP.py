import csv
from pathlib import Path

base = Path("results/ecological_strategies")

sources = {
    "Osmotic": (
        "osmotic_adaptation/tables/osmotic_strategy_final_373_MAGs.tsv",
        "Osmotic_strategy_score"
    ),
    "Redox": (
        "redox_flexibility/tables/redox_strategy_final_v2_373_MAGs.tsv",
        "Redox_strategy_score"
    ),
    "Nitrogen": (
        "nitrogen_acquisition_transformation/tables/nitrogen_strategy_final_373_MAGs.tsv",
        "Nitrogen_strategy_score"
    ),
    "Carbon_P": (
        "carbon_phosphorus_resource/tables/carbon_phosphorus_strategy_final_373_MAGs.tsv",
        "Carbon_Phosphorus_strategy_score"
    ),
    "Iron": (
        "iron_acquisition/tables/iron_strategy_final_373_MAGs.tsv",
        "Iron_strategy_score"
    ),
    "Colonisation": (
        "colonisation_host_interaction/tables/colonisation_host_strategy_final_373_MAGs.tsv",
        "Colonisation_host_strategy_score"
    )
}

niche_file = base / "MAG_niche_breadth.tsv"

with niche_file.open() as f:
    niche = {
        r["MAG_ID"]: r
        for r in csv.DictReader(f, delimiter="\t")
    }

assert len(niche) == 373

for strategy, (filename, column) in sources.items():
    with (base / filename).open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    values = {r["Final_MAG"]: float(r[column]) for r in rows}

    assert len(values) == 373
    assert set(values) == set(niche)
    assert all(0 <= v <= 1 for v in values.values())

    for mag in niche:
        niche[mag][strategy] = values[mag]

output = base / "MAG_niche_breadth_GRESP_373.tsv"

fields = list(next(iter(niche.values())).keys())

with output.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(niche.values())

print("SUCCESS")
print("Matched MAGs:", len(niche))
print("GRESP dimensions:", len(sources))

print("\nStrategy score ranges:")
for strategy in sources:
    vals = [r[strategy] for r in niche.values()]
    print(
        f"{strategy:15s}",
        f"min={min(vals):.3f}",
        f"max={max(vals):.3f}"
    )

print("\nOutput:", output)
