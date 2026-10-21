# Chapter 2 analysis scripts

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
| [02_MAG_catalogue](scripts/02_MAG_catalogue/) | 2 |
| [03_MAG_taxonomy](scripts/03_MAG_taxonomy/) | 1 |
| [04_KO_profiles_tests](scripts/04_KO_profiles_tests/) | 14 |
| [05_osmotic_strategy](scripts/05_osmotic_strategy/) | 12 |
| [06_redox_strategy](scripts/06_redox_strategy/) | 11 |
| [07_nitrogen_strategy](scripts/07_nitrogen_strategy/) | 3 |
| [08_carbon_phosphorus_strategy](scripts/08_carbon_phosphorus_strategy/) | 3 |
| [09_iron_strategy](scripts/09_iron_strategy/) | 2 |
| [10_colonisation_strategy](scripts/10_colonisation_strategy/) | 2 |
| [11_community_strategies](scripts/11_community_strategies/) | 10 |
| [12_phylogenetic_signal](scripts/12_phylogenetic_signal/) | 2 |
| [13_phylogeny_display](scripts/13_phylogeny_display/) | 8 |
| [14_redundancy_null_models](scripts/14_redundancy_null_models/) | 12 |
| [15_carriers_turnover](scripts/15_carriers_turnover/) | 9 |
| [16_lineage_decomposition](scripts/16_lineage_decomposition/) | 1 |
| [17_MAG_sharing](scripts/17_MAG_sharing/) | 2 |
| [18_niche_breadth](scripts/18_niche_breadth/) | 4 |
| [19_zeta_core](scripts/19_zeta_core/) | 5 |
| [20_16S_assembly](scripts/20_16S_assembly/) | 27 |
| [22_quality_validation](scripts/22_quality_validation/) | 1 |
| [23_CAZyme_exploratory](scripts/23_CAZyme_exploratory/) | 2 |
| [90_preprocessing_reference](scripts/90_preprocessing_reference/) | 49 |

## Coverage and validation

All supplied script versions are retained, including pilot and backup versions. Generating scripts for the soil metadata summaries and the matched 16S–MAG turnover comparison were not available in the reviewed collection.

Python and shell/PBS syntax was checked during publication preparation. R syntax was not independently checked by this helper. Input availability, package compatibility and scientific reproducibility still require verification before execution.

The full local collection remains separate and is excluded from Git.
