# Chapter 2: statistics and analyses

This folder is being prepared for the downstream statistical and ecological analyses of the saltmarsh metagenomes. Analysis scripts, inputs, and outputs have not yet been uploaded.

## Organisation

Create a descriptive folder for each analysis, with the following contents:

| Path | Contents |
| --- | --- |
| `README.md` | Research question, method, input descriptions, software requirements, run command, and output descriptions. |
| `scripts/` | Scripts and small configuration files needed to run the analysis. |
| `inputs/` | Exact input tables used, or a documented reference to shared inputs. |
| `outputs/` | Selected tables, figures, and reports produced by the documented run. |

Place inputs used by several analyses in `shared_inputs/` once populated, and reference them explicitly in the relevant analysis READMEs.

## Reproducibility

For each analysis:

1. Identify the script and input versions used to generate the outputs.
2. Include all required input tables, or document how they are generated from upstream analyses.
3. Record software/package versions, parameters, and random seeds where relevant.
4. Use paths relative to the repository or clearly documented configurable data locations.
5. Describe validation and sensitivity analyses separately from the main analysis.

Large sequencing files, assemblies, alignments, and reference databases should be documented through suitable research storage or public archive references when required for upstream processing.

Check ignored files before committing: the repository's existing `.gitignore` excludes some analysis file types, including R objects and compressed files.
