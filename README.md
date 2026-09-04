# HEG1 Binder Design Pipeline

A Snakemake pipeline for de novo miniprotein binder design against **HEG1**
(UniProt [Q13585](https://www.uniprot.org/uniprotkb/Q13585)), an endothelial
transmembrane adhesion molecule. This is a research-tool effort to design
miniprotein binders against defined epitopes on HEG1's extracellular domain,
in support of studying whether that domain is required for mechanical flow
sensing.

## Pipeline overview

RFdiffusion → make_resfile → Rosetta fixbb → renumber → AF2 initial guess → aggregate


| Stage | Tool | What it does |
|---|---|---|
| `rfdiffusion` | [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) | Generates a binder backbone against a specified epitope (contig map + hotspot residues) |
| `make_resfile` | Custom (BioPython) | Builds a Rosetta resfile: `ALLAA` for the designed binder chain, `NATRO` for the fixed target chain |
| `rosetta_fixbb` | [Rosetta](https://www.rosettacommons.org/) `fixbb` | Designs a sequence onto the RFdiffusion backbone |
| `renumber` | Custom (BioPython) | Renumbers residues (binder chain first, then target) for downstream compatibility |
| `af2_initial_guess` | [dl_binder_design](https://github.com/nrbennet/dl_binder_design) (Baker Lab) | Validates the design by predicting the complex structure with AlphaFold2, using the designed structure as an initial guess |
| `aggregate` | Custom | Collects all per-design AF2 scores into a single ranked results table |

Designs are generated per-epitope, with multiple Rosetta sequence variants
(`nstruct`) per RFdiffusion backbone, and scored by AF2-predicted metrics
(`pae_interaction`, `plddt_binder`, `plddt_target`, RMSDs, etc.).

## Repository contents

.
├── Snakefile                     # Pipeline orchestration (Snakemake)
├── config.yaml                   # Epitope definitions, hotspots, tool paths, Rosetta/AF2 params
├── scripts/
│   ├── make_resfile.py           # Builds Rosetta resfile from an RFdiffusion backbone
│   ├── renumber.py               # Renumbers binder/target chains post-Rosetta
│   └── aggregate_results.py      # Aggregates per-design AF2 score files into one ranked table
└── README.md


This repository contains **pipeline code only** — no design outputs,
structures, logs, or intermediate results are included.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/)
- [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) (conda env, e.g. `SE3nv`)
- [Rosetta](https://www.rosettacommons.org/software/license-and-download) (`fixbb` binary)
- [dl_binder_design](https://github.com/nrbennet/dl_binder_design) (Baker Lab AF2 initial guess, conda env)
- BioPython (conda env, used by `make_resfile.py` and `renumber.py`)
- A target structure (e.g. an AlphaFold model or truncated region of the target protein)

## Configuration (`config.yaml`)

```yaml
target_pdb: /path/to/target_structure.pdb
contigmap: "[A<start>-<end>/0 <min_len>-<max_len>]"

epitopes:
  ep_name:
    hotspot_res: "A<res1>,A<res2>,..."
    n_designs: <int>

rosetta:
  nstruct: <int>
  ex_flags: "-ex1 -ex2 -use_input_sc"
  binary: /path/to/fixbb.static.linuxgccrelease

af2:
  recycle: <int>
  script: /path/to/dl_binder_design/af2_initial_guess/predict.py

envs:
  rfdiffusion: <conda_env_name>
  biopython: <conda_env_name>
  af2: <conda_env_name>

parallelism:
  rosetta_jobs: <int>
  af2_jobs: <int>
```

- **`target_pdb`** — path to the target structure to design against
- **`contigmap`** — RFdiffusion contig map defining the target region and binder length range
- **`epitopes`** — one entry per epitope: hotspot residues (as RFdiffusion `ppi.hotspot_res`) and how many designs to generate against it
- **`rosetta`** — number of sequence variants per backbone (`nstruct`), extra rotamer flags, and path to the `fixbb` binary
- **`af2`** — number of AF2 recycles and path to the `predict.py` script from `dl_binder_design`
- **`envs`** — conda environment names for each tool
- **`parallelism`** — job counts for Snakemake's `--jobs` / cluster or local parallel execution

## Usage

```bash
# Dry run to check the DAG
snakemake -n

# Run locally with N parallel jobs
snakemake --cores <N>

# Or, for GPU-bound RFdiffusion steps, tune parallelism per-rule as needed
```

Output structure:

results/
├── 00_rfdiffusion/ # RFdiffusion backbone outputs
├── 01_rosetta/ # Rosetta fixbb sequence designs + resfiles
├── 02_renumbered/ # Renumbered designs (binder chain first)
├── 03_af2/ # AF2 initial guess predictions + per-design scores
└── final_results_table.csv # Aggregated, ranked results across all designs


## Scoring

Each design is scored via AF2 initial guess on the following metrics (among others):

- **`pae_interaction`** — predicted aligned error between binder and target chains (lower is better; primary ranking metric used here)
- **`plddt_binder`** / **`plddt_target`** — per-chain confidence
- **`binder_aligned_rmsd`** / **`target_aligned_rmsd`** — structural deviation from the design model

**Note:** `pae_interaction` alone can be misleading for designs where the
target region has low confidence (`plddt_target`) — a confidently-folded
binder next to a poorly-resolved target can still produce a deceptively good
`pae_interaction` average. Cross-checking with
[ipSAE](https://github.com/DunbrackLab/IPSAE) (Dunbrack et al., 2025), which
restricts scoring to residue pairs with good interchain PAE, is recommended
before treating a design as a confident hit.

## Known limitations

- `rosetta_fixbb` and `af2_initial_guess` rules currently run one design at a
  time per Snakemake job; the `parallelism` config values are intended for
  Snakemake's own job scheduling (`--jobs`), not intra-rule batching.
- `renumber.py` assumes chain `B` is the designed binder and chain `A` is the
  target — confirm this matches your RFdiffusion output convention if
  adapting this pipeline to a different target/setup.

## License

Add a license of your choice (e.g. MIT) if you intend this repository to be
reused by others.
