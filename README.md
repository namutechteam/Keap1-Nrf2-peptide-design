# Carrying a Generative Peptide Design Workflow to Measurement at the Keap1–Nrf2 Interface

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-3776AB.svg)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/data-10.6084%2Fm9.figshare.33394336-C43B69.svg)](https://doi.org/10.6084/m9.figshare.33394336)

Hye Ree Yoon<sup>†</sup>, Kyoung Pyo Kwon<sup>†</sup>, Jin Hee Lee<sup>\*</sup> and Min Sun Yeom<sup>\*</sup>

NamuICT R&D Center, NamuICT, 41 Magok Jungang 8-ro, Seoul, 07793, Republic of Korea

<sup>†</sup> These authors contributed equally. <sup>\*</sup> Corresponding authors.

---

Code and data for the prospective design and structure-based prioritization of peptides
targeting the Keap1–Nrf2 protein–protein interaction. Manuscript submitted to the
Journal of Chemical Information and Modeling.

## Pipeline

```
01  ProteinMPNN          backbone-conditioned sequence generation          .fa
02  sequence filtering   seven position and motif filters                  .txt
03  Chai-1               complex prediction, five models per sequence      .cif
04  structure filtering  intrapeptide H-bond + PLIP hot-spot contacts      .pdb + summary
05  HADDOCK3             information-driven docking and scoring            capri_ss.tsv
```

## Contents

| | |
|---|---|
| [`workflow/`](workflow/) | the five pipeline stages, each with `scripts/` and, where needed, `configs/` |
| [`example/`](example/) | the 12-residue run with the input for every stage, so any stage can be run on its own |
| [`analysis/data/`](analysis/data/) | the 45 measured peptides, the sequence set at each retained checkpoint, and Chai-1 confidence values |
| [`analysis/code/`](analysis/code/) | scripts producing the manuscript tables and figures from that data |
| [`analysis/figures/`](analysis/figures/) | the cascade-composition figure and the table-of-contents graphic |
| [`environment/`](environment/) | conda environment specifications and `setup.sh` |

Parameters for each stage are in its own README:
[01](workflow/01_proteinmpnn/README.md) ·
[02](workflow/02_sequence_filtering/README.md) ·
[03](workflow/03_chai1/README.md) ·
[04](workflow/04_chai_filtering/README.md) ·
[05](workflow/05_haddock3/README.md)

## Installation

```bash
git clone https://github.com/namutechteam/Keap1-Nrf2-peptide-design.git
cd Keap1-Nrf2-peptide-design
bash environment/setup.sh
```

Two environments are created: `keap1nrf2` for stages 01–04 and the analysis scripts,
`keap1nrf2-haddock` for stage 05. Two are needed because `chai_lab` pins numpy 1.x and
`haddock3` pins numpy 2.x. ProteinMPNN is cloned separately; point `MPNN_DIR` at it.

## Usage

Commands, stage by stage, are in [example/](example/).

## Data availability

The full structural output reported in the manuscript is deposited at Figshare
(DOI: [10.6084/m9.figshare.33394336](https://doi.org/10.6084/m9.figshare.33394336)).

## Third-party software

| Software | Version | License |
|---|---|---|
| [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) | v1.0.1, weights `v_48_020` | MIT |
| [Chai-1](https://github.com/chaidiscovery/chai-lab) (`chai_lab`) | 0.5.2 | Apache 2.0 |
| [PLIP](https://github.com/pharmai/plip) | 2.4.0 | GPL-2.0 |
| [HADDOCK3](https://github.com/haddocking/haddock3) | 2024.10.0b7 | Apache 2.0 |

## License

[MIT](LICENSE)
