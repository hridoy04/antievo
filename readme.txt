# antievo

> Predict epistasis‑aware antigenicity from a single protein sequence and map site‑wise signals onto 3D structure.

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](#requirements) [![License](https://img.shields.io/badge/license-MIT-informational.svg)](./LICENSE) [![PRs welcome](https://img.shields.io/badge/PRs-welcome-success.svg)](#contributing)

---

## Overview

This repo scores each residue in a query sequence for **antigenicity** (with **epistatic context** and **surface accessibility**) relative to a reference, then writes:

* a **per‑site table** of antigenicity / change likelihood vs the reference,
* **colored PDBs** where residue color encodes antigenicity or Δ‑likelihood,
* optional plots or summaries (WIP).

Core scripts (from this repo):

* `EpistaticSurfaceAccessibilityVirus_V4.py` – main scoring pipeline (latest).
* `EpistaticSurfaceAccessibilityVirus_V2.py` – earlier version.
* `generate_colored_pdbs.py` – map per‑site scores to a structure (PDB/mmCIF) and write colored models.

> Input: **one sequence** (+ a **reference** sequence), optional **PDB** for structural mapping.

---

## Architecture

```mermaid
flowchart LR
    subgraph Input[Inputs]
        QFA[Query sequence (FASTA / string)]
        RFA[Reference sequence]
        PDB[(Structure PDB/mmCIF)]
    end

    subgraph Feats[Feature/Context]
        EPI[Epistasis model\n(pairwise / higher‑order)]
        ACC[Surface accessibility\n(SASA / exposure)]
        ALIGN[Alignment & site mapping]
    end

    subgraph Score[Scoring]
        DIFF[Per‑site Δ vs ref]
        AGG[Combine epistasis + accessibility →\nantigenicity score]
    end

    subgraph Out[Outputs]
        TAB[CSV/TSV per‑site table]
        COLPDB[Color residues by score\n→ PDB/MMCIF]
        PLOT[Summary plots (WIP)]
    end

    QFA --> ALIGN
    RFA --> ALIGN
    ALIGN --> EPI --> AGG
    ALIGN --> ACC --> AGG
    AGG --> DIFF --> TAB
    AGG --> COLPDB
    DIFF --> COLPDB
    TAB --> PLOT
    PDB --- COLPDB
```

---

## Project structure

```
.
├── EpistaticSurfaceAccessibilityVirus_V4.py   # main pipeline
├── EpistaticSurfaceAccessibilityVirus_V2.py   # legacy
├── generate_colored_pdbs.py                   # structure coloring
├── readme.ipynb                               # draft notes
└── README.md                                  # this file
```

---

## Quickstart

### 1) Environment

We recommend a fresh virtualenv or conda env.

```bash
# conda (example)
conda create -n antievo python=3.10 -y
conda activate antievo

# or venv
python -m venv .venv && source .venv/bin/activate

# install deps (edit as needed)
pip install -r requirements.txt  # if present
# otherwise install common scientific stack used here
pip install numpy pandas biopython matplotlib
```

> If you plan to color structures, install a PDB writer/viewer stack as needed. The included script uses standard Python + Biopython I/O; PyMOL/ChimeraX are optional for viewing.

### 2) Run scoring (sequence → per‑site table)

```bash
python EpistaticSurfaceAccessibilityVirus_V4.py \
  --query path/to/query.fasta \
  --reference path/to/reference.fasta \
  --out results/ \
  --id QUERY1 \
  --format csv
```

**Expected output** (under `results/`):

* `QUERY1_sites.csv` with columns like: `site,index,ref_aa,query_aa,delta,epistasis,accessibility,antigenicity`

> Flags are indicative. If your script exposes different names, adjust accordingly (see `--help`).

### 3) Map scores to structure (per‑site table → colored PDB)

```bash
python generate_colored_pdbs.py \
  --pdb path/to/structure.pdb \
  --scores results/QUERY1_sites.csv \
  --score-column antigenicity \
  --out colored_pdbs/
```

This writes a PDB (or mmCIF) where each residue’s B‑factor (or alt field) is repurposed for the chosen score and/or a color palette is baked into the output. Load the file in **PyMOL** or **ChimeraX** and color by B‑factor to visualize hot/cold antigenic sites.

---

## CLI (proposed)

If your scripts don’t yet have CLIs, consider adding these canonical flags:

* `EpistaticSurfaceAccessibilityVirus_V4.py`

  * `--query <fasta|seq>`: query sequence (FASTA file or raw AA string)
  * `--reference <fasta|seq>`: reference sequence
  * `--model <path|name>`: epistasis model spec (if external)
  * `--exposure <method>`: surface accessibility method (e.g., ASA from structure or predicted)
  * `--out <dir>`: output directory
  * `--format {csv,tsv}`: table format

* `generate_colored_pdbs.py`

  * `--pdb <file>`: input structure (chain can be specified via `--chain`)
  * `--scores <csv>`: per‑site table
  * `--score-column <name>`: which column to map
  * `--palette <name>`: e.g., `RdBu`, `viridis`
  * `--out <dir>`: output directory

---

## Data alignment & indexing

Residue indices must match between the **per‑site table** and the **structure** (chain/offset). If your PDB numbering differs, use `--offset`/`--chain` options or pre‑align sequences to map indices correctly. Gaps/insertions should be handled explicitly.

---

## Examples

```bash
# score a single HA sequence against a reference and write a colored model
python EpistaticSurfaceAccessibilityVirus_V4.py \
  --query data/HA_query.fasta \
  --reference data/HA_ref.fasta \
  --out runs/ha_query

python generate_colored_pdbs.py \
  --pdb data/HA_structure.pdb --chain A \
  --scores runs/ha_query/QUERY_sites.csv \
  --score-column antigenicity \
  --out runs/ha_query/structure
```

---

## Requirements

* Python 3.10+
* `numpy`, `pandas`, `biopython` (I/O & sequence/structure utilities)
* Optional: `matplotlib` (plots), `scipy`

> Create a `requirements.txt` reflecting your exact imports once finalized.

---

## Testing

* Add unit tests for: alignment → index map, epistasis aggregation, accessibility normalization, and PDB coloring.
* Edge cases: non‑standard residues, gaps/insertions, missing chains, mixed numbering.

Example skeleton (pytest):

```python
# tests/test_scoring.py
import pandas as pd
from pathlib import Path

# def test_epistasis_accessibility_integration():
#     ...
```

---

## Contributing

Issues and PRs are welcome. Please include minimal repro sequences and a small PDB snippet when reporting structural mapping bugs.

---

## License

MIT (or update to your preferred license).

---

## Citation

If this code supports a paper or preprint, add a BibTeX entry here.

---

## 👀 Quick Visuals (Flowcharts, Demo Figs & Overview)

> Download the single‑image overview and place it at `docs/figs/overview.png`.
>
> ![Overview](docs/figs/overview.png)

---

### A. One‑Glance Flow

### A. One‑Glance Flow

```mermaid
graph TD
    A[User provides Query Sequence
(+ optional PDB) ] --> B[Align to Reference]
    B --> C[Epistasis Context
(pairwise / higher‑order)]
    B --> D[Surface Accessibility
(exposure / SASA)]
    C --> E[Combine Signals →
Per‑site Antigenicity]
    D --> E
    E --> F[Δ vs Reference
(change likelihood)]
    E --> G[Per‑site Table (CSV/TSV)]
    F --> G
    G --> H[Color Structure (PDB/mmCIF)]
    H --> I[Visualize in PyMOL/ChimeraX]
```

### B. What the Outputs Mean

```mermaid
flowchart LR
    subgraph Table[Per‑site Table]
        T1[site/index]
        T2[ref_aa → query_aa]
        T3[epistasis_score]
        T4[accessibility]
        T5[antigenicity]
        T6[delta_vs_ref]
    end
    subgraph Structure[Colored Structure]
        S1[Residue color encodes
antigenicity]
        S2[B‑factor stores value]
    end
    Table -->|map by index/chain| Structure
```

> **Legend**
>
> * Higher color intensity ⇒ higher predicted antigenicity.
> * Δ vs ref > 0 ⇒ more likely impactful change relative to reference at that site.
> * Ensure chain/index mapping is correct (use `--chain`/`--offset`).

### C. End‑to‑End (Expanded) — For New Readers

```mermaid
sequenceDiagram
    participant U as You
    participant Q as Query FASTA
    participant R as Reference FASTA
    participant S as Structure (PDB)
    participant AE as antievo
    participant O as Outputs

    U->>AE: run EpistaticSurfaceAccessibilityVirus_V4.py
    Q-->>AE: sequence
    R-->>AE: reference
    AE->>AE: align & index map
    AE->>AE: compute epistasis context
    AE->>AE: add surface accessibility
    AE->>O: write per‑site table
    U->>AE: run generate_colored_pdbs.py
    S-->>AE: PDB/mmCIF
    O-->>AE: table (antigenicity/Δ)
    AE->>O: colored structure
    U->>O: open in PyMOL/ChimeraX
```

### D. Demo Figures (placeholders)

> Put the following images in `docs/figs/` and they will render here when committed. Replace with real screenshots from your outputs.

| Figure                                                | Description                                                                               |
| ----------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| ![Per‑site table preview](docs/figs/demo_table.png)   | A few rows highlighting `site`, `ref_aa`, `query_aa`, `antigenicity`, and `delta_vs_ref`. |
| ![Colored structure](docs/figs/colored_structure.png) | Ribbon view colored by antigenicity (blue→white→red).                                     |
| ![Δ vs ref barplot](docs/figs/delta_barplot.png)      | Top‑N sites with largest positive Δ.                                                      |

**Tip (PyMOL coloring):**

```bash
# after loading colored PDB, color by b‑factor in PyMOL
spectrum b, blue_white_red, minimum=MIN, maximum=MAX
as cartoon
cartoon putty
```

**Tip (ChimeraX coloring):**

```bash
# use attribute coloring if antigenicity is stored in B‑factor
color byattribute bfactor palette blueswhiteorange reverse false
cartoon
```
