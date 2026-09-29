# Egyptian H5Nx and H9N2 haemagglutinin antigenic-site analysis, 2010–2026

Analysis code and derived results for:

> Khalil AA. *Haemagglutinin antigenic-site conservation within clades and a single clade
> replacement in Egyptian H5Nx and H9N2 avian influenza viruses: implications for sequence-based
> surveillance of vaccine seed match, 2010–2026.*

**Data cut-off: 28 September 2026.** Both source repositories grow continuously, so a retrieval run
later will return more records than are reported in the paper. The cut-off is bounded on record
publication date for GenBank, which makes that component exactly reproducible.

---

## What is and is not in this repository

| | |
|---|---|
| **Included** | All analysis and figure code; the column-to-H3 numbering maps; every derived result as JSON; accession lists with source, clade and geographic assignment; the published figures |
| **Not included** | Raw sequence data |

**Sequences are deliberately absent.** GISAID data may not be redistributed under the EpiFlu terms
of use. The GISAID component is instead identified by its citable EPI_SET:

- H5: `EPI_SET_260928dp` — https://doi.org/10.55876/gis8.260928dp (1,442 viruses)
- H9: `EPI_SET_260928rz` — https://doi.org/10.55876/gis8.260928rz (780 viruses)

The GenBank component is redistributable but is reconstructed exactly by running the retrieval
scripts, so it is not duplicated here. `data/SuppTableS1_accessions.tsv` lists every accession used,
its repository of origin, whether it was held by both, and its clade assignment.

---

## Reproducing the analysis

### 1. Retrieve

GenBank, on any machine with network access to NCBI E-utilities:

```powershell
cd scripts/retrieval
.\fetch_genbank_egypt.ps1            # Egyptian H5 and H9, HA and NA
.\fetch_genbank_global_context.ps1   # non-Egyptian H5 HA for the contextual tree
```

Both bound the query on publication date (`$CUTOFF`), so re-running after the cut-off returns the
same record set. Edit `$email` to your own address before running; NCBI requires a contact.

GISAID, in the EpiFlu web interface (no API is available for EpiFlu):

- Type **A**; H **5** (then repeat with **9**); **N left blank**; no clade filter
- Location **Africa / Egypt**; collection date **2010-01-01 to 2026-09-28**
- Required segment **HA**
- Download **Isolates as XLS** and **Sequences (DNA) as FASTA**

Leaving the neuraminidase unrestricted is essential, not incidental. An H5 query pinned to N1
returns no H5N8 at all and therefore cannot observe clade 2.3.4.4b during the years it established
in Egypt — see Methods §2.1.

### 2. Analyse

Run in order from the repository root:

```bash
python scripts/analysis/01_consolidate.py          # merge, geographic filter, strain-level dedup
python scripts/analysis/02_translate_qc.py         # 3-frame translation, ORF selection, QC
mafft --auto out/H5_protein.faa > out/H5_aln.faa   # and likewise for H9
python scripts/analysis/03_h3_numbering_map.py     # column -> H3 map by alignment to PDB 2FK0
python scripts/analysis/04_master_analysis.py      # conservation, substitutions, RBS, entropy
python scripts/analysis/05_entropy_rarefaction.py  # rarefied entropy + breakpoint scan
python scripts/analysis/06_random_forest.py        # RF vs persistence and random baselines
python scripts/analysis/07_na_phylogeny.py         # N1 and N2 trees
python scripts/analysis/08_global_context_tree.py  # Egyptian 2.3.4.4b in global context
```

Clade assignment uses Nextclade v3.23.0 against the community `moncla-lab/iav-h5` all-clades
dataset (2026-04-14 build):

```bash
nextclade run --input-dataset <dataset> --output-tsv out/nextclade_H5.tsv out/H5_consolidated.fasta
```

### 3. Figures

```bash
python scripts/figures/fig1_4_6.py        # Figs 1, 2, 3, 4, 6
python scripts/figures/fig5_phylogeny.py  # Fig 5
python scripts/figures/fig7_structure.py  # Fig 7
python scripts/figures/fig8_na_trees.py   # Fig 8
```

---

## Three methodological points that changed results

These are documented here because each materially altered a finding, and each is a trap that is easy
to fall into with public influenza sequence data.

**1. Retrieve without restricting the neuraminidase subtype.** Our first retrieval was inadvertently
pinned to H5N1 and returned no H5N8. Clade 2.3.4.4b was carried almost entirely on N8 in Egypt from
2016 to 2020, so that retrieval placed the clade transition in 2021 rather than 2016 — five years
late.

**2. Deduplicate by virus strain, not by accession.** 87.5% of the GenBank H5 strain set and 87.0%
of H9 are also in GISAID. Deduplicating by accession enters the same virus twice and biases
residue-frequency and entropy estimates toward double-deposited viruses.

**3. Assign H3-equivalent positions by alignment, not by a scalar offset.** A single offset
calibrated at the receptor-binding site is wrong elsewhere, because the H5/H3 offset changes across
indel regions. Under a scalar offset, positions 50, 58 and 128 disagreed with the deposited
structure while 226 and 228 agreed by construction. `03_h3_numbering_map.py` instead aligns the
consensus to PDB 2FK0 chain A — an H5 sequence deposited in H3 numbering — and validates at all
fifteen positions checked.

A fourth, smaller trap: the Entrez term `"NA"[All Fields]` matches free text as well as the gene
symbol, so a segment 6 query also returns segment 4 records. `07_na_phylogeny.py` filters on the
structured segment qualifier and on the record title.

---

## Layout

```
scripts/retrieval/   PowerShell retrieval from NCBI E-utilities
scripts/analysis/    Consolidation through to the phylogenetic tests, numbered in run order
scripts/figures/     Figure generation
data/                Accession lists and per-strain metadata (no sequences)
results/             Every derived number as JSON, including all values quoted in the paper
figures/             Published figures at 400 dpi
reference/           PDB 2FK0 chain A Cα coordinates in H3 numbering
```

`results/NUMBERS.json` is the single source of truth for the figures quoted in the manuscript.

---

## Requirements

Python ≥ 3.11 with `biopython`, `numpy`, `scipy`, `pandas`, `scikit-learn`, `matplotlib`, `xlrd`;
plus MAFFT v7.505, FastTree 2 and Nextclade v3.23.0 on the path. See `requirements.txt`.

---

## Licence and citation

Code is released under the MIT licence (`LICENSE`). If you use it, please cite the paper and this
repository; `CITATION.cff` carries the metadata.

Sequence data remain subject to the terms of their source repositories. Users of the GISAID
component must acknowledge the originating and submitting laboratories; the EPI_SET identifiers
above resolve to the full contributor list.
