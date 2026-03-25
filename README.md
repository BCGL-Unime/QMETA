# QMETA — QIIME2 Interactive Metagenomics Pipeline

An interactive, step-by-step Python pipeline for 16S/ITS amplicon metagenomics analysis built on top of [QIIME2](https://qiime2.org/). The pipeline guides the user through the entire workflow — from raw FASTQ import to differential abundance analysis — with interactive prompts at each step and the ability to resume from existing intermediate outputs.

## Features

- Supports both **paired-end** and **single-end** reads
- Optional primer/adapter trimming with **Cutadapt**
- Denoising and ASV generation with **DADA2**
- Taxonomic classification with a **pre-trained QIIME2 classifier** (e.g. SILVA 138)
- Automatic filtering of **Eukaryota** and **Unassigned** features
- Taxa collapse and export at levels **2 (phylum), 5 (family), 6 (genus)**
- **Alpha diversity** analysis (Chao1, Shannon, Simpson, ACE, Pielou's e, Simpson's e, observed features)
- **Core metrics** (non-phylogenetic)
- **Alpha rarefaction** curves
- **Differential abundance** analysis with ANCOM-BC
- Smart **resume support**: if intermediate outputs already exist, the pipeline asks whether to reuse them
- Full **log file** saved to `qiime2_out/pipeline.log`

---

## Requirements

### Software

| Tool | Notes |
|------|-------|
| [QIIME2](https://docs.qiime2.org/) | Must be installed and activated (`qiime` must be in PATH) |
| [biom-format](https://biom-format.org/) | For BIOM → TSV conversion (`biom convert`) |
| Python ≥ 3.8 | Standard library only — no additional Python packages required |

> The pipeline has been developed and tested with **QIIME2 2024.5**.

### QIIME2 plugins required

The following QIIME2 plugins must be installed in your environment:

- `qiime tools` (core)
- `qiime demux`
- `qiime cutadapt`
- `qiime dada2`
- `qiime feature-classifier`
- `qiime feature-table`
- `qiime taxa`
- `qiime metadata`
- `qiime diversity`
- `qiime composition` (for ANCOM-BC)

All of the above are included in standard QIIME2 distributions.

---

## Input files

| File | Description |
|------|-------------|
| **Manifest file** (TSV, V2 format) | Maps sample IDs to FASTQ file paths. Required for import. |
| **FASTQ files** | Raw sequencing reads — single-end or paired-end, Phred33-encoded. |
| **Metadata file** (TSV/TXT) | QIIME2-compatible sample metadata with a header row and a `sample-id` column. |
| **Pre-trained classifier** (`.qza`) | A QIIME2 Naive Bayes classifier trained on a reference database (e.g. SILVA 138 99% nb-classifier). |

### Manifest file format (V2)

**Paired-end:**
```
sample-id	forward-absolute-filepath	reverse-absolute-filepath
sample1	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz
sample2	/path/to/sample2_R1.fastq.gz	/path/to/sample2_R2.fastq.gz
```

**Single-end:**
```
sample-id	absolute-filepath
sample1	/path/to/sample1.fastq.gz
sample2	/path/to/sample2.fastq.gz
```

### Metadata file format

```
sample-id	group	timepoint	...
sample1	control	T0	...
sample2	treatment	T1	...
```

The `sample-id` column must match the sample IDs in the manifest. Additional columns are used for diversity grouping and differential abundance analysis.

### Pre-trained classifier

The classifier must be compatible with the QIIME2 version installed in your environment. Pre-trained classifiers for SILVA 138 are available from the [QIIME2 data resources page](https://docs.qiime2.org/2024.5/data-resources/).

The default path expected by the pipeline is:
```
/mnt/hdd2/database/qiime_pretrained_2024-5/silva-138-99-nb-classifier.qza
```
This can be changed interactively at runtime.

---

## Usage

```bash
# Activate your QIIME2 conda environment first
conda activate qiime2-2024.5

# Run the pipeline
python QMETA.py
```

The pipeline is fully interactive. At each step it will ask for input parameters and file paths. You can type `exit` at any prompt to safely quit.

---

## Pipeline steps

| Step | Description |
|------|-------------|
| **1 — Import** | Imports FASTQ reads into a QIIME2 artifact (`.qza`) using a V2 manifest file. Generates a demux summary visualization. |
| **2 — Cutadapt** *(optional)* | Trims primers and adapters using `qiime cutadapt trim-paired` or `trim-single`. Requires forward (and reverse, for paired-end) primer sequences. |
| **3 — DADA2** | Denoises reads and generates an ASV feature table (`table.qza`), representative sequences (`rep-seqs.qza`), and denoising statistics. Asks for truncation and trimming parameters. |
| **4 — Taxonomy** | Classifies ASVs using a pre-trained Naive Bayes classifier (`classify-sklearn`). Alternatively, an existing `taxonomy.qza` can be provided. |
| **5 — Taxonomic filter** | Removes Eukaryota and Unassigned features. Generates filtered feature table, rep-seqs, and taxa barplot visualizations. |
| **6 — Taxa collapse & export** | Collapses the feature table at taxonomic levels **2** (phylum), **5** (family), and **6** (genus). Exports each level as BIOM and TSV. Also generates relative frequency tables. |
| **7 — Diversity** | Runs non-phylogenetic `core-metrics`, per-level alpha diversity (7 metrics), alpha-group-significance tests, and alpha rarefaction curves. |
| **8 — ANCOM-BC** | Runs ANCOM-BC differential abundance analysis on collapsed tables at levels 2, 5, and 6. Supports specifying reference levels per comparison. Generates `da-barplot` visualizations. |

---

## Output structure

All outputs are saved in the `qiime2_out/` directory:

```
qiime2_out/
├── pipeline.log
├── import.qza
├── summary_import.qzv
├── import_cutadapt.qza             # if trimming was run
├── summary_import_cutadapt.qzv
├── table.qza
├── rep-seqs.qza
├── denoising-stats.qza / .qzv
├── table.qzv
├── table-no-eukaryota-no-unassigned.qza / .qzv
├── rep-seqs.qzv
├── core-metrics-results/
│   ├── ...
│   └── alpha-rarefaction.qzv
└── classification/
    ├── taxonomy_denoised.qza
    ├── taxonomy_denoised_metadata.qzv
    ├── taxa-bar-plots.qzv
    ├── level2/
    │   ├── collapsed_table.qza
    │   ├── level2_relative.qza / _tabulated.qzv
    │   ├── feature-table.biom
    │   ├── level2.tsv
    │   ├── alpha-diversity/
    │   └── DA/
    ├── level5/
    │   └── ...
    └── level6/
        └── ...
```

---

## Notes

- The pipeline does **not** perform phylogenetic analyses (no phylogenetic tree is built).
- Taxonomic collapse is performed at SILVA levels 2, 5, and 6 (phylum, family, genus). These can be changed by editing the `TAXONOMY_LEVELS` variable at the top of the script.
- At step 7, the sampling depth for rarefaction must be chosen based on the feature table summary (`table-no-eukaryota-no-unassigned.qzv`).
- ANCOM-BC reference levels must follow the format `column::value` (e.g. `group::control`).
