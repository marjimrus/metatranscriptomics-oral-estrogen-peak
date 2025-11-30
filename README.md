# Oral Metatranscriptomics Analysis: Microbial Resilience to Estradiol Surges

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.3-blue)](https://www.r-project.org/)
[![DESeq2](https://img.shields.io/badge/DESeq2-1.44-green)](https://bioconductor.org/packages/DESeq2/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17416734.svg)](https://doi.org/10.5281/zenodo.17416734)

Analysis pipeline for oral metatranscriptomic sequencing (RNA-seq) of subgingival plaque samples from oocyte donors before and after controlled ovarian hyperstimulation (COH).

> **Publication**: *Microbial Resilience and Functional Gene Adaptation to Acute Estradiol Surges in the Human Oral Cavity* (manuscript under review)

## Overview

This repository contains the bioinformatics pipeline from raw reads to differential gene expression analysis, demonstrating functional resilience of the oral microbiome to acute estrogen fluctuations, with taxon-specific responses in *Prevotella* species.

- **Preprocessing pipeline** with FastQC, fastp, Bowtie2 (host removal), and SortMeRNA (rRNA filtering)
- **Microbial alignment** using Bowtie2 + featureCounts with HOMD PROKKA annotations
- **Statistical analysis** in R with DESeq2 (paired design)

## Study Design

| Parameter | Description |
|-----------|-------------|
| Study focus | Effect of acute estradiol surges on oral microbiome gene expression |
| Model | Controlled Ovarian Hyperstimulation (COH) for oocyte donation |
| Cohort | 10 healthy women (18-35 years) undergoing COH |
| Sample type | Subgingival plaque (gingival margin - GM niche) |
| Timepoints | T1 (baseline, before COH) vs T2 (peak estradiol, oocyte retrieval day) |
| Design | Paired samples (same individual at T1 and T2) |
| Estradiol levels | T1: ~50 pg/mL (baseline) → T2: ~2000-4000 pg/mL (peak) |
| Platform | Illumina NovaSeq 6000 (paired-end 150bp RNA-seq) |

Controlled ovarian hyperstimulation provides a unique model for studying acute hormonal effects: standardized intervention, 40-80 fold increase in serum estradiol over 10-14 days, healthy subjects without underlying pathology, and precise timing of baseline and peak hormone sampling.

This project uses RNA-seq rather than DNA-based metagenomics to capture what the microbiome is actively doing (functional activity), with high temporal sensitivity to environmental stimuli like hormones.

## Key Findings

### Main Result: Functional Resilience

The oral microbiome demonstrates remarkable functional stability during acute estradiol fluctuations:

| Analysis | Result |
|----------|--------|
| Genes tested | 172,851 (after filtering: ≥10 counts in ≥5 samples) |
| Significant DEGs | **0** (padj < 0.05 after Benjamini-Hochberg correction) |
| Nominal DEGs | 8,893 genes with p < 0.05 (before FDR correction) |
| Interpretation | Functional resilience at community level |

### Taxon-Specific Responses

Despite overall community resilience, *Prevotella* species showed taxon-specific transcriptional responses:

- *P. oris* and *P. melaninogenica*: Gene expression correlated with salivary estradiol changes
- Biological relevance: *Prevotella* species are known estrogen metabolizers (beta-glucuronidase activity)

### PCA Clustering

- Samples cluster by SubjectID rather than Timepoint
- High inter-individual variation dominates the signal
- Justifies the paired statistical design

## Repository Structure

```
metatranscriptomics-oral-estrogen-peak/
├── README.md
├── LICENSE
├── scripts/
│   ├── preprocessing/
│   │   ├── 01_quality_control.sh        # FastQC + MultiQC
│   │   ├── 02_fastp_trimming.sh         # Quality trimming
│   │   ├── 03_host_removal.sh           # Human DNA/RNA removal
│   │   └── 04_rrna_filtering.sh         # rRNA removal with SortMeRNA
│   ├── alignment/
│   │   └── 05_bowtie2_microbiome.sh     # Microbial alignment + featureCounts
│   ├── 06_clinical_analysis.Rmd         # Clinical data analysis
│   └── 07_metatranscriptome_analysis.Rmd # Differential gene expression
└── data/
    ├── README_data.md                   # Data access instructions
    ├── metadata.tsv                     # Sample metadata (anonymized)
    └── clinical_data.txt                # Clinical variables (anonymized)
```

## Pipeline Overview

```
Raw FASTQ files (Illumina paired-end RNA-seq)
        │
        ▼
┌───────────────────────────────────────┐
│  PREPROCESSING                        │
│  ├── FastQC (quality assessment)      │
│  ├── fastp (Q20, 70bp min, cut_front) │
│  ├── Bowtie2 (human removal, GRCh38)  │
│  └── SortMeRNA (rRNA filtering)       │
└───────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────┐
│  ALIGNMENT                            │
│  └── Bowtie2 + featureCounts          │
│      └── HOMD PROKKA annotations      │
└───────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────┐
│  STATISTICAL ANALYSIS (R)             │
│  └── DESeq2 (paired design)           │
│      └── ~ SubjectID + Timepoint      │
└───────────────────────────────────────┘
```

## Methods

### 1. Preprocessing

| Step | Tool | Parameters |
|------|------|------------|
| Quality assessment | FastQC v0.12.1 | Default |
| Report aggregation | MultiQC v1.15 | Default |
| Quality trimming | fastp v0.23.4 | Q20, min-length 70bp, cut_front 10bp, adapter auto-detection |
| Host removal | Bowtie2 v2.5.4 | --very-sensitive, GRCh38 reference |
| rRNA filtering | SortMeRNA v4.3.7 | SILVA 138.1 + RFAM 14.9 databases |

### 2. Alignment

| Step | Tool | Parameters |
|------|------|------------|
| Microbiome alignment | Bowtie2 v2.5.4 | HOMD RefSeq v10.1 reference, --very-sensitive |
| Gene quantification | featureCounts v2.0.6 | PROKKA v1.14.6 GFF annotations, paired-end, -p -B -C |

### 3. Statistical Analysis (R)

| Analysis | Method | Purpose |
|----------|--------|---------|
| Normalization | VST (DESeq2) | Variance-stabilizing transformation for visualization |
| DGE analysis | DESeq2 v1.44.0 | Paired design: ~ SubjectID + Timepoint |
| Gene filtering | ≥10 counts in ≥5 samples | Remove low-expression genes |
| Multiple testing | Benjamini-Hochberg | FDR correction (padj < 0.05) |

## Requirements

### Key Dependencies

**Preprocessing:**
- FastQC v0.12.1
- MultiQC v1.15
- fastp v0.23.4
- Bowtie2 v2.5.4
- SortMeRNA v4.3.7
- samtools v1.18

**Alignment:**
- featureCounts (Subread v2.0.6)

**Statistical Analysis (R ≥4.3):**
```r
library(DESeq2)         # v1.44.0 - Differential expression
library(tidyverse)      # Data manipulation
library(pheatmap)       # Heatmaps
library(EnhancedVolcano) # Volcano plots
library(ggpubr)         # Publication plots
```

### Database Requirements

1. **Human genome** (for host removal):
   ```bash
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
   bowtie2-build --threads 8 GRCh38.fna GRCh38
   ```

2. **HOMD database** (for microbial alignment):
   ```bash
   # Download HOMD V10.1 genomes with PROKKA annotations from: https://www.homd.org/
   bowtie2-build --threads 8 HOMD_genomes.fna HOMD_v10.1
   ```

## Quick Start

```bash
# Preprocessing
bash scripts/preprocessing/01_quality_control.sh
bash scripts/preprocessing/02_fastp_trimming.sh
bash scripts/preprocessing/03_host_removal.sh
bash scripts/preprocessing/04_rrna_filtering.sh

# Alignment
bash scripts/alignment/05_bowtie2_microbiome.sh
```

```r
# Statistical analysis
rmarkdown::render("scripts/06_clinical_analysis.Rmd")
rmarkdown::render("scripts/07_metatranscriptome_analysis.Rmd")
```

## Results Interpretation

The absence of significant differentially expressed genes after FDR correction is a positive finding:

1. **Functional resilience**: The oral microbiome maintains stable gene expression despite extreme hormonal fluctuations
2. **Homeostatic mechanisms**: Suggests robust compensatory mechanisms at the community level
3. **Clinical relevance**: COH may not disrupt oral microbiome function, reassuring for fertility patients

While the community shows resilience, individual taxa (*Prevotella*) respond to estrogen, supporting the "estrobolome" concept of microbes metabolizing estrogen.

## Data Access

Raw RNA sequences and metatranscriptomic metadata have been deposited in the European Nucleotide Archive (ENA) under accession number [PRJEB101296](https://www.ebi.ac.uk/ena/browser/view/PRJEB101296).

See `data/README_data.md` for detailed data access instructions.

## Citation

If you use this pipeline or data, please cite:

> *Microbial Resilience and Functional Gene Adaptation to Acute Estradiol Surges in the Human Oral Cavity* (manuscript under review)

**Code and Supplementary Data**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17416734.svg)](https://doi.org/10.5281/zenodo.17416734)

## Ethics

This study was approved by the Research Ethics Committee of the Virgen Macarena and Virgen del Rocío University Hospitals (CEI-24/2022). All participants provided written informed consent.

## Funding

This work was supported by the Spanish Ministry of Science and Innovation (Grant PID2020-118557GA-I00).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Maria J. Rus**
[![ORCID](https://img.shields.io/badge/ORCID-0000--0003--3659--2821-green)](https://orcid.org/0000-0003-3659-2821)
Email: marjimrus@gmail.com

## Acknowledgments

- Universidad de Sevilla
- University of Leeds (collaboration)
- Human Oral Microbiome Database (HOMD)
