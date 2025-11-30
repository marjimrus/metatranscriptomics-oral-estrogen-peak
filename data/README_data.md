# Data Access

## Raw Sequencing Data

Raw metatranscriptomic sequences have been deposited in the European Nucleotide Archive (ENA).

- **Accession**: [PRJEB101296](https://www.ebi.ac.uk/ena/browser/view/PRJEB101296)
- **Platform**: Illumina (paired-end RNA-seq)
- **Samples**: 20 samples (10 subjects × 2 timepoints)
- **Sample type**: Subgingival plaque (gingival margin)

### Download sequences from ENA

```bash
# Option 1: Download directly from ENA using wget/curl
# Visit https://www.ebi.ac.uk/ena/browser/view/PRJEB101296
# and download FASTQ files from the "Read Files" tab

# Option 2: Using enaBrowserTools
pip install enaBrowserTools
enaDataGet -f fastq PRJEB101296

# Option 3: Using Aspera (faster for large files)
# Install Aspera Connect, then use the Aspera links from ENA
```

## Sample Information

### Study Design

| Parameter | Description |
|-----------|-------------|
| Cohort | 10 healthy oocyte donors |
| Timepoints | T1 (before treatment), T2 (after controlled ovarian stimulation) |
| Sample type | Subgingival plaque from gingival margin (GM niche) |
| Design | Paired samples (same subject at T1 and T2) |

### Sample naming convention

Sample IDs follow the pattern: `ESTnnn_Tn`

- `EST001-EST010`: Subject identifiers
- `T1`: Before hormonal treatment (baseline)
- `T2`: After controlled ovarian stimulation (peak estrogen)

Example: `EST007_T1` = Subject 7, Timepoint 1 (before treatment)

### Metadata columns

| Column | Description |
|--------|-------------|
| SampleID | Unique sample identifier (ESTnnn.n format for analysis) |
| SubjectID | Subject identifier (ESTnnn) |
| Timepoint | 1 (T1/baseline) or 2 (T2/post-treatment) |

### Clinical data columns

| Column | Description |
|--------|-------------|
| subject_id | Subject identifier |
| age | Patient age at sampling |
| days_btw | Days between T1 and T2 sampling |
| d1_e2 | Day 1 blood estradiol (pg/mL) |
| t2_e2 | Final blood estradiol at T2 (pg/mL) |
| estrogen_sa_t1 | Salivary estrogen at T1 (pg/mL) |
| estrogen_sa_t2 | Salivary estrogen at T2 (pg/mL) |
| e2_change | Change in blood E2 (t2_e2 - d1_e2) |
| e2_change_normalized | Change per day (e2_change / days_btw) |

## Reference Databases

### Human Genome (for host removal)

- **Version**: GRCh38.p14 (GCF_000001405.40)
- **Source**: NCBI RefSeq
- **Purpose**: Remove human DNA/RNA contamination

```bash
# Download
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

# Build Bowtie2 index
bowtie2-build --threads 8 GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38
```

### GENCODE Annotation (for human alignment)

- **Version**: GENCODE v47
- **Source**: https://www.gencodegenes.org/human/
- **Purpose**: Human transcriptome quantification with STAR

```bash
# Download annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
```

### SortMeRNA Databases (for rRNA filtering)

- **SILVA 138.1**: 16S/18S and 23S/28S rRNA
- **RFAM 14.9**: 5S and 5.8S rRNA
- **Source**: https://github.com/biocore/sortmerna

```bash
# Databases included with SortMeRNA installation
# Or download separately from SILVA and RFAM
```

### HOMD Database (for microbial analysis)

- **Version**: Human Oral Microbiome Database V10.1
- **Source**: https://www.homd.org/
- **Content**: ~770 oral bacterial species
- **Annotations**: PROKKA gene predictions (GFF format)
- **Purpose**: Microbial transcriptome alignment and quantification

The HOMD database includes:
- `HOMD_genomes.fna`: Concatenated genome sequences
- `ALL_genomes.gff`: PROKKA gene annotations (GFF3 format)
- `ALL_genomes.faa`: Protein sequences for functional annotation
- `gene_strain_table.tsv`: Mapping of genes to species

#### Generating gene_strain_table.tsv

The `gene_strain_table.tsv` file (~1 GB) maps PROKKA gene IDs to their source species. This file is too large for GitHub but can be generated from the HOMD PROKKA annotations:

```bash
# Download HOMD genomes with PROKKA annotations from https://www.homd.org/

# Generate gene-to-species mapping from PROKKA GFF files
# Each GFF file contains genes for one genome/strain
for gff in HOMD_PROKKA/*.gff; do
    strain=$(basename "$gff" .gff)
    # Extract gene IDs and add strain name
    grep -P "^(?!#)" "$gff" | grep "CDS" | \
    awk -F'\t' -v strain="$strain" '{
        split($9, attrs, ";");
        for (i in attrs) {
            if (attrs[i] ~ /^ID=/) {
                gsub("ID=", "", attrs[i]);
                print attrs[i] "\t" strain
            }
        }
    }'
done > gene_strain_table.tsv
```

Alternatively, place the file in `data/gene_strain_table.tsv` after obtaining it from the corresponding author.

### eggNOG Database (for functional annotation)

- **Version**: eggNOG 5.0
- **Source**: http://eggnog5.embl.de/
- **Purpose**: KEGG, GO, COG annotation of bacterial genes

```bash
# Download eggNOG databases
download_eggnog_data.py --data_dir eggnog_data -y
```

## Output Files Description

### From Preprocessing

| File | Description |
|------|-------------|
| `*_trimmed_R1/R2.fastq.gz` | Quality-trimmed reads |
| `*_nonhost_R1/R2.fastq.gz` | Reads after human removal |
| `*_mRNA_R1/R2.fastq.gz` | mRNA-enriched reads (rRNA filtered) |
| `fastp_summary.tsv` | Trimming statistics per sample |
| `host_removal_summary.tsv` | Host removal statistics |
| `sortmerna_summary.tsv` | rRNA filtering statistics |

### From Alignment

| File | Description |
|------|-------------|
| `*_human.bam` | Human-aligned reads (STAR) |
| `*_homd.bam` | HOMD-aligned reads (Bowtie2) |
| `gene_counts_matrix.tsv` | Gene count matrix (featureCounts) |
| `merged_gene_counts.tsv` | Merged human gene counts |

### From Functional Annotation

| File | Description |
|------|-------------|
| `filtered_gene_ids.txt` | Gene IDs passing expression filter |
| `filtered_genes.faa` | Protein sequences for annotation |
| `eggnog_output.emapper.annotations` | Functional annotations |
| `eggnog_summary.txt` | Annotation statistics |

### From Statistical Analysis

| File | Description |
|------|-------------|
| `gsea_ko_annotated_keggrest.tsv` | GSEA results with KEGG annotations |
| `cytoscape_nodes_ko.tsv` | Node file for Cytoscape network |
| `cytoscape_edges_ko_similarity.tsv` | Edge file for Cytoscape network |
| `Fig_PCA.pdf` | PCA visualization |
| `Fig_Volcano.pdf` | Volcano plot of DGE results |
| `Fig_GSEA_KO.pdf` | GSEA enrichment plot |

## Quality Metrics

### Expected values for oral metatranscriptomic samples

| Metric | Expected Range | Notes |
|--------|----------------|-------|
| Q30 rate | >85% | High quality RNA-seq |
| Host retention | 30-70% | Depends on sample type |
| rRNA content | 40-80% | Normal for metatranscriptomics |
| HOMD alignment | 20-60% | Depends on microbial content |
| Gene detection | 10,000-50,000 | After filtering |

## Notes

- All patient data is anonymized
- Clinical metadata available upon reasonable request to corresponding author
- Sample collection approved by appropriate ethics committees
- Informed consent obtained from all participants
