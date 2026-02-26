# Tabula Muris/Senis scRNAseq Reanalysis with Expanded Gene Annotations

Pipeline for reanalyzing [Tabula Muris](https://tabula-muris.ds.czbiohub.org/)/Tabula Senis single-cell RNA sequencing data with expanded gene annotations (Gencode M38), particularly to include lncRNAs. Currently limited to SmartSeq2/Plate-seq data from mice at 3, 18, and 24 months of age.

## Table of Contents

- [Processed Data](#processed-data)
- [Pipeline](#pipeline)
  - [Environment Setup](#environment-setup)
  - [Retrieve Data and Metadata](#retrieve-data-and-metadata)
  - [Retrieve and Prepare Annotation Files](#retrieve-and-prepare-annotation-files)
  - [Build Salmon Index and Align Reads](#build-salmon-index-and-align-reads)
  - [R Analysis Variables](#r-analysis-variables)
  - [TPMs to Cell x Gene Matrix](#tpms-to-cell-x-gene-matrix)
  - [QC and Normalize](#qc-and-normalize)
  - [Build SingleCellExperiment](#build-singlecellexperiment)
  - [Per-Cell-Class Analysis](#per-cell-class-analysis)
  - [Global PCA and UMAP](#global-pca-and-umap)
- [Plotting](#plotting)
- [References](#references)

---

## Processed Data

Processed SingleCellExperiment R data files are available for download on Zenodo:

**[https://zenodo.org/records/18765151](https://zenodo.org/records/18765151)**

---

## Pipeline

### Environment Setup

Environment YAML files are in the [`envs/`](envs/) folder. Create them with:

```bash
micromamba create -f envs/scrna_seq.yml -y
micromamba create -f envs/scrna_seq_R_env.yml -y
```

### Retrieve Data and Metadata

This grabs the SmartSeq2/Plate-seq data for 3, 18, and 24-month-old mice.

```bash
micromamba activate scrna_seq

mkdir Plate_seq && cd Plate_seq

for AGE in 3_month 18_month 24_month; do
  aws s3 sync s3://czb-tabula-muris-senis/Plate_seq/${AGE}/ ./${AGE}/ \
    --no-sign-request --exclude "*" --include "*.fastq.gz"
done

cd ..

aws s3 sync s3://czb-tabula-muris-senis/Metadata/ ./Metadata/ --no-sign-request
```

### Retrieve and Prepare Annotation Files

Get MM10PLUS files from the original publication to extract transgene and spike-in sequences.

```bash
mkdir annotations_genomes/ && cd annotations_genomes/

aws s3 cp s3://czb-tabula-muris-senis/reference-genome/MM10-PLUS.tgz . --no-sign-request
tar -xzf ./MM10-PLUS.tgz
grep -v "chr" MM10-PLUS/genes/genes.gtf > non_chromosomal_genes.gtf
awk '/^>/ {p = ($0 !~ /chr/)} p' MM10-PLUS/fasta/genome.fa > non_chromosomal_genes.fa
```

Get the latest Gencode release and mouse genome.

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.transcripts.fa.gz && gunzip gencode.vM38.transcripts.fa.gz

gunzip *.gz
```

Generate the combined transcriptome.

```bash
# cleanup gencode header
sed -i 's/|.*//' gencode.vM38.transcripts.fa

cat gencode.vM38.transcripts.fa non_chromosomal_genes.fa > gencode.vM38plus.fa
```

Make a tx2gene mapping file for downstream use.

```bash
grep -v '^#' gencode.vM38.primary_assembly.annotation.gtf | awk -F'\t' '$3 == "transcript" {
    split($9, attributes, "; ");
    txid = ""; gname = "";
    for (i in attributes) {
        if (attributes[i] ~ /^transcript_id/) {
            match(attributes[i], /"([^"]+)"/, arr);
            txid = arr[1];
        } else if (attributes[i] ~ /^gene_name/) {
            match(attributes[i], /"([^"]+)"/, arr);
            gname = arr[1];
        }
    }
    if (txid != "" && gname != "") {
        print txid "\t" gname;
    }
}' > gencodevM38plus.tx2gene.tsv

awk -F'\t' 'OFS="\t" { print $1, $1}' non_chromosomal_genes.gtf >> gencodevM38plus.tx2gene.tsv
```

### Build Salmon Index and Align Reads

Build the index and quantify. Pass `--geneMap` with a GTF to streamline downstream counting.

```bash
salmon index -t gencode.vM38plus.fa -p 16 -i gencode.vM38plus.index

salmon quant \
  -i annotations_genomes/gencode.vM38plus.index \
  -l A \
  -1 ${path1}_R1_001.fastq.gz \
  -2 ${path1}_R2_001.fastq.gz \
  -o ${outdir} \
  --validateMappings \
  --gcBias \
  --seqBias \
  --numBootstraps 30 \
  -p 1
```

### R Analysis Variables

Set threading variables (may be needed on some systems).

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
```

Set paths according to setup.

```bash
ENV=scrna_seq_R_env
SCRIPT_DIR=scripts
TX2GENE=annotations_genomes/gencodevM38plus.tx2gene.tsv
OUTPATH=counts/
INPATH=salmon_out/
METADATA_DIR=Metadata/
```

### TPMs to Cell x Gene Matrix

Convert Salmon quant.sf TPM files to a gene-count sparse matrix via tximport (batched in groups of 1000 samples).

```bash
micromamba run -p ${ENV} \
  Rscript ${SCRIPT_DIR}/single_cell_gene_counts.R \
  --tx2gene ${TX2GENE} \
  --quant_dir ${INPATH} \
  --outdir ${OUTPATH}
```

### QC and Normalize

Filter cells (default: min 500 genes, 50k reads) and normalize (log(1+CPM)).

```bash
RDS=${OUTPATH}/gene_counts.sparse.rds

micromamba run -p ${ENV} \
  Rscript ${SCRIPT_DIR}/QC_normalize_counts.R \
  --counts_rds ${RDS} \
  --outdir ${OUTPATH}
```

### Build SingleCellExperiment

Combine raw counts, logcounts, and published cell metadata into a SingleCellExperiment object.

```bash
micromamba run -p ${ENV} \
  Rscript ${SCRIPT_DIR}/build_sce_from_counts_and_metadata.R \
  --counts_rds    ${OUTPATH}/gene_counts.filtered.raw.rds \
  --logcounts_rds ${OUTPATH}/gene_counts.normalized_logcp.rds \
  --fastq_manifest ${METADATA_DIR}/tabula-muris-senis-facs-official-raw-obj__cell-metadata__cleaned_ids__read1_read2.csv \
  --metadata       ${METADATA_DIR}/tabula-muris-senis-facs-official-raw-obj__cell-metadata.csv \
  --outdir         ${OUTPATH}/ \
  --output_rds     sce_plate_seq.rds
```

### Per-Cell-Class Analysis

Per-cell-class HVG selection, PCA, and UMAP. Use `--list_cells` to enumerate available classes, then `--cell_classes` to select specific ones.

```bash
micromamba run -p ${ENV} \
  Rscript ${SCRIPT_DIR}/sce_per_cellclass_hvg_pca_umap.R \
  --sce_rds ${OUTPATH}/sce_plate_seq.rds \
  --list_cells

micromamba run -p ${ENV} \
  Rscript ${SCRIPT_DIR}/sce_per_cellclass_hvg_pca_umap.R \
  --sce_rds ${OUTPATH}/sce_plate_seq.rds \
  --cell_classes "microglial cell,B cell,hematopoietic stem cell,granulocyte,macrophage,monocyte,CD4+ aB T cell,CD8+ aB T cell,NK cell,myeloid dendritic cell,neutrophil" \
  --outdir ${OUTPATH} \
  --umap_min_dist 0.5 \
  --umap_neighbors 10
```

### Global PCA and UMAP

Gene filtering, HVG selection, PCA, UMAP, and Louvain clustering.

```bash
micromamba run -p ${ENV} \
  Rscript ${SCRIPT_DIR}/global_umap_clusters.R \
  --sce_rds ${OUTPATH}/sce_plate_seq.rds \
  --outdir  ${OUTPATH}/full_data_umap \
  --min_frac_cells 0.0000000001 \
  --min_mean_logexpr 0 \
  --n_hvg 3000 \
  --n_pcs 50 \
  --umap_neighbors 30 \
  --umap_min_dist 0.5 \
  --cluster_method louvain
```

---

## Plotting

Plotting utility functions are in `plotting_functions/plot_umap.R`. These are meant to be sourced interactively in an R session.

Example outputs:

![Global UMAP](figures/Global_labeled.png)

![Zeb2 Expression](figures/Zeb2_global.png)

---

## References

The Tabula Muris Consortium. Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. *Nature* 562, 367–372 (2018). [https://doi.org/10.1038/s41586-018-0590-4](https://doi.org/10.1038/s41586-018-0590-4)

The Tabula Muris Consortium. A single-cell transcriptomic atlas characterizes ageing tissues in the mouse. *Nature* 583, 590–595 (2020). [https://doi.org/10.1038/s41586-020-2496-1](https://doi.org/10.1038/s41586-020-2496-1)

Malekos, E., Smaliy, V., & Carpenter, S. (2026). Integrative multiomics analysis and CRISPR screening identify functional noncanonical translation loci in the mouse immune system. *bioRxiv* 2026.02.23.707583. [https://doi.org/10.64898/2026.02.23.707583](https://doi.org/10.64898/2026.02.23.707583)
