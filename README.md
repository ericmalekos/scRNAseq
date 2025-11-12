# scRNAseq workflow

Using this repo to collect notes/pipeline for reanalyzing Tabula Muris/Tabula Senis data to expanded gene annotations, particularly to include lncRNAs.  
Currently limited to SmartSeq2/Plateseq pipeline and mice at 3 months.

## Quantification Steps

### set up environment 
```
    micromamba create -n scrna_seq \
    -c bioconda -c conda-forge "awscli>2" samtools salmon -y
```

### Retrieve data and metadata
This grabs the SmartSeq2/Plateseq data for 3 month old mice.  
Total size: ~20 Terabytes
```
    micromamba activate scrna_seq

    mkdir Plate_seq && cd Plate_seq

    aws s3 sync s3://czb-tabula-muris-senis/Plate_seq/3_month/ ./3_month/ --no-sign-request --exclude "*" --include "*.fastq.gz"

    cd ..

    aws s3 sync s3://czb-tabula-muris-senis/Metadata/ ./Metadata/ --no-sign-request

```


### Retrieve latest gencode annotation and prepare annotation files
Get MM10PLUS files from publication. We want to extract the transgenes and spike in sequences from these.

```
    mkdir annotations_genomes/ && cd annotations_genomes/

    aws s3 cp s3://czb-tabula-muris-senis/reference-genome/MM10-PLUS.tgz . --no-sign-request
    tar -xzf ./MM10-PLUS.tgz
    grep -v "chr" MM10-PLUS/genes/genes.gtf > non_chromosomal_genes.gtf
    awk '/^>/ {p = ($0 !~ /chr/)} p' MM10-PLUS/fasta/genome.fa > non_chromosomal_genes.fa
```

Get latest Gencode release and mouse genoem
```
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.annotation.gtf.gz
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.annotation.gtf.gz
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.transcripts.fa.gz && gunzip gencode.vM38.transcripts.fa.gz

    gunzip *.gz
```
Generate new transcriptome
```
    # cleanup gencode header
    sed -i 's/|.*//' gencode.vM38.transcripts.fa
    cat gencode.vM38.transcripts.fa non_chromosomal_genes.fa > gencode.vM38plus.fa
```

Make a tx2gene txt mapping file for later use
```
    grep -v '^#' gencode.vM38.primary_assembly.annotation.gtf | awk -F'\t' '$3 == "transcript" {
        # Split the attribute field (column 9) into key-value pairs
        split($9, attributes, "; ");

        # Initialize variables
        txid = "";
        gname = "";

        # Loop through the array and extract the specific attributes
        for (i in attributes) {
            if (attributes[i] ~ /^transcript_id/) {
                # Extract content between quotes for transcript_id
                match(attributes[i], /"([^"]+)"/, arr);
                txid = arr[1];
            } else if (attributes[i] ~ /^gene_name/) {
                # Extract content between quotes for gene_name
                match(attributes[i], /"([^"]+)"/, arr);
                gname = arr[1];
            }
        }

        # Print the TSV line
        if (txid != "" && gname != "") {
            print txid "\t" gname;
        }
    }' > gencodevM38plus.tx2gene.tsv

    awk -F'\t' 'OFS="\t" { print $1, $1}' non_chromosomal_genes.gtf >> gencodevM38plus.tx2gene.tsv
```


### Build Salmon index and align reads
Build index and process in parallel. System specific, basic Salmon command below.
```
    salmon index -t gencode.vM38plus.fa -p 16 -i gencode.vM38plus.index

    micromamba run -n scrna_seq \
    salmon quant \
    -i /hive/users/emalekos/tabula_muris_facs_bams/annotations_genomes/gencode.vM38plus_nodecoy.index \
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

Build index and process in parallel. System specific, basic Salmon command below.


## R single cell analysis steps
R scripts used in this section are in the scripts folder in this repo.

### Create a new environment

```
    micromamba create -y -n scrna_seq_R_env \
    -c conda-forge -c bioconda \
    r-base=4.4.* \
    bioconductor-tximport \
    bioconductor-scran \
    bioconductor-scater \
    bioconductor-scuttle \
    bioconductor-singlecellexperiment \
    bioconductor-summarizedexperiment \
    bioconductor-bluster \
    r-matrix r-data.table r-readr r-ggplot2 r-igraph r-uwot r-pkgconfig r-BiocManager
```

### Set variables and resources

### TPMs to CellxGene matrix
I had to set these due to multithreading errors on the server, you may not.
```
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export BLIS_NUM_THREADS=1
```

Set these according to your setup.
```
    ENV=scrna_seq_R_env
    SCRIPT_DIR=scripts
    TX2GENE=annotations_genomes/gencodevM38plus.tx2gene.tsv
    OUTPATH=counts/
    INPATH=salmon_out/
    METADATA_DIR=/Metadata/
```
### Read in salmon file and convert TPM to counts.
This will take while.  
Outputs include gene-count matrices.

```
    micromamba run -n ${ENV} \
    Rscript ${SCRIPT_DIR}/single_cell_gene_counts.R \
    --tx2gene ${TX2GENE} \
    --quant_dir ${INPATH} \
    --outdir ${OUTPATH}
```

### QC and normalize
See output for default filtering

```
    RDS=${OUTPATH}/gene_counts.sparse.rds
    micromamba run -n ${ENV} \
    Rscript ${SCRIPT_DIR}/QC_normalize_counts.R \
    --counts_rds ${RDS} \
    --outdir ${OUTPATH}

        Loading counts matrix from: /hive/users/emalekos/tabula_muris_facs_bams/Plate_seq/count_out//gene_counts.sparse.rds
        Raw counts matrix dimensions (before filtering): 77048 genes x 52315 cells
        Computing per-cell QC metrics (before filtering)...
        Filtering cells using thresholds: min_genes >= 500, min_counts >= 50000
        Keeping 38600 / 52315 cells (73.78%)
        Dropped 13715 cells; details in dropped_cells_qc.tsv
        Counts matrix dimensions (after filtering): 77048 genes x 38600 cells
        Computing gene detection stats (after filtering)...
        Gene detection summary (after filtering):
                            metric n_genes frac_genes
                            <char>   <int>      <num>
        1:                zero_cells    3268 0.04241512
        2:  expressed_in_<1pct_cells   45418 0.58947669
        3:  expressed_in_<5pct_cells   55843 0.72478195
        4: expressed_in_<10pct_cells   60104 0.78008514
        Computing cell detection stats (after filtering)...
        Cell detection summary (after filtering):
                                metric n_cells frac_cells pct_cells
                                <char>   <int>      <num>     <num>
        1:  cells_with_<50_genes_detected       0          0         0
        2: cells_with_<500_genes_detected       0          0         0
        Normalizing counts (CP1e+06 + log1p)...
```

### Merge with previously published cell annotation info

```
    micromamba run -n  ${ENV} \
    Rscript ${SCRIPT_DIR}/build_sce_from_counts_and_metadata.R \
    --counts_rds    ${OUTPATH}/gene_counts.filtered.raw.rds \
    --logcounts_rds ${OUTPATH}/gene_counts.normalized_logcp.rds \
    --fastq_manifest ${METADATA_DIR}/tabula-muris-senis-facs-official-raw-obj__cell-metadata__cleaned_ids__read1_read2.csv \
    --metadata       ${METADATA_DIR}/tabula-muris-senis-facs-official-raw-obj__cell-metadata.csv \
    --outdir         ${OUTPATH}/ \
    --output_rds     sce_plate_seq.rds
```

### PCA and UMAP 
see output for details
```
micromamba run -n ${ENV} \
  Rscript ${SCRIPT_DIR}/global_umap_clusters.R \
  --sce_rds ${OUTPATH}/sce_plate_seq.rds \
  --outdir  ${OUTPATH}/full_data_umap \
  --min_frac_cells 0.001 \
  --min_mean_logexpr 0 \
  --n_hvg 3000 \
  --n_pcs 50 \
  --umap_neighbors 30 \
  --umap_min_dist 0.5 \
  --cluster_method louvain


    Loading SCE from: /hive/users/emalekos/tabula_muris_facs_bams/Plate_seq/count_out//sce_plate_seq.rds
    Initial SCE dimensions: 77048 genes x 38600 cells
    Applying global gene filter...
    Keeping 46576 / 77048 genes (>=0.001 frac cells, mean logexpr > 0.000)
    After filtering (if any): 46576 genes x 38600 cells
    Computing HVGs...
    Selected 3000 HVGs.
    Running PCA on HVGs (50 PCs)...
    PCA completed.
    Running UMAP (n_neighbors = 30, min_dist = 0.5, name = 'UMAP_refined')...
    UMAP completed.
    Building SNN graph (k = 20)...
    Clustering with method: louvain
    Clustering completed: 34 clusters (stored in colData(sce)$cluster_louvain).
```

### Tissue level
If preference is for only tissue level analysis:

```
    micromamba run -n ${ENV} \
    Rscript ${SCRIPT_DIR}/sce_per_tissue_hvg_pca_umap.R \
    --sce_rds ${OUTPATH}/sce_plate_seq.rds \
    --outdir  ${OUTPATH}/tissue_umaps \
    --umap_neighbors 15 \
    --umap_min_dist 0.2 \
    --tissues "Brain_Myeloid,Spleen,Marrow,Brain_Non_Myeloid,Thymus,Liver"
```

### Genereate figures
Functions for generating figures are in `scripts/sce_per_tissue_hvg_pca_umap.R`   
These are meant to be interactive and don't have a CLI (yet, and maybe never).  
Example outputs:

![Global](figures/Global_labeled.png)

![Zeb2](figures/Zeb2_global.png)
