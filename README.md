# meta-wgs-ont-pipeline
Pipeline for processing and assembling metagenomic WGS ONT data

## Construction outline:

### 1. Preprocessing
- filtlong/fastp
- PlasFlow to ID plasmid reads (?)

### 2. Assembly
- Flye, Raven, miniasm (run all three)
- Cluster and obtain consensus with Trycycler
- No subsetting: don't have enough coverage for low-abundance species

### 3. Polishing (avoiding reliance on Illumina reads)
Important choice: which/how many of these to use, and how many times to run them.
- Medaka
- Racon
- FMLRC2
- Clair3 (actually a variant caller)
- ntEdit
- NextPolish

### 4. Taxonomic ID with Kraken

### 5. Orientation and Annotation
- MUMmer to view alignment against reference (also good for evaluating quality)
- Prokka for gene annotation

### 6. Assessment
- QUAST for general assembly quality
- Prodigal --> predicted proteins, shorter length = more stop codons = more errors
