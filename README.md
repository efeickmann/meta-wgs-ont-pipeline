# meta-wgs-ont-pipeline
Pipeline for processing and assembling metagenomic WGS ONT data

## File descriptions
- snakefile: Main pipeline file. Contains directions for assembly.
- meta_config.yml: Configuration file with instructions to create the necessary Conda environment.

## Usage
1. Make a local clone of the repo:
```
git clone https://github.com/efeickmann/meta-wgs-ont-pipeline.git \
&& cd meta-wgs-ont-pipeline
```

2. Use the provided config file (```meta_pipe.yml```) to create a Conda environment (Recommend using Mamba. Instructions for installation at https://mamba.readthedocs.io/en/latest/):
```
mamba create -n meta_pipe_env --file meta_pipe.yml
```
Activate the environment with ```mamba activate meta_pipe_env```

3. Open ```snakefile``` in your favorite text editor and edit the following parameters:
- RUN_DATE: The date your ONT data was produced, in the form YYYYMMDD. Used to ID the sequencing run. Use (or not) as necessary for personal organization.
- FLOW_CELL: The ID of the flowcell used for the run. Use (or not) as necessary for personal organization.
- ONT_MODEL: The chemistry, flowcell type, and basecaller used to generate your data. Should be given in the form ```r1041_min_high_g303``` for data produced with R10.4.1 chemistry on a minION flowcell and basecalled with Guppy v3.0.3 set to high accuracy basecalling. Use ```medaka medaka tools list\_models``` to view a list of all supported combinations.
- ARCHIVE_DIR: The directory containing the basecalled fastq files from your run. 
- BARCODE_FILE: A .tsv file with a list of the samples you'd like to process, along with their barcodes. Repo comes with an example file (```selected_samples.tsv```). If you use your own, it should be in the form:
```
sample   barcode
31278   barcode01
...   ...
```
- MIN_NANO_LEN: The minimum read length you'd like to keep.
- Threads: Each rule, if computationally expensive, has a 'threads' parameter. Lowering this parameter will cause rules to run more slowly and in some cases fail (e.g. flye requires >= 4 threads). Feel free to raise the number of threads used, at the cost of reduced parallelism in the case of processing multiple samples.

4. In the terminal, run ```snakemake -c[CORES]```, replacing [CORES] with however many CPUs you have available.

5. Wait.

6. Visualize your assembly quality using the .html summary file generated by QUAST.

## Construction outline:

### 1. Preprocessing with Filtlong
- Relevant parameters: MIN_NANO_LEN
- Recommend MIN_NANO_LEN > 200

### 2. Assembly using Flye

### 3. Polishing (avoiding reliance on Illumina reads) with Racon

### 4. Polishing with Medaka

### 5. (Recommended) Assessment with metaQUAST [not included]
