# meta-wgs-ont-pipeline
Pipeline for processing and assembling metagenomic WGS ONT data

## Construction outline:

### 1. Preprocessing with Filtlong
- Relevant parameters: MIN_NANO_LEN
- Recommend MIN_NANO_LEN > 200


### 2. Assembly using Flye

### 3. Polishing (avoiding reliance on Illumina reads)
- Two rounds of Medaka
- Generate polished versions with:
  - Racon

To add:
- Clair3 (actually a variant caller)
- ntEdit
- NextPolish
- FMLRC2

### 4. Obtain consensus between various polishers (final assembly)

### 5. Assessment with QUAST
