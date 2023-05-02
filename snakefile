#Pipeline to create assemblies from a sequenced metagenome
#Remember to activate proper environment with conda activate [env]

import os, sys
import pandas as pd

### RUN DATA ###
RUN_DATE="20230413"
ONT_MODEL=""

### DIRECTORIES ###
CWD = os.getcwd()
OUT_DIR=f"{CWD}/out/{RUN_DATE}.{FLOW_CELL}"
ARCHIVE_DIR=f"""/media/uhlemannlab/Nano_Data/Nanopore/Metagenomics/\
{RUN_DATE}_MG/{CALL_DATE}_{ONT_MODEL}_{FLOW_CELL}_{HEX_ID}"""
BARCODE_FILE=f"{CWD}/selected_samples.csv"

samples_df = pd.read_csv(BARCODE_FILE, sep='\t', header=None)
samples_df.columns=['sample', 'barcode']
samples_df = samples_df.set_index("sample", drop=False)

samples = list(samples_df.index)
barcodes = list(samples_df['barcode'])

### SOFTWARE PARAMS ###
MIN_NANO_LEN=1000

rule all:
	input:
		expand(f"{OUT_DIR}/{{sample}}.assemblies/{{sample}}.final_assembly.fasta", 
		sample=samples)
		
#Loose QC		
rule qc:
	input:
		expand(f"{ARCHIVE_DIR}/{{sample}}.fastq.gz", sample=samples)
	output:
		f"{OUT_DIR}/{{sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz"
	shell:
		f"""filtlong --min_length {MIN_NANO_LEN} --keep_percent 95 {ARCHIVE_DIR}/\
{{wildcards.sample}}.fastq.gz | gzip > {{output}}"""

#Assemble with flye
rule flye:
	input:
		f"{OUT_DIR}/{{sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz"
	output:
		f"{OUT_DIR}/{{sample}}.assemblies/assembly.fasta"
	threads: 4
	shell:
		f"""flye --nano-hq {OUT_DIR}/{{wildcards.sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz \
-o {OUT_DIR}/{{wildcards.sample}}.assemblies \
--threads 4 --meta --debug --read-error 0.03"""

#Maybe remove short contigs using Seqtk?
#Maybe cut out duplicated sequences in circular contigs?

###Polishing: run several polishers in parallel, obtain consensus from results
###Apply two rounds of Medaka to consensus (Medaka doesn't rely on reads)

rule read_contig_overlap:
	input:
		f"{OUT_DIR}/{{sample}}.assemblies/assembly.fasta"
	output:
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/overlap.sam"
	threads: 2
	shell:
		f"""minimap2 -a -t 2 {{input}} \
{OUT_DIR}/{{wildcards.sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz > {{output}}"""

rule racon:
	input:
		f"{OUT_DIR}/{{sample}}.assemblies/assembly.fasta",
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/overlap.sam"
	output:
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/racon.fasta"
	threads: 2
	shell:
		f"""racon -t 2 {OUT_DIR}/{{wildcards.sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz \
{OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/overlap.sam \
{OUT_DIR}/{{wildcards.sample}}.assemblies/assembly.fasta"""

rule mini_align: 
	input:
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/racon.fasta"
	output:
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/align_out.bam"
	threads: 2
	shell:
		f"""mini_align -i {OUT_DIR}/{{sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz \
-r {{input}} -p {{output}} -t 2"""

rule med_consensus:
	input:
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/align_out.bam"
	output:
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/contig_consensus.hdf"
	threads: 2 
	shell:
		f"""medaka consensus {{input}} {{output}} --model {ONT_MODEL} --batch 200 \
--threads 2"""

rule med_stitch:
	input:
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/contig_consensus.hdf"
	output:
		f"{OUT_DIR}/{{sample}}.assemblies/{{sample}}.final_assembly.fasta"
	shell:
		f"medaka stitch {{input}} {{output}}"
		
#FMLRC2 to be implemented only if extra time

# rule bwt:
# 	input:
# 		f"{OUT_DIR}/{{sample}}.assemblies/assembly.fasta"
# 	output:
# 		
# 
# rule FMLRC2:
# 	input:
# 		f"{OUT_DIR}/{{sample}}.assemblies/assembly.fasta"
# 	output"
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/racon.fasta"
# 	threads:2
# 	shell:
# 		f"""
		
	
		
	
