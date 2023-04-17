#Remember to activate snake env

import os, glob, sys
import pandas as pd

#Constants

RUN_DATE="20230413"
FLOW_CELL="FAR73398"
DSBIT_USER="efe2108"

CWD = os.getcwd()
TEMP = f"{CWD}/temp"
TEMP_RUN_DIR=f"{TEMP}/{RUN_DATE}.{FLOW_CELL}"
ARCHIVE_DIR=f"{CWD}/archives"
#ILLUMINA_READS_DIR=f"{TEMP_RUN_DIR}/illumina"
DATA_DIR=f"{CWD}/fastq_pass"
BARCODE_FILE=f"{CWD}/samples_barcodes.txt"

samples_df = pd.read_csv(BARCODE_FILE, sep='\t', header=None)
samples_df.columns=['sample', 'barcode']
samples_df = samples_df.set_index("sample", drop=False)

samples = list(samples_df.index)
barcode = list(samples_df['barcode'])

assemblers = ["flye", "raven", "miniasm", "canu"]

#Hardware parameters
NUM_THREADS_PIPELINE=12

#Software parameters
MIN_NANO_LEN=1000
MAX_NANO_HOMOP=20
TRYCYCLER_SUBSAMP_NUM=12

rule all:
	input:
		expand("{TEMP_RUN_DIR}/{sample}.{RUN_DATE}/post/cluster", sample = samples)
		
#Combine all fastqs from a given barcode into one file
rule concatenate:
	input: 
		lambda wc: glob.glob("{DATA_DIR}/{wc.barcode}/*.fastq.gz")
	output:
		"{ARCHIVE_DIR}/{RUN_DATE}.{FLOW_CELL}.{barcode}.{sample}.fastq.gz"
	shell:
		"cat {input} > {output}"

#Loose QC		
rule qc:
	input:
		lambda wc: "{ARCHIVE_DIR}/{RUN_DATE}.{FLOW_CELL}.{wc.barcode}.{wc.sample}.fastq.gz"
	output:
		"{TEMP_RUN_DIR}/{sample}.{barcode}.filter_len.{MIN_NANO_LEN}.fastq.gz"
	shell:
		"filtlong --min_length {MIN_NANO_LEN} --keep_percent 95 {input} | gzip > {output}"

#Multiple assemblies so Trycycler has something to work with	
rule flye:
	input:
		filt_fq = lambda wc: "{TEMP_RUN_DIR}/{wc.sample}.{wc.barcode}.filter_len.{MIN_NANO_LEN}.fastq.gz"
	output:
		"{TEMP_RUN_DIR}/{sample}.flye.{RUN_DATE}/assembly.fasta"
	shell:
		"flye --nano-hq {input.filt_fq} -o {TEMP_RUN_DIR}/{sample}.flye.{RUN_DATE} \ "
		"--threads {NUM_THREADS_PIPELINE} --meta --debug --read-error 0.03"
	
rule raven:
	input:
		lambda wc: "{TEMP_RUN_DIR}/{wc.sample}.{wc.barcode}.filter_len.{MIN_NANO_LEN}.fastq.gz"
	output:
		"{TEMP_RUN_DIR}/{sample}.raven.{RUN_DATE}/assembly.fasta"
	shell:
		"raven -t {NUM_THREADS_PIPELINE} {input} > {output}"

rule miniasm:
	input:
		lambda wc: "{TEMP_RUN_DIR}/{wc.sample}.{wc.barcode}.filter_len.{MIN_NANO_LEN}.fastq.gz"
	output:
		"{TEMP_RUN_DIR}/{sample}.miniasm.{RUN_DATE}/assembly.fasta"
	shell:
		"overlaps=$(mktemp)\".paf\"\n"
		"unpolished_assembly=$(mktemp)\".gfa\"\n"
		
		"minimap2 -x ava-ont -t {NUM_THREADS_PIPELINE} {input} {input} "
		"> \"$overlaps\"\n"
		
		"miniasm -f {input} \"$overlaps\" > \"$unpolished_assembly\"\n"
		
		"minipolish --threads {NUM_THREADS_PIPELINE} {input} \"$unpolished_assembly\" "
		"> {output}\n"
		
		"rm -f \"$overlaps\" \"$unpolished_assembly\""

rule canu:
	input:
		lambda wc: "{TEMP_RUN_DIR}/{wc.sample}.{wc.barcode}.filter_len.{MIN_NANO_LEN}.fastq.gz"
	output:
		"{TEMP_RUN_DIR}/{sample}.canu.{RUN_DATE}/assembly.fasta"
	shell:
		"canu -d {TEMP_RUN_DIR}/{sample}.canu.{RUN_DATE}/ -p {sample}.canu.{RUN_DATE} \ "
		"--correctedErrorRate 0.06 genomeSize=3.87m -nanopore -fast -corrected {input} \ "
		" > {output}"

#May need to unzip read file for this part
rule try_cluster:
	input:
		expand("{TEMP_RUN_DIR}/{sample}.{assemb}.{RUN_DATE}/assembly.fasta", assemb=assemblers)
	output:
		expand("{TEMP_RUN_DIR}/{sample}.{assemb}.{RUN_DATE}/cluster/contigs.phylip")
	shell:
		"trycycler cluster --assemblies {TEMP_RUN_DIR}/{sample}.{assemb}.{RUN_DATE}/assembly.fasta "
		"--reads {TEMP_RUN_DIR}/{wc.sample}.{wc.barcode}.filter_len.{MIN_NANO_LEN}.fastq.gz "
		"--out-dir {TEMP_RUN_DIR}/{sample}.{assemb}.{RUN_DATE}/cluster --threads {NUM_THREADS_PIPELINE}"
		
	
