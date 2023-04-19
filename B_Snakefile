#Second of two Snakefiles for metagenomic assembly pipeline
#Picks up after Trycycler clustering and relevant human intervention
#Use FigTree to visualize clusters
#Use Bandage to judge circularity

#Requires cluster csv: in the format below
#001	"--linear"
#002	""
#...	...
#Where the left column is the cluster number and the right column is its linearity 

#Make sure appropriate environment (same as for first Snakefile) is active

import os, glob, sys
import pandas as pd

#Constants
RUN_DATE="20230413"
FLOW_CELL="FAR73398"
DSBIT_USER="efe2108"
ONT_MODEL=""

CWD = os.getcwd()
TEMP = f"{CWD}/temp"
TEMP_RUN_DIR=f"{TEMP}/{RUN_DATE}.{FLOW_CELL}"
ARCHIVE_DIR=f"{CWD}/archives"
#ILLUMINA_READS_DIR=f"{TEMP_RUN_DIR}/illumina"
DATA_DIR=f"{CWD}/fastq_pass"
BARCODE_FILE=f"{CWD}/samples_barcodes.txt"
CLUSTER_FILE=f"{CWD}/clusters.csv"

clusters_df = pd.read_csv(CLUSTER_FILE, sep='\t' header=None)
clusters_df.columns = ['clusters', 'linearity']
clusters = clusters_df['clusters'].tolist()

samples_df = pd.read_csv(BARCODE_FILE, sep='\t', header=None)
samples_df.columns=['sample', 'barcode']
samples_df = samples_df.set_index("sample", drop=False)

samples = list(samples_df.index)
barcode = list(samples_df['barcode'])

#Hardware parameters — could be replaced by CL args
NUM_THREADS_PIPELINE=12

rule all:
	input:
	
rule try_reconcile:
	input:
		expand("{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}", 
		clus=clusters, sample=samples)
	output:
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/2_all_seqs.fasta"
	params:
		circ=clusters_df.loc[clusters_df['clusters']==clus, 'linearity']	
		bc=samples_df.loc[samples_df['sample']==sample, 'barcode']
	shell:
		"trycycler reconcile --reads {TEMP_RUN_DIR}/{sample}.{params.bc}.filter_len.{MIN_NANO_LEN}.fastq.gz "
		"--threads {NUM_THREADS_PIPELINE} {params.circ} --verbose --cluster_dir {input}"
		
rule try_msa:
	input:
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/2_all_seqs.fasta"
	output:
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/3_msa.fasta"
	shell:
		"trycycler msa --cluster_dir {TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus} "
		"--threads {NUM_THREADS_PIPELINE}"
		
rule try_partition:
	input:
		expand("{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/2_all_seqs.fasta", 
		clus=clusters, sample=samples)
	params:
		bc=samples_df.loc[samples_df['sample']==sample, 'barcode']
	output:
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/4_reads.fastq"
	shell:
		"trycyler partition --reads {TEMP_RUN_DIR}/{sample}.{params.bc}.filter_len.{MIN_NANO_LEN}.fastq.gz "
		"--threads {NUM_THREADS_PIPELINE} --cluster_dirs "
		expand("{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}", clus=clusters)
		
rule try_consensus:
	input:
		expand("{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/3_msa.fasta",
		clus=clusters, sample=samples),
		expand("{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/4_reads.fastq",
		clus=clusters, sample=samples)
	output:
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/7_final_consensus.fasta"
	params:
		circ=clusters_df.loc[clusters_df['clusters']==clus, 'linearity']	
	shell:
		"trycycler consensus --cluster_dir {TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus} "
		"{params.circ} --verbose --threads {NUM_THREADS_PIPELINE}"
		
rule medaka:
	input:
		expand("{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/8_racon.fasta",
		clus=clusters, sample=samples)
	output:
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/9_medaka.fasta"
	shell:
		"medaka_consensus -i {TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/4_reads.fastq "
		"{input} -o {TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/ -m 
		"{ONT_MODEL} -t {NUM_THREADS_PIPELINE}\n"
		
		"mv {TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/consensus.fasta {output}"
		
rule racon:
	input:
		expand("{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/7_final_consensus.fasta",
		clus=clusters, sample=samples),
		expand("{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/overlap.sam",
		clus=clusters, sample=samples),
	output:
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/8_racon.fasta"
	shell:
		"racon -t {NUM_THREADS_PIPELINE} "
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/4_reads.fastq "
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/overlap.sam"
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/7_final_consensus.fasta"
		
rule read_contig_overlap:
	input:
		expand("{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/7_final_consensus.fasta",
		clus=clusters, sample=samples)
	output:
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/overlap.sam"
	shell:
		"minimap2 -a -t {NUM_THREADS_PIPELINE} {input} "
		"{TEMP_RUN_DIR}/{sample}.{RUN_DATE}.clusters/cluster_{clus}/4_reads.fastq "
		"> {output}"
	
#May not make sense to use, since other polishers may cause differences
#between reads and assembly. Instead: multiple rounds of medaka/racon?
#Consult Medini
rule FMLRC2:

	
		
		
		
		
	