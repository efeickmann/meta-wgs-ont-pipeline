#Pipeline to create assemblies from a sequenced metagenome
#Remember to activate proper environment with conda activate [env]

import os, sys
import pandas as pd

### RUN DATA ###
RUN_DATE="20230501"
FLOW_CELL="..."
ONT_MODEL="r941_min_sup_g507"

### DIRECTORIES ###
CWD = os.getcwd()
OUT_DIR=f"{CWD}/SRR17913199"
#OUT_DIR=f"{CWD}/out/{RUN_DATE}.{FLOW_CELL}"
ARCHIVE_DIR=f"{CWD}/reads"
BARCODE_FILE=f"{CWD}/selected_samples.tsv"

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
	# conda:
# 		"filtlong_env"
	shell:
		f"""filtlong --min_length {MIN_NANO_LEN} --keep_percent 95 {ARCHIVE_DIR}/\
{{wildcards.sample}}.fastq.gz | gzip > {{output}}"""

#Assemble with flye
rule flye:
	input:
		f"{OUT_DIR}/{{sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz"
	output:
		f"{OUT_DIR}/{{sample}}.assemblies/assembly.fasta"
	# conda:
# 		"flye_env"
	threads: 14
	shell:
		f"""flye --nano-hq {OUT_DIR}/{{wildcards.sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz \
-o {OUT_DIR}/{{wildcards.sample}}.assemblies \
--threads {{threads}} --meta --read-error 0.03"""

###Polishing: pass results from Racon to Medaka

#Racon rules
rule read_contig_overlap:
	input:
		f"{OUT_DIR}/{{sample}}.assemblies/assembly.fasta"
	output:
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/overlap.sam"
	# conda:
# 		"racon_env"
	threads: 14
	shell:
		f"""minimap2 -a -t {{threads}} {{input}} \
{OUT_DIR}/{{wildcards.sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz > {{output}}"""

rule racon:
	input:
		f"{OUT_DIR}/{{sample}}.assemblies/assembly.fasta",
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/overlap.sam"
	output:
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/{{sample}}.racon.fasta"
	# conda:
# 		"racon_env"
	threads: 14
	shell:
		f"""racon -t {{threads}} {OUT_DIR}/{{wildcards.sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz \
{OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/overlap.sam \
{OUT_DIR}/{{wildcards.sample}}.assemblies/assembly.fasta > {{output}}"""

#Medaka
rule medaka:
	input:
		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/{{sample}}.racon.fasta"
	output:
		f"{OUT_DIR}/{{sample}}.assemblies/{{sample}}.final_assembly.fasta"
	# conda:
# 		"medaka_env"
	threads: 14
	shell:
		f"""medaka_consensus -i {OUT_DIR}/{{wildcards.sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz \
-d {{input}} -o {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp -t {{threads}} -m {ONT_MODEL}
mv {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/consensus.fasta {{output}}"""

### Entering Area Under Construction ###

# #ntEdit rules
# rule nthits:
# 	input:
# 		f"{OUT_DIR}/{{sample}}.assemblies/assembly.fasta"
# 	output:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/solidBF_k40.bf"
# 	threads: 6
# 	shell:
# 		f"""cd {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/
# nthits -c 1 --outbloom -p solidBF -b 36 -k 40 -t {{threads}} {{input}}
# cd {CWD}"""
# 
# rule ntedit:
# 	input:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/solidBF_k40.bf"
# 	output:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/{{sample}}.ntedit.fasta"
# 	threads: 6
# 	shell:
# 		f"""cd {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/
# ntedit -f {OUT_DIR}/{{wildcards.sample}}.assemblies/assembly.fasta \
# -r {{input}} -t {{threads}}
# mv assembly.fasta_k40_z100_rsolidBF_k40.bf_i4_d5_m0_edited.fa {{output}}
# cd {CWD}"""

# NextPolish has dependencies which are incompatible with other tools
# Consequently, it will need to be added later
# #NextPolish rules
# rule nextpolish:
# 	input:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/{{sample}}.medaka.fasta"
# 	output:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/{{sample}}.nextpol.fasta"
# 	threads: 4
# 	shell:
# 		f"""ls {OUT_DIR}/{{wildcards.sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz > 
# {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/lgs.fofn
		

# #Consensus rules
# rule count_tigs:
# 	input: 
# 		f"{OUT_DIR}/{{sample}}.assemblies/assembly.fasta"
# 	output:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/tig_count.txt"
# 	shell:
# 		f"""wc -l {OUT_DIR}/{{wildcards.sample}}.assemblies/assembly_info.txt > \
# {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/tig_count.txt"""
# 
# rule cat_polish:
# 	input:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/{{sample}}.racon.fasta",
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/{{sample}}.ntedit.fasta",
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/{{sample}}.medaka.fasta",
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/tig_count.txt"
# 	output:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/clusters/cluster_1/2_all_seqs.fasta"
# 	shell:
# 		f"count=$(cut -d \" \" -f1 {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/tig_count.txt)\n"
# 		"for clstr in $(seq 1 $count)\n"
# 		"do\n"
# 		f"	mkdir {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/" 
# 		"cluster_${{clstr}}\n"
# 		f"	mkdir {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/" 
# 		"cluster_${{clstr}}/1_contigs\n"
# 		"done\n"
# 
# 		f"{CWD}/cluster.sh {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/{{wildcards.sample}}.racon.fasta \
# \'{OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_\'\n"
# 
# 		f"{CWD}/cluster.sh {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/{{wildcards.sample}}.medaka.fasta \
# \'{OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_\'\n"
# 
# 		f"{CWD}/cluster.sh {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/{{wildcards.sample}}.ntedit.fasta \
# \'{OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_\'\n"
# 
# 		"for clstr in $(seq 1 $count)\n"
# 		"do\n"
# 		f"	cat {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/" "cluster_${{clstr}}"
# 		f"""/1_contigs/* > {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_$count/2_all_seqs.fasta
# done"""
# 	
	# echo \"contig_\"$clstr > tmp.txt
# 	python {CWD}/faSomeRecords.py --fasta {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/{{wildcards.sample}}.racon.fasta \
# --list tmp.txt -o {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_$clstr/1_contigs/{{wildcards.sample}}.racon.fasta
# 	python {CWD}/faSomeRecords.py --fasta {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/{{wildcards.sample}}.ntedit.fasta \
# --list tmp.txt -o {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_$clstr/1_contigs/{{wildcards.sample}}.ntedit.fasta
# 	python {CWD}/faSomeRecords.py --fasta {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/{{wildcards.sample}}.medaka.fasta \
# --list tmp.txt -o {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_$clstr/1_contigs/{{wildcards.sample}}.medaka.fasta
# 	rm -f tmp.txt
#done
# cat {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_$count/1_contigs/* \
# > {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_$count/2_all_seqs.fasta"""

# rule try_msa:
# 	input:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/clusters/cluster_1/2_all_seqs.fasta"
# 	output:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/clusters/cluster_1/3_msa.fasta"
# 	threads: 12
# 	shell:
# 		f"""for clstr in $(seq 1 $count)
# do
# 	trycycler msa --cluster_dir {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_$clstr \
# --threads {{threads}}
# done"""
# 
# rule try_partition:
# 	input:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/clusters/cluster_1/2_all_seqs.fasta"
# 	output:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/clusters/cluster_1/4_reads.fastq"
# 	threads: 12
# 	shell:
# 		f"""trycycler partition --reads {OUT_DIR}/{{wildcards.sample}}.filter_len.{MIN_NANO_LEN}.fastq.gz \
# --threads {{threads}} --cluster_dirs {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_*"""
# 
# rule try_consensus:
# 	input:
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/clusters/cluster_1/4_reads.fastq",
# 		f"{OUT_DIR}/{{sample}}.assemblies/polish_temp/clusters/cluster_1/3_msa.fasta"
# 	output:
# 		f"{OUT_DIR}/{{sample}}.assemblies/{{sample}}.final_assembly.fasta"
# 	threads: 12
# 	shell:
# 		f"""for clstr in $(seq 1 $count)
# do
# 	trycycler consensus --cluster_dir {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_$clstr\
# --linear --verbose --threads {{threads}}
# done
# cat {OUT_DIR}/{{wildcards.sample}}.assemblies/polish_temp/clusters/cluster_*/7_final_consensus.fasta \
# > {OUT_DIR}/{{wildcards.sample}}.assemblies/{{wildcards.sample}}.final_assembly.fasta"""
		
#FMLRC2 to be implemented only if extra time — it is written in Rust which is a whole
#other beast
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
		
	
		
	
