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
ILLUMINA_READS_DIR=f"{TEMP_RUN_DIR}/illumina"
DATA_DIR=f"{CWD}/fastq_pass"
BARCODE_FILE=f"{CWD}/samples_barcodes.txt"

samples_df = pd.read_csv(BARCODE_FILE, sep='\t', header=None)
samples_df.columns=['sample', 'barcode']
samples_df = samples_df.set_index("sample", drop=False)

samples = list(samples_df.index)
barcode = list(samples_df['barcode'])

#Hardware parameters
NUM_THREADS_PIPELINE=12

#Software parameters
MIN_NANO_LEN=1000
MAX_NANO_HOMOP=20
TRYCYCLER_SUBSAMP_NUM=12


