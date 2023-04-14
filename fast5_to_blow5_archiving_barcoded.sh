######## ONLY CHANGE THESE VARIABLES ##################
######## Run 'conda activate slow5tools' first ########

NANO_RUN="20211022_U54_LUC"
FLOW_CELL="FAR05758"
NUM_THREADS=10
FAST5_DIR="./fast5_pass"
FASTQ_DIR="./fastq_pass"
BLOW5_DIR="./blow5_pass" # must not already exist
ARCHIVE_DIR_NAME="archives"

#######################################################


######## Step 1. Defining a function to check that BLOW5 conversion is accurate ###################
######## Based on counting the number of reads in BLOW5 file and comparing to FAST5 file ##########

sanity_check_fast5_num_reads(){

	FAST5_BARCODE_DIR=$1
	NUM_SLOW5_READS=$2

	test -d  ${FAST5_BARCODE_DIR} || die "$FAST5_BARCODE_DIR not found"

	if parallel --version > /dev/null
	then
		NUM_READS=$(find $FAST5_BARCODE_DIR -name '*.fast5' | parallel -I% --max-args 1 strings % | grep "read_.*-.*-.*" | wc -l | awk 'BEGIN {count=0;} {count=count+$0} END {print count;}')
	else
		NUM_READS=$(find $FAST5_BARCODE_DIR -name '*.fast5' | xargs --max-args 1 strings | grep "read_.*-.*-.*" | wc -l | awk 'BEGIN {count=0;} {count=count+$0} END {print count;}')
	fi

	if [ ${NUM_READS} -ne ${NUM_SLOW5_READS} ]
	then
		echo "ERROR: Sanity check has failed. $NUM_READS reads in $FAST5_BARCODE_DIR FAST5, but $NUM_SLOW5_READS reads in SLOW5/BLOW5" >> ../${ARCHIVE_DIR_NAME}/${NANO_RUN}.sanity_check_output.txt
		exit 1
	else
		echo "Sanity check has passed. $NUM_READS reads in $FAST5_BARCODE_DIR FAST5, $NUM_SLOW5_READS reads in SLOW5/BLOW5" >> ../${ARCHIVE_DIR_NAME}/${NANO_RUN}.sanity_check_output.txt
	fi

}

######## Step 2. Collating FASTQ files per barcode and a run FASTQ file for archiving #############

mkdir ${ARCHIVE_DIR_NAME}

cd ${FASTQ_DIR}

for FASTQ_BARCODE_DIR in */; do
  echo "Gzipping FASTQ files for ${FASTQ_BARCODE_DIR}"
  cd ${FASTQ_BARCODE_DIR}
  BARCODE=${FASTQ_BARCODE_DIR%*/}
  gzip "${FLOW_CELL}"*.fastq
  cat "${FLOW_CELL}"*.fastq.gz > ../../${ARCHIVE_DIR_NAME}/"${NANO_RUN}"."${FLOW_CELL}"."${BARCODE}".archive.fastq.gz
  cd ..
done

cd ..

cat ./${ARCHIVE_DIR_NAME}/"${NANO_RUN}"."${FLOW_CELL}"*.archive.fastq.gz > ./${ARCHIVE_DIR_NAME}/${NANO_RUN}.archive.fastq.gz

######## Step 3. Converting FAST5 files per barcode into BLOW5, merging, and running sanity check #########

slow5tools f2s --retain ${FAST5_DIR} -d ${BLOW5_DIR} -p ${NUM_THREADS}

cd ${BLOW5_DIR}

for BLOW5_BARCODE_DIR in */; do
    BARCODE=${BLOW5_BARCODE_DIR%*/}
    FAST5_BARCODE_DIR=../${FAST5_DIR#"./"}/${BARCODE%*/}
    echo "Generating BLOW5 archive file for ${BARCODE}"
    slow5tools merge ${BLOW5_BARCODE_DIR} -o ../${ARCHIVE_DIR_NAME}/${NANO_RUN}_${BARCODE}.archive.blow5 -t ${NUM_THREADS}
    echo "Sanity check for FAST5 & SLOW5/BLOW5 read counts in ${BARCODE}"
    NUM_SLOW5_READS=$(slow5tools stats ../${ARCHIVE_DIR_NAME}/${NANO_RUN}_${BARCODE}.archive.blow5 | grep "number of records" | awk '{print $NF}')
    sanity_check_fast5_num_reads ${FAST5_BARCODE_DIR} ${NUM_SLOW5_READS}
done

cd ..

slow5tools merge ./${ARCHIVE_DIR_NAME}/ -o ./${ARCHIVE_DIR_NAME}/${NANO_RUN}.archive.blow5 -t ${NUM_THREADS}

######## Step 4. Check if all is well and if you can delete fast5 folder ###############

NUM_BARCODES=$(find ./${ARCHIVE_DIR_NAME} -name '*barcode*.blow5' | wc -l)
NUM_BARCODES_PASSED=$(awk '/passed/ {count++} END{print count}' ./${ARCHIVE_DIR_NAME}/${NANO_RUN}.sanity_check_output.txt)

echo "${NUM_BARCODES} barcodes detected." 
echo "${NUM_BARCODES_PASSED} barcodes passed conversion to BLOW5 (based on number of reads in FAST5 and BLOW5 files)."

if [ ${NUM_BARCODES} -ne ${NUM_BARCODES_PASSED} ]
	then
		echo "STOP AND FIX BLOW5 CONVERSION"
		exit 1
	else
		echo "All set to delete"
	fi

read -p "Delete FAST5 files? (Enter yes or no): " response

case "$response" in
    yes) echo "Ok, deleting"
         rm -rf ${FAST5_DIR};;
    no) echo "Ok, will keep the FAST5s";;
    *) echo "Enter yes or no" ;;
esac
