#!/bin/bash

################################################################################
#                                                                              #
#    -NEW ChIP-Seq processing pipeline-   by Nel Marin (13/04/2023)            #
#                Based on the script of Nuria Serna                            #    
#                                                                              #
################################################################################

########## REQUIREMENTS

: '
Tool installation & reference genome index
Tools: bowtie2, samtools, deeptools, bedgraphtowig.pl
 '

########## ARGUMENT PARSING

POSITIONAL_ARGS=()

#Set the default values
SAMPLE=""
FOLDER=""
MERGE=""
BINSIZE=1
THREADS=10

while [[ $# -gt 0 ]]; do
  case $1 in
    -s|--samples)
      SAMPLE="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -p|--path)
      FOLDER="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -m|--mergeSamples)
      MERGE="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -b|--binSize)
      BINSIZE="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -t|--threads)
      THREADS="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -h|--help)
      echo " 
Usage: bash ChIPseq_processing.sh [options] -s [files] -p [path]

Required arguments:
  -s | --samples       S  Name of the sample files        
  -p | --path          S  Path to store all of the intermediate files

Optional arguments:  
  -m | --mergeSamples  S  Histones or histone marks to merge (not L1/L2)
  -b | --binSize       N  Size of the bins for bamCompare (default: 1)
  -t | --threads       N  Number of threads/processors (default: 10)
  -h | --help             Help

L1 and L2 files refer to the same sample (same sample ID) but it has been
divided in order to sequence it correctly by the lab in charge. This script 
will merge these files automatically without asking the user for the input. 
However, if different samples (different sample IDs) want to be merged, the
user has to indicate the name of the histone variant or the histone mark that
appears in the filename (e.g. H1_2 or K27me3) using the -m, --mergeSamples
option.
"
      exit 1
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1")  #Save positional arguments
      shift     #Past argument
      ;;
  esac
done

#Check for required variables
if [[ -z "$SAMPLE" ]]
then
  echo "ERROR: -s|--samples is a required argument"
  exit 1
fi 

if [[ -z "$FOLDER" ]]
then
  echo "ERROR: -p|--path is a required argument"
  exit 1
fi

set -- "${POSITIONAL_ARGS[@]}"  #Restore positional parameters



########## 1. ALIGNMENT TO REFERENCE GENOME
#Bowtie2 aligns sequencing reads to long reference sequences

mkdir -p ${FOLDER}/Results

for file in ${FOLDER}/FASTQ_Files/${SAMPLE}  	
do
  bowtie2 -x /home/ajvlab/Softwares/hg19_rDNA/hg19-rDNA_v1.0 \
          -U $file \
	        -S ${file%%.*}.aligned.sam \
          -p $THREADS
  mv ${file%%.*}.aligned.sam ${FOLDER}/Results
  echo -e "\e[1;33m### Alignment of $file to reference genome DONE ###\e[0m"
  echo -e "\e[1;32m### New file: $(basename "$file") ###\e[0m"
done 
 


########## 2. SAM to BAM conversion
#Samtools reads, writes, edits, indexes, views SAM/BAM/CRAM format files

### 2.1. BAM sorting
#Sorting of the alignments and SAM to BAM direct conversion

for file in ${FOLDER}/Results/${SAMPLE}.aligned.sam
do
  samtools sort $file \
             -o ${file%%.*}.sorted.bam 
  echo -e "\e[1;33m### $file from SAM to BAM and sorted DONE ###\e[0m"
  echo -e "\e[1;32m### New file: $(basename "$file") ###\e[0m"
done

#Merge the files from different lanes of a flow cell (L1/L2)
for file in ${FOLDER}/Results/${SAMPLE}.sorted.bam
do
  if [[ $file == *"_L1_"* ]]
  then
    samtools merge "${file/_L1_/_L1L2Merged_}" "${file}" "${file/_L1_/_L2_}"
    mv ${file} ${file}.no
    mv ${file/_L1_/_L2_} ${file/_L1_/_L2_}.no
    echo -e "\e[1;33m### BAM files merged in ${file/_L1_/_L1L2Merged_} DONE ###\e[0m"
    echo -e "\e[1;32m### New file: $(basename "$file") ###\e[0m"
  fi
done 

#Merge the files that share the pattern asked to the user
if [[ $MERGE != "" ]]
then
  firstfile=${FOLDER}/Results/$(ls ${FOLDER}/Results | grep "${MERGE}.*\.sorted\.bam" | head -n 1)
  secondfile=${FOLDER}/Results/$(ls ${FOLDER}/Results | grep "${MERGE}.*\.sorted\.bam" | tail -n 1)
  samtools merge "${firstfile/NoDox_*_T2T/NoDox_Merged_T2T}" \
                 "$firstfile" \
	         "$secondfile"
  mv ${firstfile} ${firstfile}.no
  mv ${secondfile} ${secondfile}.no
  echo -e "\e[1;33m### $MERGE BAM files merged DONE ###\e[0m"
  echo -e "\e[1;32m### New file: $(basename "${firstfile/NoDox_*_T2T/NoDox_Merged_T2T}") ###\e[0m"
fi

#Remove the alignment files (SAM)
rm ${FOLDER}/Results/*.aligned.sam

#Remove the unmerged files
rm ${FOLDER}/Results/*.no


### 2.2. BAM filtering
#Filtering of the alignments with a spcific bit set in the flag field

for file in ${FOLDER}/Results/${SAMPLE}.sorted.bam
do
  samtools view -F 4 $file \
	              -q 4 \
                -o ${file%%.*}.filtered.bam
  echo -e "\e[1;33m### $file filtered DONE ###\e[0m"
  echo -e "\e[1;32m### New file: $(basename "$file") ###\e[0m"
done

### 2.3. BAM indexing
#Indexing of the BAM files for downtream steps

for file in ${FOLDER}/Results/${SAMPLE}.filtered.bam
do
  samtools index $file -b ${file%%.*}.filtered.bai
  echo -e "\e[1;33m### $file indexed DONE ###\e[0m"
  echo -e "\e[1;32m### New file: $(basename "$file") ###\e[0m"
done

#Remove the first binary files (BAM)
rm ${FOLDER}/Results/*.sorted.bam



########## 3. BAM to BEDGraph conversion & input substraction
#DeepTools bamCompare compares two BAM files while normalizing

conda init bash
source activate Python37
for file in ${FOLDER}/Results/${SAMPLE}.filtered.bam
do
  #Change this if (not working)
  if [[ $file == "Input_" ]]
  then
    continue
  fi
  bamCompare -b1 $file \
             -b2 ${FOLDER}/Results/Input_*.filtered.bam \
             -of bedgraph \
             -o  ${file%%.*}.bdg \
  	         --operation subtract \
  	         --normalizeUsing CPM \
             --scaleFactorsMethod None \
             --binSize $BINSIZE \
             -p $THREADS
	     #binSize to reduce the bin size and have more resolution (default 50)
	     #p to use more threads (default 2)
  echo -e  "\e[1;33m### $file from BAM to BEDGraph conversion & input substraction DONE ###\e[0m"
  echo -e "\e[1;32m### New file: $(basename "$file").bdg ###\e[0m"
done
conda deactivate



########## 4. BEDGraph to WIG conversion
#It converts bedGraph to fixedStep wig files (defined step size)

for file in ${FOLDER}/Results/${SAMPLE}.bdg
do 
  /home/ajvlab/Softwares/bedgraph_to_wig.pl --bedgraph $file \
                                            --wig ${file%%.*}.wig \
                              					    --step 50
  echo -e "\e[1;33m### $file BEDGraph to WIG conversion DONE ###\e[0m"
  echo -e "\e[1;32m### New file: $(basename "$file") ###\e[0m"
done
