#!/bin/bash

################################################################################
#                                                                              #
#           -ChIP-Seq mapping pipeline-   by Nel Marin (24/04/2023)            #
#                 Based on the pipeline of Nuria Serna                         #    
#                                                                              #
################################################################################

########## REQUIREMENTS

: '
Tool installation
Tools: bedtools
 '

########## ARGUMENT PARSING

POSITIONAL_ARGS=()

#Set the default values
SAMPLE=""
FOLDER=""
REFERENCE=""
OUTPUT=""

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
    -r|--reference)
      REFERENCE="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -o|--output)
      OUTPUT="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -h|--help)
      echo " 
Usage: bash ChIPseq_mapping.sh -s [files] -p [path] -r [reference] -o [output]

Required arguments:
  -s | --samples       S  Name of the sample files        
  -p | --path          S  Path to store all of the intermediate files
  -r | --reference     S  Name of the file storing genome features
  -o | --output        S  Additional name to the output

Optional arguments:  
  -h | --help             Help

In order to select all the samples in a folder, it has to be specified as \"*\" and not * only.
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

if [[ -z "$REFERENCE" ]]
then
  echo "ERROR: -r|--reference is a required argument"
  exit 1
fi

if [[ -z "$OUTPUT" ]]
then
  echo "ERROR: -o|--output is a required argument"
  exit 1
fi

set -- "${POSITIONAL_ARGS[@]}"  #Restore positional parameters


########## CHROMOSOME REORDERING (27/03/2023)
#Chromosome lexicographically reordering 

for file in ${FOLDER}/Results/${SAMPLE}.bdg
do
  if [[ $(cut -f1 $file | uniq | head -n 2 | tail -n 1) == "chr2" ]] 
  then	  
    sort -k1,1 -k2,2n $file > temp
    mv temp ${FOLDER}/Results/"$(basename "$file")"
    echo -e "\e[1;33m### $file chromosome reordering DONE ###\e[0m"
    echo -e "\e[1;32m### New file: $(basename "$file") ###\e[0m"
  fi
    echo -e "\e[1;34m### $file is already reordered ###\e[0m"
done

: '
It loops through all the files in the current directory with the extension
".bdg" and sorts the chromosomes lexicographically. The sorted content is 
written to a temporary file called "temp". Finally, the temporary file is 
renamed to the original file name and overwritten with the sorted content.
 '



########## CHIP-SEQ BEDGRAPH MAPPING (22/03/2023)
#Mapping of the values of one file to the qualitative features of another one

for file in ${FOLDER}/Results/${SAMPLE}.bdg
do
  bedtools map -a ${REFERENCE} \
               -b $file \
   	           -c 4 \
  	           -o mean \
  	           > ${file%%.*}_vs${OUTPUT}.txt
  echo -e "\e[1;33m### $file mapped against the cytobands DONE ###\e[0m"
  echo -e "\e[1;32m### New file: ${SAMPLE}_vs${OUTPUT}.txt ###\e[0m"
done

: '
The output is basically the same as the A file with a new
column where the mean calculated is added.
 '


