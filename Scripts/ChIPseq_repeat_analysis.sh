#!/bin/bash

################################################################################
#                                                                              #
#       -Abundance analysis in the repeats-   by Nel Marin (03/05/2023)        #
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
FILE=""
FOLDER=""
REPEATS=""

while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--files)
      FILE="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -p|--path)
      FOLDER="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -r|--repeats)
      REPEATS="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -h|--help)
      echo " 
Usage: bash ChIPseq_repeat_analysis.sh -f [files] -p [path] -r [repeats]

Required arguments:
  -f | --files         S  Name of the files        
  -p | --path          S  Path to store all of the intermediate files
  -r | --repeats       S  Name of the repeats file

Optional arguments:  
  -h | --help             Help
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
if [[ -z "$FILE" ]]
then
  echo "ERROR: -f|--files is a required argument"
  exit 1
fi 

if [[ -z "$FOLDER" ]]
then
  echo "ERROR: -p|--path is a required argument"
  exit 1
fi

if [[ -z "$REPEATS" ]]
then
  echo "ERROR: -r|--repeats is a required argument"
  exit 1
fi

set -- "${POSITIONAL_ARGS[@]}"  #Restore positional parameters




########## 1. Dataframe construction

#Read the data and create a dataframe with all the samples.
#The repeats would be indicated in the rows; while the samples, in the columns.
#The columns should be: shared columns, class, family, group, sample data (x10).



### 1.1. Shared columns

#Create a new folder for this analysis
mkdir -p ${FOLDER}/Analysis

#Extract the shared columns
for file in ${FOLDER}/Results/${FILE}
do
  cut -f 1-8 $file > ${FOLDER}/Analysis/shared.txt
  break
done
echo -e "\e[1;33m### Shared columns extracted DONE ###\e[0m"
echo -e "\e[1;32m### New file: ${FOLDER}/Analysis/shared.txt ###\e[0m"

#Check if the files are the same (cmp does not produce an output if true)
cut -f 1-8 ${FOLDER}/Results/ChIP_H1_2_NoDox_AJV26F_T47D6c2_T2T_vsTEs.txt \
    > ${FOLDER}/Analysis/test.txt
cmp ${FOLDER}/Analysis/shared.txt ${FOLDER}/Analysis/test.txt



### 1.2. Taxonomy columns

#The F option is used in order to allow awk to read tab-delimited files
awk -F '\t' '{
    match($9, /gene_id "([^"]+)"; transcript_id "([^"]+)"; \
               family_id "([^"]+)"; class_id "([^"]+)"/, ids);
    print ids[1] > "gene_id.txt";
    print ids[2] > "transcript_id.txt";
    print ids[3] > "family_id.txt";
    print ids[4] > "class_id.txt"}' ${REPEATS}
echo -e "\e[1;33m### Taxonomy columns extracted DONE ###\e[0m"

#Move the created files to the proper folder
mv *_id.txt ./ChIP_T47D_6c1_6c2_6c4_6c6_6c7/Analysis
echo -e "\e[1;32m### New files: ${FOLDER}/Analysis/gene_id.txt \
	                        ${FOLDER}/Analysis/transcript_id.txt \
	                        ${FOLDER}/Analysis/family_id.txt \
	                        ${FOLDER}/Analysis/class_id.txt ###\e[0m"

#Paste the files to the shared.txt file (not transcript_id.txt)
paste -d '\t' ${FOLDER}/Analysis/shared.txt \
              ${FOLDER}/Analysis/class_id.txt \
              ${FOLDER}/Analysis/family_id.txt \
              ${FOLDER}/Analysis/gene_id.txt \
              > ${FOLDER}/Analysis/merged.txt
echo -e "\e[1;33m### Taxonomy and shared columns merged DONE ###\e[0m"
echo -e "\e[1;32m### New file: ${FOLDER}/Analysis/merged.txt ###\e[0m"



### 1.3. Data columns 

#Extract the last column of each data file (histone & histone marks)
for file in ${FOLDER}/Results/${FILE}
do
  awk '{print $NF}' $file \
      > ${FOLDER}/Analysis/$(basename "$file" .txt)_last_column.txt
  echo -e "\e[1;33m### Data columns extracted DONE ###\e[0m"
  echo -e "\e[1;32m### New file: \ 
          ${FOLDER}/Analysis/$(basename "$file" .txt)_last_column.txt ###\e[0m"
done

#Paste all the data files to the merged.txt file 
paste -d '\t' ${FOLDER}/Analysis/merged.txt \
              ${FOLDER}/Analysis/*_last_column.txt \
              > ${FOLDER}/Analysis/thedataframe.tsv
echo -e "\e[1;33m### merged.txt and data columns merged DONE ###\e[0m"
echo -e "\e[1;32m### New file: ${FOLDER}/Analysis/thedataframe.tsv ###\e[0m"




########## 2. Blacklist regions removal

#Remove the repeats that overlap with the blacklist regions.
#To do this, the operation "count" of bedtools map must be used
#It counts the number of overlaps between two files (interested in 0s)

#Chromosome reordering
sort -k1,1 -k2,2n ${FOLDER}/Analysis/T2T.excluderanges_from_R_Bioconductor.bed \
     > ${FOLDER}/Analysis/T2T.excluderanges_from_R_Bioconductor_sorted.bed
echo -e "\e[1;33m### Chromosome reordered DONE ###\e[0m"
echo -e "\e[1;32m### New file: \
       	${FOLDER}/Analysis/T2T.excluderanges_from_R_Bioconductor_sorted.bed ###\e[0m"


#Mapping
#I could not do it directly to the dataframe file because there were issues
#with the tabs I assume
bedtools map \
     -a ${REPEATS} \
     -b ${FOLDER}/Analysis/T2T.excluderanges_from_R_Bioconductor_sorted.bed \
     -c 4 \
     -o count \
     > ${FOLDER}/Analysis/repeat_count_overlapping.txt
echo -e "\e[1;33m### Mapping DONE ###\e[0m"

#Extract the last column with the overlapping indexes
awk -i inplace '{print $NF}' FS='\t' ${FOLDER}/Analysis/repeat_count_overlapping.txt
echo -e "\e[1;32m### New file: \
        ${FOLDER}/Analysis/repeat_count_overlapping.txt ###\e[0m"

#Paste the 0/1 column to the dataframe
paste -d '\t' ${FOLDER}/Analysis/thedataframe.tsv \
              ${FOLDER}/Analysis/repeat_count_overlapping.txt \
              > ${FOLDER}/Analysis/thedataframe_blacklisted.tsv
echo -e "\e[1;33m### Blacklist index added to the dataframe DONE ###\e[0m"
echo -e "\e[1;32m### New file: \
       	${FOLDER}/Analysis/thedataframe_blacklisted.tsv ###\e[0m"

#Extract the repeats with 0 in the last column 
awk -F '\t' '$NF == 0 && NF == 22' \
	${FOLDER}/Analysis/thedataframe_blacklisted.tsv \
	> ${FOLDER}/Analysis/thedataframe_filtered.tsv
echo -e "\e[1;33m### Blacklist regions removal DONE ###\e[0m"
echo -e "\e[1;32m### New file: \
        ${FOLDER}/Analysis/thedataframe_filtered.tsv ###\e[0m"




########## 3. No abundance removal

#Remove the repeats that bedtools map was not able to calculate the mean ('.')
#I created "ChIPseq_dots_check.sh" to check if there are dots in the data
#All 0 columns removed using R




########## 4. Median abundance calculation

#Remove those repeats with "DNA?" in the class level
awk -F '\t' '$9 != "DNA?" {print}' \
       	${FOLDER}/Analysis/thedataframe_filtered.tsv \
       	> ${FOLDER}/Analysis/thedataframe_noDNA?.tsv
echo -e "\e[1;33m### DNA? class removal DONE ###\e[0m"
echo -e "\e[1;32m### New file: \
        ${FOLDER}/Analysis/thedataframe_complete.tsv ###\e[0m"

#Remove the last column (all 0)
cut -d $'\t' -f 1-21 /home/ajvlab/Nel/ChIP_T47D_6c1_6c2_6c4_6c6_6c7/Analysis/thedataframe_noDNA?.tsv > /home/ajvlab/Nel/ChIP_T47D_6c1_6c2_6c4_6c6_6c7/Analysis/thedataframe_complete.tsv
