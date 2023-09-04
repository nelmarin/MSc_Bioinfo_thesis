#!/bin/bash

################################################################################
#                                                                              #
#              -ChIPseq data check-   by Nel Marín (08/05/2023)                #
#                                                                              #
################################################################################

########## ARGUMENT PARSING

POSITIONAL_ARGS=()

#Set the default values
FILE=""

while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--files)
      FILE="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -h|--help)
      echo " 
Usage: bash ChIPseq_dots_check.sh -f [file] -p [path]

Required arguments:
  -f | --file         S  Name of the file        

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
  echo "ERROR: -f|--file is a required argument"
  exit 1
fi 

set -- "${POSITIONAL_ARGS[@]}"  #Restore positional parameters



########## Check for empty data

#Count the number of columns in the input file
n_columns=$(awk -F '\t' '{print NF; exit}' $FILE)

#Count the number of rows in the input file
n_rows=$(wc -l $FILE | awk '{print $1}')

#Summary of the file

if [ $n_columns == 1 ]
then
  echo "Input file has $n_columns column"
  echo "and $n_rows rows."
  pattern=$(grep '^\.' $FILE)
  if [[ -z "$pattern" ]]
  then
    echo "There are NO dots (\.) in the input file"
  else 
    echo "There are dots (\.) in the input file"
  fi
else
  echo "Input file has $n_columns columns"
  echo "and $n_rows rows."
  pattern=$(awk '{print $12,$13,$14,$15,$16,$17,$18,$19,$20,$21}' $FILE | grep '^\.')
  if [[ -z "$pattern" ]]
  then
    echo "There are NO dots (\.) in the input file"
  else 
    echo "There are dots (\.) in the input file"
  fi
fi  


