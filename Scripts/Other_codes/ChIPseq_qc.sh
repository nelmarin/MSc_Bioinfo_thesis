#!/bin/bash

mkdir /home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/Quality_control


###### FASTQC
for file in /home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/r2/FASTQ_Files/* 
do
  fastqc $file -o /home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/Quality_control 
  echo "### Quality control of $(basename "$file") DONE ###"
done


for file in /home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/r3/FASTQ_Files/*
do
  fastqc $file -o /home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/Quality_control
  echo "### Quality control of $(basename "$file") DONE ###"
done


for file in /home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/PTMs/FASTQ_Files/*
do
  fastqc $file -o /home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/Quality_control
  echo "### Quality control of $(basename "$file") DONE ###"
done


###### MULTIQC
multiqc /home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/Quality_control/. \
	-o /home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/Quality_control/multiqc
echo "### Final quality control DONE ###"


