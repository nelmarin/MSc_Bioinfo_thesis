#!/bin/bash

################################################################################
#                                                                              #
#                -rDNA analysis-   by Nel Marin (13/07/2023)                   #
#                                                                              #
################################################################################

########## REQUIREMENTS

: '
Tool installation
Tools: deeptools
 '

##### 1. Extract the chromosome R lines in other bedgraph files
#Step not needed if working with bigwig samples 
mkdir -p /home/ajvlab/Nel/ChIP_T47D_6c1_6c2_6c4_6c6_6c7/Profile_rDNA_hg19/Profile
for file in /home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/Profile/*.bdg
do
  awk '$1 == "chrR"' $file > "$file"_chrR
  echo "### Other than chromosome R removed from $file DONE ###"
done


##### 2. Create the chrom.sizes file (only chrR) 
#Step not needed if working with bigwig samples 
echo -e "chrR\t44838" > /home/ajvlab/Softwares/REFs_hg19_rDNA/hg19.chrom.size.chrR


##### 3. Convert BedGraph files to BigWig ones
#Step not needed if working with bigwig samples 
for file in /home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/Profile/*.bdg_chrR
do
  filename=$(basename "$file")
  output_file="/home/ajvlab/Nel/ChIP_T47D_Samples_for_rDNA/Profile/${filename%.bdg}.bw"
  /home/ajvlab/Softwares/KentUtilities_UCSC/bedGraphToBigWig $file /home/ajvlab/Softwares/REFs_hg19_rDNA/hg19.chrom.size.chrR $output_file
  echo "### BedGraph to BigWig conversion of $filename DONE ###"
done


##### 4. Extract the chromosome R lines in other BED file
awk '$1 == "chrR"' /home/ajvlab/Softwares/REFs_hg19_rDNA/hg19-rDNA_v1.0.bed > /home/ajvlab/Softwares/REFs_hg19_rDNA/hg19-rDNA_v1.0_chrR.bed
##### Additional. Create a BED file
echo -e "chrR\t1\t44838" > /home/ajvlab/Softwares/REFs_hg19_rDNA/hg19-rDNA_complete_chrR.bed


##### 5. Use computeMatrix
conda init bash
source activate Python37
computeMatrix scale-regions \
        -S /home/ajvlab/Nel/ChIP_T47D_6c1_6c2_6c4_6c6_6c7/Profile_rDNA_hg19/First_profile/ChIP_H1*.bw \
        -R /home/ajvlab/Softwares/REFs_hg19_rDNA/hg19-rDNA_complete_chrR.bed \
        --skipZeros \
        --binSize 1 \
	-p 10 \
        -o /home/ajvlab/Nel/ChIP_T47D_6c1_6c2_6c4_6c6_6c7/Profile_rDNA_hg19/Profile/ComputedMatrix_rDNA_v8.gz


##### 6. Use plotProfile
plotProfile -m /home/ajvlab/Nel/ChIP_T47D_6c1_6c2_6c4_6c6_6c7/Profile_rDNA_hg19/Profile/ComputedMatrix_rDNA_v8.gz \
            -o /home/ajvlab/Nel/ChIP_T47D_6c1_6c2_6c4_6c6_6c7/Profile_rDNA_hg19/Profile/chrR_hg19-rDNA_profile_6c_H1_variants.png \
	    --perGroup \
            --startLabel "0" \
	    --endLabel "44838 nt" \
	    --legendLocation "upper-right" \
	    --samplesLabel "H1.0" "H1.2" "H1.3" "H1.4" "H1.5" "H1X"
conda deactivate


