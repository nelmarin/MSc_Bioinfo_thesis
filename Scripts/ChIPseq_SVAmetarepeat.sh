#!/bin/bash

################################################################################
#                                                                              #
#            -MetaRepeat analysis-   by Nel Marin (20/06/2023)                 #
#                                                                              #
################################################################################

########## REQUIREMENTS

: '
Tool installation
Tools: deeptools
 '

########## ARGUMENT PARSING

POSITIONAL_ARGS=()

#Set the default values
FILE=""
FOLDER=""
REPEATS=""
AFTER=3000
BEFORE=3000
BINSIZE=1

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
	-a|--afterEndRegion)
      AFTER="$2"
      shift     #Past argument
      shift     #Past value
      ;;
	-b|--beforeStartRegion)
      BEFORE="$2"
      shift     #Past argument
      shift     #Past value
      ;;
	-s|--binsize)
      BINSIZE="$2"
      shift     #Past argument
      shift     #Past value
      ;;
    -h|--help)
      echo " 
Usage: bash ChIPseq_metarepeat.sh -f [files] -p [path] -r [repeats]

Required arguments:
  -f | --files         S  Name of the files        
  -p | --path          S  Path to store all of the intermediate files
  -r | --repeats       S  Name of the repeats file (BED/GTF format)

Optional arguments:
  -b | --beforeStartRegion     N  Distance upstream of the regions start site (default: 3000)
  -a | --afterEndRegion        N  Distance downstream of the regions end site (default: 3000)
  -s | --binsize               N  Length in bases for score computing (default: 1)
  -h | --help                     Help

The optional arguments are applied in the computeMatrix function.
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



##### 1. Remove the entries corresponding to the mitochondral chomosome
#Create a new folder
mkdir -p ${FOLDER}/Heatmaps
#Remove chrM entries
for file in ${FOLDER}/Results_v0/${FILE}
do
  awk '$1 != "chrM"' $file > "$file"_noM
  echo "### Chromosome M removed from $file DONE ###"
done


##### 2. Convert BedGraph files to BigWig ones
for file in ${FOLDER}/Results_v0/*.bdg_noM
do
  filename=$(basename "$file")
  output_file="${FOLDER}/Heatmaps/${filename%.bdg}.bw"
  /home/ajvlab/Softwares/KentUtilities_UCSC/bedGraphToBigWig $file /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_chrom.sizes $output_file
  echo "### BedGraph to BigWig conversion of $filename DONE ###"
done


##### 3. Get only the SVAs from the GTF file
awk -F '\t' '$NF ~ /Retroposon/'${REPEATS} > /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon.gtf
count_TE=$(wc -l ${REPEATS} | cut -d ' ' -f 1)
count_SVA=$(wc -l /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon.gtf | cut -d ' ' -f 1)
echo "Summary:"
echo "The number of repeats is $count_TE"
echo "The number of Retroposon repeats (SVAs) is $count_SVA"


##### 4. Separate depending on the gene_id
for X in A B C D E F
do 
  grep -w 'gene_id "SVA_'"$X"'"' /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon.gtf > /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_"$X".gtf
done
count_A=$(wc -l /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_A.gtf | cut -d ' ' -f 1)
count_B=$(wc -l /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_B.gtf | cut -d ' ' -f 1)
count_C=$(wc -l /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_C.gtf | cut -d ' ' -f 1)
count_D=$(wc -l /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_D.gtf | cut -d ' ' -f 1)
count_E=$(wc -l /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_E.gtf | cut -d ' ' -f 1)
count_F=$(wc -l /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_F.gtf | cut -d ' ' -f 1)
count_notA=$(awk -F '\t' '$NF ~ /SVA_[BCDEF]/' /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_A.gtf | wc -l)
count_notB=$(awk -F '\t' '$NF ~ /SVA_[ACDEF]/' /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_B.gtf | wc -l)
count_notC=$(awk -F '\t' '$NF ~ /SVA_[ABDEF]/' /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_C.gtf | wc -l)
count_notD=$(awk -F '\t' '$NF ~ /SVA_[ABCEF]/' /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_D.gtf | wc -l)
count_notE=$(awk -F '\t' '$NF ~ /SVA_[ABCDF]/' /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_E.gtf | wc -l)
count_notF=$(awk -F '\t' '$NF ~ /SVA_[ABCDE]/' /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_F.gtf | wc -l)
echo "Summary:"
echo "The number of SVAs A: $count_A ($count_notA)"
echo "The number of SVAs B: $count_B ($count_notB)"
echo "The number of SVAs C: $count_C ($count_notC)"
echo "The number of SVAs D: $count_D ($count_notD)"
echo "The number of SVAs E: $count_E ($count_notE)"
echo "The number of SVAs F: $count_F ($count_notF)"


##### 5. Modify the duplicated lines
for file in /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_?.gtf
do
  awk -F'\t' '{
    if (!seen[$NF]) {
      seen[$NF] = 1;
    } else {
      seen[$NF]++;
    }
    count = seen[$NF];
    gsub(/"; family_id/, "_" count "&");
    print;
  }' "$file" > "${file%.*}_ready.gtf"
done


##### 6. Use computeMatrix
#Hide the two K9me3 samples (not used)
#mkdir ${FOLDER}/Heatmaps/Hiden
#mv ${FOLDER}/Heatmaps/ChIP_K9me3_NoDox_AJV7_T47D6c1_L1L2Merged_T2T.bdg_noM.bw ${FOLDER}/Heatmaps/Hiden/ChIP_K9me3_NoDox_AJV7_T47D6c1_L1L2Merged_T2T.bdg_noM.bw
#mv ${FOLDER}/Heatmaps/ChIP_K9me3_NoDox_AJV84X_T47D6c4_T2T.bdg_noM.bw ${FOLDER}/Heatmaps/Hiden/ChIP_K9me3_NoDox_AJV84X_T47D6c4_T2T.bdg_noM.bw

#Use computeMatrix
conda init bash
source activate Python37

computeMatrix scale-regions \
        -S ${FOLDER}/Heatmaps/*.bw \
        -R /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_A_ready.gtf /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_B_ready.gtf /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_C_ready.gtf /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_D_ready.gtf /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_E_ready.gtf /home/ajvlab/Softwares/REFs_T2TCHM13v2.0_Files/T2TCHM13v2.0_rmsk_TE_Retroposon_F_ready.gtf \
        -a ${AFTER} \
        -b ${BEFORE} \
        --skipZeros \
        --binSize ${BINSIZE} \
        --transcriptID exon \
        --transcript_id_designator transcript_id \
        -o ${FOLDER}/Heatmaps/ComputedMatrix_v1.gz



##### 7. Use plotHeatmap
plotHeatmap -m ${FOLDER}/Heatmaps/ComputedMatrix_v1.gz \
            -z "SVA A" "SVA B" "SVA C" "SVA D" "SVA E" "SVA F" \
	        --legendLocation upper-right \
	        --colorList blue,white,red \
            --startLabel "" \
	        --endLabel "" \
	        --samplesLabel "H1.0" "H1.2" "H1.3" "H1.4" "H1.5" "H1X" "K9me3" "K27me3" \
	        -x "Meta-repeat" \
	        -y "Repeats" \
	        -o ${FOLDER}/Heatmaps/MetaRepeat_v1.png
conda deactivate

