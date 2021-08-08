#!/bin/bash 
#SBATCH -A mgowda
#SBATCH -p standard
#SBATCH --job-name=insilicoRD
#SBATCH --cpus-per-task=20
#SBATCH --time=15:00:00

main(){
module load gcc/7.1.0 openmpi/3.1.4
module load anaconda/2019.10-py3.7
module load bedtools/2.26.0
pip install --user pybedtools



ENZ="MspI_65"
python insilicoRD.py
python sw.py


#to compute number of cytosine of interest are covered in each sliding window fragment range
module load gcc/7.1.0
module load bedtools/2.26.0

sed -i '1d' sw.csv
>temp.csv
for START in $(cut -f 1 sw.csv)
do
awk -v S=$START '{ FS="\t"; OFS="\t"; if ($7 >= S && $7 < S+200) print $1,$2,$3,$4,$5,$6,$7;}' intersected_COV.bed > "Filtered_fragments.bed"
sort -k 1,1 -k2,2n Filtered_fragments.bed > Filtered_fragments_srtd.bed
bedtools intersect -sorted -wa -a Filtered_fragments_srtd.bed -b merged_features.bed | wc -l >> temp.csv	
done	
echo -e 'Start\tEnd\tC_per_read\tFeatures_of_interest\tTotal_C\tFragments' > $ENZ"_output.txt"
paste --delimiters="\t" temp.csv sw.csv | awk '{OFS="\t"; print $2,$3,$4,$1,$5,$6}' >> $ENZ"_output.txt"
rm temp.csv digested_COORD.bed Filtered_fragments.bed Filtered_fragments_srtd.bed	

#R insilicoRD_plots.R

}

main
exit
