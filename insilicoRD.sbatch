#!/bin/bash 
#SBATCH -A mgowda
#SBATCH -p standard
#SBATCH --job-name=insilicoRD
#SBATCH --cpus-per-task=20
#SBATCH --time=32:00:00

main(){
module load gcc/7.1.0
module load bedtools/2.26.0
module load gcc/9.2.0 mvapich2/2.3.3
module load anaconda/2019.10-py3.7

python get_coordinates.py
# output = digested_COORD.bed

module load gcc/7.1.0
module load bedtools/2.26.0

# merge features
bedtools merge -d -1 -i features2.bed > merged_features.bed
# remove header of outputs 
sed -i '1d' digested_COORD.bed 
# sort outputs 
sort -k 1,1 -k2,2n digested_COORD.bed > sorted_digested_COORD.bed
# intersect - find overlap
bedtools intersect -wa -a sorted_digested_COORD.bed -b merged_features.bed > intersected.bed
# merge again 
bedtools merge -c 4,5,6,7 -o distinct -d -1 -i intersected.bed > merged_intersected.bed

python cytosines_covered.py # outputs = intersected_COV.bed, C_coverage.txt
python sw.py # outputs = sw.csv, maximum_c.txt
}

main
exit
