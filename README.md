# insilicoRD

* Determines the optimal combination of restriction enzymes for the first part of the Reduced Representation Bisulfite Sequencing (RRBS-seq) protocol, in which genomic DNA is digested by restriction enzymes. 
* The user enters the names of the restriction enzymes for which they want to compute cytosine coverage
* This program outputs (1) a bed file containing information (chromosome, starting base pair, ending base pair, and the number of each dinucleotide combination) for every resultant DNA fragment from the restriction enzyme digestion, (2) a txt file containing information on the maximum cytosine count per fragment and the corresponding base pair range, and (3) a txt file with the total cytosine coverage.

## Goals
* Identify the coordinates of restriction enzyme recognition sites in a given sequence 
* Compute cytosine coverage for each restriction enzyme
* Identify the optimal combination of restriction enzymes for the RRBS-seq protocol
* Ultimately, optimize the RRBS-seq protocol 

## Usage

### Required inputs
1. ```enzyme```: the enzyme string. 
2.  ```enzyme2```: a second enzyme string. Default = None. 
3.  ```x```: sequencing read length. Default = 60. 
4.  ```fragment_range_start```: smallest fragment length acceptable. 
5.  ```fragment_range_end```: largest fragment length acceptable. Default = 2000. 
6.  fasta file for each chromosome

### Outputs 
1. BED file with fields: ```chr```, ```start```, ```end```, ```length```, ```CG```, ```CA```, ```CT```, ```CC```, ```TG```, ```AG```, ```GG```
2. maximum_c.txt: with the range with the maximum cytosine count per fragment 
3. {enzyme}_C_coverage.txt: cytosine coverage
4. Fragment Length Density of Restriction Enzymes Plot
