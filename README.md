# MobiVision-M

## Introduction
MobiVision-M is a linux based software design specifically for single-microbe RNA sequencing analysis.
Start from raw data, MobiVision-M preforms barcoding, adaptor-removing, alignment, re-assignment of multi-loci mapped reads, statistical analysis of data and making a summary report.

## Hardware
- Processor: An 8-core Intel or AMD processor, x86 architecture (16 cores or more is recommended)  
- Memory: 32GB (64GB or more is recommended)  
- Operating System: Linux operating system, such as 64-bit CentOS 7, ubuntu:22.04, or higher versions are recommended 

## Environment
MobiVision-M was tested in an environment where the following software is available in PATH:
1) STAR (2.7.10b)
2) samtools (1.12)
3) fastqc (v0.12.1)
4) cutadapt (3.5)
5) fastp (0.23.4)
6) python (3.8.12)

MobiVision-M was tested in an environment where the following packages installed:
1) pandas	1.4.2
2) numpy	1.24.4
3) scanpy	1.9.8
4) rtoml	0.9.0
5) structlog	25.3.0
6) loguru	0.7.3
7) mako	1.3.0
8) matplotlib	3.6.0
9) sklearn	0.24.2
10) scipy	1.10.1

## Quick start
### Making an index
```bash
mobivision-M mkindex \
-n CP001363 \
-f /share/home/sc/Projects/microbeRNA-seq/reference_20231128_CP001363/CP001363.fasta \
-g /share/home/sc/Projects/microbeRNA-seq/reference_20231128_CP001363/CP001363.gtf \
-o /share2/Data/sc/Mobi-RNA_seq_project/final_test/A-1/CP001363
```
### Quantify a data
```bash
o_dir=/share/home/sc/Projects/Mobi-RNA_seq_project/final_test
raw_data_dir=/share/home/sc/Projects/Mobi-RNA_seq_project/final_test/raw_datas/test_data
ID=C-1
mobivision-M quantify \
-f $raw_data_dir/240131SW_240129B-S-SW-E01 \
-t 12 \
-i /share/home/sc/Projects/microbeRNA-seq/reference_20240112_GCA_003019295/CP028101 \
-s 240131SW_240129B-S-SW-E01 \
-o /share/home/sc/Projects/Mobi-RNA_seq_project/final_test/$ID/240131SW_240129B-S-SW-E01
```
### Re-call microbes form an analysised result
```bash
o_dir=/share/home/sc/Projects/Mobi-RNA_seq_project/final_test
ID=C-1
mkdir /share/home/sc/Projects/Mobi-RNA_seq_project/final_test/$ID
mobivision-M rcmicrobe \
-i /share/home/sc/Projects/Mobi-RNA_seq_project/final_test/B-1/240131SW_240129B-S-SW-E01/240131SW_240129B-S-SW-E01 \
-o $o_dir/$ID/$ID're_call'\
-t 8 \
--cr2.2
```
## Detailed docs
- [Introduction](docs/MobiVision-M_Introduction.md)
- [User manual](docs/MobiVision-M_v1.3_User_Manual.md)
- [Analysis custom library](docs/MobiVision-M_Custom_Library_Structure_Analysis_Pipeline.md)


## Author/Support
Shan Chao, 

## Limitations/License

## Funding