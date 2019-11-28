# sa
smRNA Aligner Based on Needleman-Wunsch
 
## Introduction  
SA is a software for small RNA sequence align against to reference,   
such as miRNA, piRNA, siRNA and so on. It was based on seed search   
and Needleman-Wunsch global comparison, which makes it more fast   
and accurate for smRNA alignment.   

In order to avoid the seed searching faiture due to the mismatch,   
we take the number of mismatch plus one seed searching strategy.  

## Getting started  
	git clone https://github.com/zhuqianhua/sa.git   
	make   
	cd example   
	sh run.sh   

### Input format
To reduce computation, the input files is Fasta format, thus data  
with Fastq format from smRNA-seq must be transformed to Fasta files  
before running SA, and TFast.pl could convert Fastq format to Fasta 
format.

### Output format
The output files are automatically generated in SA:Z alignment tag  
included in sequence alignment/map format (SAM) format, value 4  
of the second column in which indicated unmapping, while 0 and 16  
indicated positive-strand alignment/forward alignment and negative-strand  
alignment/reverse alignment, respectively. If output results were  
set as the entirety (parameter: -a or --all), all the results were  
presented in SA:Z alignment tag. Particularly, SA also provided statistic  
information of mapping rates for further analysis.

## Note
It relies on multithreaded execution and your system needs to   
support -lpthread. If you have questions about SA, you may send  
the questions to zhuqianhua@bgi.com.  

## Citing
