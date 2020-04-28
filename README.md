# sa
**Small RNA Aligner Based on Needleman-Wunsch**
 
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

## Options
Options | Type | Description
--- | --- | ---
-q, --query | str | type of small rna data, f for fasta and q for fastq, default f
-s, --seq | *str | small rna data in fasta/fastq format
-r, --ref | *str | reference sequence in fasta format
-p, --prefix | str | prefix of output, default ./sa
-m, --mismatch | int | number of mismatch, default 1
-l, --seed | int | length of seed, default 9
-t, --thread | int | number of threads, default 1
-a, --all | bool | report all best alignments, default random one
-b, --debug | bool | print the debug information
-h, --help | bool | print this information

*note * indicates required parameters*

## Output
The outputs contain two files, an alignment file in SAM format and  
a statistic file. The alignment results are automatically generated 
in SA:Z alignment tag  included in sequence alignment/map  
format (SAM) format, value 4 of the second column in which indicated  
unmapping, while 0 and 16 indicated positive-strand alignment/forward  
alignment and negative-strand alignment/reverse alignment, respectively.  
If output results were set as the entirety (parameter: -a or --all), all the  
results were presented in SA:Z alignment tag. Particularly, SA also provided  
statistic information of mapping rates for further analysis.

## Note
It relies on multithreaded execution and your system needs to   
support -lpthread. If you have questions about SA, you may send  
the questions to zhuqianhua@bgi.com.  

## Reference
1) Desvignes, T., et al. (2014) Expanding the annotation of zebrafish microRNAs based on small RNA sequencing, Gene, 546, 386-389.
2) Helge, G. and Witold, F. (2008) Molecular biology: the expanding world of small RNAs, Nature, 451, 414.
3) Jing, G., et al. (2014) Comprehensive analysis of human small RNA sequencing data provides insights into expression profiles and miRNA editing, RNA Biology, 11, 1375-1385.
4) Needleman, S.B. and Wunsch, C.D. (1970) A general method applicable to the search for similarities in the amino acid sequence of two proteins, Journal of Molecular Biology, 48, 443-453.
5) Zhang, C. (2009) Novel functions for small RNA molecules, Current Opinion in Molecular Therapeutics, 11, 641-651.
6) Ziemann, M., Kaspi, A. and El-Osta, A. (2016) Evaluation of microRNA alignment techniques, RNA, 22, 1120-1138