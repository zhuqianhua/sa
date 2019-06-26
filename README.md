# sa
smRNA Aligner Based on Needleman-Wunsch

## Getting started  
	git clone https://github.com/zhuqianhua/sa.git  
	make  
	cd example  
	sh run.sh  
  
## Introduction  
SA is a software for small RNA sequence align against to reference,   
such as miRNA, piRNA, siRNA and so on. It was based on seed search   
and Needleman-Wunsch global comparison, which makes it more fast   
and accurate for smRNA alignment.   

In order to avoid the seed searching faiture due to the mismatch,   
we take the number of mismatch plus one seed searching strategy.  

## Note
It relies on multithreaded execution and your system needs to   
support -lpthread. If you have questions about SA, you may send  
the questions to zhuqianhua@bgi.com.  

## Citing
