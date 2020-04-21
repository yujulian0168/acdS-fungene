# Making a fungene resource
## acdS example

### Download seed sequences
- acdS fungene repository 
- http://fungene.cme.msu.edu/hmm_detail.spr?hmm_id=564

```
minimum score of 425
minimum aa size of 338
minimum HMM coverage of 95% 
do not display the environmental sequences
#parse down to 4,310 sequences
```
```
Dereplicate using RDPtools Dereplicate
# 1,059 unique reference sequences
```

### Alignment
```
Input seed sequences (AA) into clustal omega for alignment, can upload a fasta file. Choose protein for type of sequences, and output format: Stockholm
```

### Build hmm model
	
- Move to agave fungene resources


```
scp acds.stockholm julianyu@agave.asu.edu:/home/julianyu/fungene/fungene_pipeline/resources/acdS/

module load hmmer/3.1b2

hmmbuild acdS.hmm acds.stockholm

#open .hmm file to confirm it was written properly. 


To make framebot.fasta file, move a copy of the fasta file containing seed sequences to resources\acdS 
```

>If you would prefer to have framebot.fasta be an alignment then you can redo the seed sequence alignment with the same
sequences and instead select PEARSON/FASTA output. (AFS is unsure why there is a difference in some of the framebot.fasta files 
in the different functional gene resource folders). 

