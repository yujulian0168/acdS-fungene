#Fungene Pipeline
##*acdS* Amplicon Workflow 

###On Agave server

scp -r ~/Downloads/Raw/ julianyu@agave.asu.edu:/home/julianyu/

> Count reads in fastq files

```
for i in *_L001_R1_001.fastq;
do expr $(cat $i|wc -l) / 4;
done

for i in *_L001_R2_001.fastq;
do expr $(cat $i|wc -l) / 4;
done
```

- Merge paired-end reads

```
interactive -n 20
cd /home/julianyu/Raw/fastq/
module load casper/0.8.2

casper 1B_GGAGCTAC-TCGACTAG_L001_R1_001.fastq 1B_GGAGCTAC-TCGACTAG_L001_R2_001.fastq -o out_1B -j
casper 1C_GCGTAGTA-TCGACTAG_L001_R1_001.fastq 1C_GCGTAGTA-TCGACTAG_L001_R2_001.fastq -o out_1C -j
casper 1D_CGGAGCCT-TCGACTAG_L001_R1_001.fastq 1D_CGGAGCCT-TCGACTAG_L001_R2_001.fastq -o out_1D -j
casper 1E_TACGCTGC-TCGACTAG_L001_R1_001.fastq 1E_TACGCTGC-TCGACTAG_L001_R2_001.fastq -o out_1E -j
casper 1F_ATGCGCAG-TCGACTAG_L001_R1_001.fastq 1F_ATGCGCAG-TCGACTAG_L001_R2_001.fastq -o out_1F -j
casper 2A_TAGCGCTC-TCGACTAG_L001_R1_001.fastq 2A_TAGCGCTC-TCGACTAG_L001_R2_001.fastq -o out_2A -j
casper 2B_ACTGAGCG-TCGACTAG_L001_R1_001.fastq 2B_ACTGAGCG-TCGACTAG_L001_R2_001.fastq -o out_2B -j
casper 2C_CCTAAGAC-TCGACTAG_L001_R1_001.fastq 2C_CCTAAGAC-TCGACTAG_L001_R2_001.fastq -o out_2C -j
casper 2D_CGATCAGT-TCGACTAG_L001_R1_001.fastq 2D_CGATCAGT-TCGACTAG_L001_R2_001.fastq -o out_2D -j
casper 2E_TGCAGCTA-TCGACTAG_L001_R1_001.fastq 2E_TGCAGCTA-TCGACTAG_L001_R2_001.fastq -o out_2E -j
casper 2F_TCGACGTC-TCGACTAG_L001_R1_001.fastq 2F_TCGACGTC-TCGACTAG_L001_R2_001.fastq -o out_2F -j
casper 3A_ACTCGCTA-TTCTAGCT_L001_R1_001.fastq 3A_ACTCGCTA-TTCTAGCT_L001_R2_001.fastq -o out_3A -j
casper 3B_GGAGCTAC-TTCTAGCT_L001_R1_001.fastq 3B_GGAGCTAC-TTCTAGCT_L001_R2_001.fastq -o out_3B -j
casper 3C_GCGTAGTA-TTCTAGCT_L001_R1_001.fastq 3C_GCGTAGTA-TTCTAGCT_L001_R2_001.fastq -o out_3C -j
casper 3D_CGGAGCCT-TTCTAGCT_L001_R1_001.fastq 3D_CGGAGCCT-TTCTAGCT_L001_R2_001.fastq -o out_3D -j
casper 3E_TACGCTGC-TTCTAGCT_L001_R1_001.fastq 3E_TACGCTGC-TTCTAGCT_L001_R2_001.fastq -o out_3E -j
casper 3F_ATGCGCAG-TTCTAGCT_L001_R1_001.fastq 3F_ATGCGCAG-TTCTAGCT_L001_R2_001.fastq -o out_3F -j
casper BRR1_TGCAGCTA-CCTAGAGT_L001_R1_001.fastq BRR1_TGCAGCTA-CCTAGAGT_L001_R2_001.fastq -o out_BRR1 -j
casper BRR2_TCGACGTC-CCTAGAGT_L001_R1_001.fastq BRR2_TCGACGTC-CCTAGAGT_L001_R2_001.fastq -o out_BRR2 -j
casper BRR3_ACTCGCTA-GCGTAAGA_L001_R1_001.fastq BRR3_ACTCGCTA-GCGTAAGA_L001_R2_001.fastq -o out_BRR3 -j
casper BRR4_GGAGCTAC-GCGTAAGA_L001_R1_001.fastq BRR4_GGAGCTAC-GCGTAAGA_L001_R2_001.fastq -o out_BRR4 -j
casper DRR1_GCGTAGTA-GCGTAAGA_L001_R1_001.fastq DRR1_GCGTAGTA-GCGTAAGA_L001_R2_001.fastq -o out_DRR1 -j
casper DRR2_CGGAGCCT-GCGTAAGA_L001_R1_001.fastq DRR2_CGGAGCCT-GCGTAAGA_L001_R2_001.fastq -o out_DRR2 -j
casper DRR3_TACGCTGC-GCGTAAGA_L001_R1_001.fastq DRR3_TACGCTGC-GCGTAAGA_L001_R2_001.fastq -o out_DRR3 -j
casper DRR4_ATGCGCAG-GCGTAAGA_L001_R1_001.fastq DRR4_ATGCGCAG-GCGTAAGA_L001_R2_001.fastq -o out_DRR4 -j
casper s52_TAGCGCTC-TTCTAGCT_L001_R1_001.fastq s52_TAGCGCTC-TTCTAGCT_L001_R2_001.fastq -o out_s52 -j
casper s53_ACTGAGCG-TTCTAGCT_L001_R1_001.fastq s53_ACTGAGCG-TTCTAGCT_L001_R2_001.fastq -o out_s53 -j
casper s54r_CCTAAGAC-TTCTAGCT_L001_R1_001.fastq s54r_CCTAAGAC-TTCTAGCT_L001_R2_001.fastq -o out_s54r -j
casper s55r_CGATCAGT-TTCTAGCT_L001_R1_001.fastq s55r_CGATCAGT-TTCTAGCT_L001_R2_001.fastq -o out_s55r -j
casper s56_TGCAGCTA-TTCTAGCT_L001_R1_001.fastq s56_TGCAGCTA-TTCTAGCT_L001_R2_001.fastq -o out_s56 -j
casper s57_TCGACGTC-TTCTAGCT_L001_R1_001.fastq s57_TCGACGTC-TTCTAGCT_L001_R2_001.fastq -o out_s57 -j
casper s58_ACTCGCTA-CCTAGAGT_L001_R1_001.fastq s58_ACTCGCTA-CCTAGAGT_L001_R2_001.fastq -o out_s58 -j
casper s59_GGAGCTAC-CCTAGAGT_L001_R1_001.fastq s59_GGAGCTAC-CCTAGAGT_L001_R2_001.fastq -o out_s59 -j
casper s60r_GCGTAGTA-CCTAGAGT_L001_R1_001.fastq s60r_GCGTAGTA-CCTAGAGT_L001_R2_001.fastq -o out_s60r -j
casper s61_CGGAGCCT-CCTAGAGT_L001_R1_001.fastq s61_CGGAGCCT-CCTAGAGT_L001_R2_001.fastq -o out_s61 -j
casper s62_TACGCTGC-CCTAGAGT_L001_R1_001.fastq s62_TACGCTGC-CCTAGAGT_L001_R2_001.fastq -o out_s62 -j
casper s63_ATGCGCAG-CCTAGAGT_L001_R1_001.fastq s63_ATGCGCAG-CCTAGAGT_L001_R2_001.fastq -o out_s63 -j
casper s65_ACTGAGCG-CCTAGAGT_L001_R1_001.fastq s65_ACTGAGCG-CCTAGAGT_L001_R2_001.fastq -o out_s65 -j
casper s66r_CCTAAGAC-CCTAGAGT_L001_R1_001.fastq s66r_CCTAAGAC-CCTAGAGT_L001_R2_001.fastq -o out_s66r -j
casper s67r_CGATCAGT-CCTAGAGT_L001_R1_001.fastq s67r_CGATCAGT-CCTAGAGT_L001_R2_001.fastq -o out_s67r -j

```

- convert fastq to fasta

```
module load seqtk/1.3
seqtk seq -a out_1A.fastq > out_1A.fasta
seqtk seq -a out_1B.fastq > out_1B.fasta
seqtk seq -a out_1C.fastq > out_1C.fasta
seqtk seq -a out_1D.fastq > out_1D.fasta
seqtk seq -a out_1E.fastq > out_1E.fasta
seqtk seq -a out_1F.fastq > out_1F.fasta
seqtk seq -a out_2A.fastq > out_2A.fasta
seqtk seq -a out_2B.fastq > out_2B.fasta
seqtk seq -a out_2C.fastq > out_2C.fasta
seqtk seq -a out_2D.fastq > out_2D.fasta
seqtk seq -a out_2E.fastq > out_2E.fasta
seqtk seq -a out_2F.fastq > out_2F.fasta
seqtk seq -a out_3A.fastq > out_3A.fasta
seqtk seq -a out_3B.fastq > out_3B.fasta
seqtk seq -a out_3C.fastq > out_3C.fasta
seqtk seq -a out_3D.fastq > out_3D.fasta
seqtk seq -a out_3E.fastq > out_3E.fasta
seqtk seq -a out_3F.fastq > out_3F.fasta
seqtk seq -a out_BRR1.fastq > out_BRR1.fasta
seqtk seq -a out_BRR2.fastq > out_BRR2.fasta
seqtk seq -a out_BRR3.fastq > out_BRR3.fasta
seqtk seq -a out_BRR4.fastq > out_BRR4.fasta
seqtk seq -a out_DRR1.fastq > out_DRR1.fasta
seqtk seq -a out_DRR2.fastq > out_DRR2.fasta
seqtk seq -a out_DRR3.fastq > out_DRR3.fasta
seqtk seq -a out_DRR4.fastq > out_DRR4.fasta
seqtk seq -a out_s52.fastq > out_s52.fasta
seqtk seq -a out_s53.fastq > out_s53.fasta
seqtk seq -a out_s54r.fastq > out_s54r.fasta
seqtk seq -a out_s55r.fastq > out_s55r.fasta
seqtk seq -a out_s56.fastq > out_s56.fasta
seqtk seq -a out_s57.fastq > out_s57.fasta
seqtk seq -a out_s58.fastq > out_s58.fasta
seqtk seq -a out_s59.fastq > out_s59.fasta
seqtk seq -a out_s60r.fastq > out_s60r.fasta
seqtk seq -a out_s61.fastq > out_s61.fasta
seqtk seq -a out_s62.fastq > out_s62.fasta
seqtk seq -a out_s63.fastq > out_s63.fasta
seqtk seq -a out_s65.fastq > out_s65.fasta
seqtk seq -a out_s66r.fastq > out_s66r.fasta
seqtk seq -a out_s67r.fastq > out_s67r.fasta
 
```

## Fungene pipeline

- Before anything can be run, need to set up the 'Options.txt' file.
	-   This file specifies the output directories
- framebot_command.txt
	- file specifies the commands to be run 

```
cat framebot_cluster_commands.txt 

dereplicate	unaligned	false
chimera_check
refresh_mapping chimera_filtered
framebot	glocal	100 0.3	false
refresh_mapping framebot_filtered
align
distance	0.3	0.20 #=GC_RF
cluster	complete	0.01
jaccard_sorensen    0.1    0.0
shannon_chao
rarefaction
explode_mapping filtered_sequences

cat options-acds-DanM.txt 

acdS
/home/julianyu/Raw/output-acds-DanM/
/home/julianyu/Raw/output-acds-DanM/pipeline-job
julianyu@asu.edu
/home/julianyu/Raw/output-acds-DanM/status.txt
/home/julianyu/Raw/output-acds-DanM/pipeline-job.tgz
/home/julianyu/Raw/output-acds-DanM/mail_message.txt
```


```
# make the base output directory specified in the options.txt file
mkdir output-acds-DanM

#load all necessary modules
module load python/3.4.3
module load biopython/1.65-python-3.4.3
module load usearch/9.2.64
module load rdptools/2.0.2
module load hmmer/3.1b2
module load infernal/1.1.2
module load ncbi-blast/2.6.0
module load r/latest


python3 ~/fungene/fungene_pipeline/fgp_wrapper.py ./options-acds-DanM.txt ./framebot_cluster_commands.txt out_1A.fasta out_1B.fasta out_1C.fasta out_1D.fasta out_1E.fasta out_1F.fasta out_2A.fasta out_2B.fasta out_2C.fasta out_2D.fasta out_2E.fasta out_2F.fasta out_3A.fasta out_3B.fasta out_3C.fasta out_3D.fasta out_3E.fasta out_3F.fasta out_BRR1.fasta out_BRR2.fasta out_BRR3.fasta out_BRR4.fasta out_DRR1.fasta out_DRR2.fasta out_DRR3.fasta out_DRR4.fasta out_s52.fasta out_s53.fasta out_s54r.fasta out_s55r.fasta out_s56.fasta out_s57.fasta out_s58.fasta out_s59.fasta out_s60r.fasta out_s61.fasta out_s62.fasta out_s63.fasta out_s65.fasta out_s66r.fasta out_s67r.fasta

```

> Running all samples will fail due to memory limits
> 
>  java.lang.OutOfMemoryError: GC overhead limit exceeded 
> 

#### Bash scripts
```
sbatch acds1.sh
Submitted batch job 2988978

#!/bin/bash
 
#SBATCH -n 20                        # number of cores
#SBATCH -N 1                                                            #number of nodes (computers)
#SBATCH -t 2-00:00                  # wall time (D-HH:MM)
##SBATCH -A julianyu             # Account hours will be pulled from (commented out with double # in front)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when the job starts, stops, or fails
#SBATCH --mail-user=julianyu@asu.edu # send-to address

module load python/3.4.3
module load biopython/1.65-python-3.4.3
module load usearch/9.2.64
module load rdptools/2.0.2
module load hmmer/3.1b2
module load infernal/1.1.2
module load ncbi-blast/2.6.0module load hmmer/3.1b2
module load r/latest

cd /home/julianyu/Raw/fastq/


python3 ~/fungene/fungene_pipeline/fgp_wrapper.py ./options-acds-DanM-1.txt ./framebot_cluster_commands.txt out_1A.fasta out_1B.fasta out_1C.fasta out_1D.fasta out_1E.fasta out_1F.fasta out_2A.fasta out_2B.fasta out_2C.fasta out_2D.fasta out_2E.fasta out_2F.fasta


sbatch acds-2.sh
Submitted batch job 2989000

#!/bin/bash
 
#SBATCH -n 20                        # number of cores
#SBATCH -N 1                                                            #number of nodes (computers)
#SBATCH -t 2-00:00                  # wall time (D-HH:MM)
##SBATCH -A julianyu             # Account hours will be pulled from (commented out with double # in front)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when the job starts, stops, or fails
#SBATCH --mail-user=julianyu@asu.edu # send-to address

module load python/3.4.3
module load biopython/1.65-python-3.4.3
module load usearch/9.2.64
module load rdptools/2.0.2
module load hmmer/3.1b2
module load infernal/1.1.2
module load ncbi-blast/2.6.0
module load r/latest

cd /home/julianyu/Raw/fastq/


python3 ~/fungene/fungene_pipeline/fgp_wrapper.py ./options-acds-DanM-2.txt ./framebot_cluster_commands.txt out_3A.fasta out_3B.fasta out_3C.fasta out_3D.fasta out_3E.fasta out_3F.fasta out_BRR1.fasta out_BRR2.fasta out_BRR3.fasta out_BRR4.fasta 


sbatch acds-3.sh
Submitted batch job 2989087

#!/bin/bash
 
#SBATCH -n 20                        # number of cores
#SBATCH -N 1                                                            #number of nodes (computers)
#SBATCH -t 2-00:00                  # wall time (D-HH:MM)
##SBATCH -A julianyu             # Account hours will be pulled from (commented out with double # in front)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when the job starts, stops, or fails
#SBATCH --mail-user=julianyu@asu.edu # send-to address

module load python/3.4.3
module load biopython/1.65-python-3.4.3
module load usearch/9.2.64
module load rdptools/2.0.2
module load hmmer/3.1b2
module load infernal/1.1.2
module load ncbi-blast/2.6.0
module load r/latest

cd /home/julianyu/Raw/fastq/


python3 ~/fungene/fungene_pipeline/fgp_wrapper.py ./options-acds-DanM-3.txt ./framebot_cluster_commands.txt out_DRR1.fasta out_DRR2.fasta out_DRR3.fasta out_DRR4.fasta out_s52.fasta out_s53.fasta out_s54r.fasta out_s55r.fasta out_s56.fasta out_s57.fasta out_s58.fasta out_s59.fasta 

sbatch acds-4.sh
Submitted batch job 2989022


#!/bin/bash
 
#SBATCH -n 20                        # number of cores
#SBATCH -N 1                                                            #number of nodes (computers)
#SBATCH -t 2-00:00                  # wall time (D-HH:MM)
##SBATCH -A julianyu             # Account hours will be pulled from (commented out with double # in front)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when the job starts, stops, or fails
#SBATCH --mail-user=julianyu@asu.edu # send-to address

module load python/3.4.3
module load biopython/1.65-python-3.4.3
module load usearch/9.2.64
module load rdptools/2.0.2
module load hmmer/3.1b2
module load infernal/1.1.2
module load ncbi-blast/2.6.0
module load r/latest

cd /home/julianyu/Raw/fastq/


python3 ~/fungene/fungene_pipeline/fgp_wrapper.py ./options-acds-DanM-4.txt ./framebot_cluster_commands.txt out_s60r.fasta out_s61.fasta out_s62.fasta out_s63.fasta out_s65.fasta out_s66r.fasta out_s67r.fasta

```

- If pipeline fails during clustering go to trace.txt and run java scripts individually 
- Clustering steps can use a lot of memory. IF the clustering scripts fail: GC memory overhead exceeded, try to increase memory
	- For example: java -Xmx3g - memory is 3Gs 	 

```
java -Xmx3g -jar /packages/7x/rdptools/2.0.2/Clustering.jar dmatrix --id-mapping /scratch/julianyu/acds/gup/output-acds-6a/pipeline-job/framebot_filtered/filtered_ids.txt --in /scratch/julianyu/acds/gup/output-acds-6a/pipeline-job/alignment/non_chimeric_prot_corr_aligned.fasta --outfile /scratch/julianyu/acds/gup/output-acds-6a/pipeline-job/dist_matrix/non_chimeric_prot_corr_aligned.fasta_matrix.bin -l 25 --dist-cutoff 0.3
java -Xmx3g -jar /packages/7x/rdptools/2.0.2/Clustering.jar cluster --method complete --id-mapping /scratch/julianyu/acds/gup/output-acds-6a/pipeline-job/framebot_filtered/filtered_ids.txt --sample-mapping /scratch/julianyu/acds/gup/output-acds-6a/pipeline-job/framebot_filtered/filtered_samples.txt --dist-file /scratch/julianyu/acds/gup/output-acds-6a/pipeline-job/dist_matrix/non_chimeric_prot_corr_aligned.fasta_matrix.bin --outfile /scratch/julianyu/acds/gup/output-acds-6a/pipeline-job/clustering/non_chimeric_prot_corr_aligned.fasta_complete.clust --step 0.01

```

### Generate OTU table 
> Transfer Fungene output to local computer
 
```
scp -r julianyu@agave.asu.edu:/home/julianyu/Raw/fastq/output-acds-DanM-3/ ~/Desktop/Gupta-MENA/Daniel-acds/
```
> Make txt file of clusters at distance cutoff 0.0 from complete cluster file

- non_chimeric_prot_corr_aligned.fasta_complete.clust
- move new 0.0-clust.txt file and framebot output to the same directory

#### Python3
> Run on each cluster and then cat all outputs to a new file

```
from collections import Counter
#Copy dist 0.0 clusters into a new file

clust_filename = input("Enter .clust filename:" )
#0.0-clust.txt
cluster =[]

with open(clust_filename, "r") as infile:
    for line in infile:
        if line[0].isdigit():
            x = line.rstrip('\n').split('\t')
            cluster.append(list(x))


framebot_filename = input('enter non_chimeric_framebot.txt file: ')
framebot = []
#non_chimeric_framebot.txt
with open(framebot_filename,'r') as infile:
    for line in infile:
        y = line.split('\t')
        framebot.append(y)

#pull out stats line from each entry
seq_taxa = []
for i in framebot:
    if "STATS" in i:
        seq_taxa.append(i[1:3])

#Match Cluster number to taxonomy
#formatting:
#Cluster# \t taxonomic lineage \t sequence id

data = []
for i in seq_taxa:
    for j in cluster:
        if i[1] in j[3]:
            data.append(j+i)


output_data = input("enter output file name: ")
with open(output_data, "w") as out_file:
    out_file.writelines('\t'.join(i) + '\n' for i in data)

## Prompt:
#Enter .clust filename:acds-3-0.0.txt
#enter non_chimeric_framebot.txt file: non_chimeric_framebot.txt
#enter output file name: acds3-0.0-out.txt
```
```
cat acds*-0.0-out.txt > acds_all-0.0-out.txt
```

> Parse file and add header to include clusterID_tax (lineage), sample, num_clust (count) 

```
cut -f 2,3,5 acds_all_out.txt > acds_all_parsed.txt
echo -e "sample\tnumclust\tAccessionID" | cat - acds_all_parsed.txt > acds_all_parsed.txt

# python to match acc id with lineage

cut -f 2,1,5 acds_all_acc_lineage.txt > acds_all_lin_parse.txt

echo -e "sample\tnumclust\tlineage" | cat - acds_all_lin_parse.txt > acds_all_lin_w_header.txt

```

* save as .csv for input into R

##### R/R-studio

```
library(reshape2)

data_raw<-read.csv('~/Desktop/Gupta-MENA/Daniel-acds/acsd_all-0.0-outR.csv', header = T)
mdata <- melt(data_raw, id=c("sample","clusterID_tax"))
dataframe<-dcast(mdata, clusterID_tax~sample) #defaulting to length

write.table(dataframe, "~/Desktop/Gupta-MENA/Daniel-acds/OTU_acds-Daniel.txt", sep = "\t")
```
