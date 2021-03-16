# Reads trimming
## Check sequence quality using FastQC
```
mkdir -p FastQCk6
mkdir -p FastQCk6/Unique
mkdir -p FastQCk6/Multi
FastQC-0.11.3/fastqc all.fastq -o FastQCk6 -k 6
```

## Adapter trimming
### Use cutadapt to cut universal sequences in adapters
```
cutadapt -a CTGTCTCTTATACAC all.fastq | \
cutadapt -a CGTATGCCGTCTTCTGCTTG - | \
cutadapt -g TGCCGTCTTCTGCTTG - | \
cutadapt -g GGTAACTTTGTGTTT - | \
cutadapt -g CTTTGTGTTTGA - | \
cutadapt -a CACTCGTCGGCAGCGTTAGATGTGTATAAG - > trimmedfastq
```

### Use Trimmomatic to trim barcode-containing adapters
```
java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads 4 -trimlog trimmomatic.log trimmedfastq alltrimmedfastq \
ILLUMINACLIP:AllAdaptors.fa:2:30:10 MAXINFO:35:0.9 MINLEN:50 
```

### Check sequence quality again using FastQC
```
FastQC-0.11.3/fastqc alltrimmedfastq -o FastQCk6 -k 6
```

## Hard-trimming to a fixed length and filter out the short reads
```
#The python2 script can be found here https://github.com/brianpenghe/fastq-scripts/blob/master/trimfastq.py
python trimfastq.py alltrimmedfastq 50 -stdout > allfastq50 
```

# Reads mapping


