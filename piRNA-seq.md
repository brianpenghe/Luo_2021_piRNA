# Reads trimming and filtering
## Check sequence quality using FastQC
```
mkdir -p FastQCk6
FastQC-0.11.3/fastqc all.fastq -o FastQCk6 -k 6
```
## Adapter trimming
```
cutadapt -a AGATCGGAAGAGCACACGTCT -m 15 all.fastq > alltrimmedfastq
```
### Check sequence quality again using FastQC
```
FastQC-0.11.3/fastqc alltrimmedfastq -o FastQCk6 -k 6
```
## Extract reads between 23 and 29 bp and FastQC
```
cat alltrimmedfastq | paste - - - - | awk -F"\t" 'length($2) >= 23 && length($2) <= 29' | sed 's/\t/\n/g' > allfastq23_29 
FastQC-0.11.3/fastqc allfastq23_29 -o FastQCk6 -k 6
```
## Extract reads of 21 bp and FastQC
```
cat alltrimmedfastq | paste - - - - | awk -F"\t" 'length($2) >= 21 && length($2) <= 21' | sed 's/\t/\n/g' > allfastq21_21 
FastQC-0.11.3/fastqc allfastq21_21 -o FastQCk6 -k 6 
```

# Map reads
## Map to genome
```
bowtie genome/bowtie-indexes/dm3 -p 8 --chunkmbs 1024 -v 0 -a -m 1 -t --sam-nh --best --strata -q --sam \
allfastq23_29 -k 1 --al dm3.23_29mer.unique.fastq --un dm3.23_29_unmapped.fastq --max dm3.23_29_multi.fastq | \
samtools-0.1.8/samtools view -bT  genome/bowtie-indexes/dm3.fa - | \
samtools-0.1.16/bin/samtools sort - dm3.23_29mer.unique.dup
samtools-0.1.16/bin/samtools index dm3.23_29mer.unique.dup.bam
```
## Remove mitochondria reads
```
samtools-0.1.16/bin/samtools view dm3.23_29mer.unique.dup.bam | egrep -v chrM | \
samtools-0.1.8/samtools view -bT genome/bowtie-indexes/dm3.fa - -o dm3.23_29mer.unique.dup.nochrM.bam
samtools-0.1.16/bin/samtools index dm3.23_29mer.unique.dup.nochrM.bam
```
## Map to vector
```
bowtie genomes/YichengVectors/42AB_UASG -p 8 --chunkmbs 1024 -v 0 -a -m 1 -t --sam-nh --best --strata -q --sam \
allfastq23_29 -k 1 --al dm3.23_29mer.42AB_UASG.vectoronly.fastq | \
samtools-0.1.8/samtools view -bT genomes/YichengVectors/42AB_UASG.fa - | \
samtools-0.1.16/bin/samtools sort - dm3.23_29mer.42AB_UASG.vectoronly.dup
samtools-0.1.16/bin/samtools index dm3.23_29mer.42AB_UASG.vectoronly.dup.bam
```
# Further process
## Calculate read counts per equal-sized bins in vectors
```
export PYTHONPATH=$PYTHONPATH:~/deepTools-2.4.2_develop/lib/python2.7/site-packages
ls *vectoronly*.bam | rev | cut -d. -f2- | rev > bams
for i in 10 100 1000
  do 
    while read bam
      do 
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam.bam -of bedgraph -bs $i -o $bam.$i.bg4
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam.bam -of bedgraph -bs $i --samFlagInclude 16 -o $bam.$i.Minus.bg4
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam.bam -of bedgraph -bs $i --samFlagExclude 16 -o $bam.$i.Plus.bg4
      done<bams
  done

for i in 10 100 1000
  do ls *.$i.*bg4 > bg4.$i.list
    while read bg4
      do 
        awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' $bg4 > signal.bed
        bedops --chop $i signal.bed | bedmap --echo --echo-map-score - signal.bed | sed -e 's/|/\t/g' > $bg4.chopped.bg4
      done<bg4.$i.list

  done
```
## Make browser-track bigWig files
```
while read bam
  do
    deepTools-2.4.2_develop/bin/bamCoverage -b $bam -of bigwig -bs 1 --samFlagInclude 16 -o $bam.1.Minus.bigWig;
    deepTools-2.4.2_develop/bin/bamCoverage -b $bam -of bigwig -bs 1 --samFlagExclude 16 -o $bam.1.Plus.bigWig;
  done<<<$(ls *.dm3.23_29mer.unique.dup.bam)
```
## Inspect ping-pong signature
```
#The signature.py program can be downloaded from http://drosophile.org/GEDlab/?page_id=730
#reference https://www.researchgate.net/publication/263054411_Computing_siRNA_and_piRNA_overlap_signatures
while read bam
  do 
    samtools view -h $bam > $bam.sam &&
    python2 signature.py $bam.sam 23 29 1 29 $bam.pingpong & 
  done <<<$(ls *23_29mer.{GFP_SV40,UBIG}.vectoronly.dup.bam)
```
## Calculate read counts per equal-sized bins in transposable elements
```
ls *23_29mer.unique*.bam > bamsgenome

for i in 10 100 1000
  do 
    while read bam
      do
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chrX:21392175-21431907 -o $bam.$i.20A.bg4
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chrX:21392175-21431907 --samFlagInclude 16 -o $bam.$i.20A.Minus.bg4
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chrX:21392175-21431907 --samFlagExclude 16 -o $bam.$i.20A.Plus.bg4
      done<bamsgenome
  done



for i in 10 100 1000
  do 
    while read bam
      do
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chr2R:2144349-2386719 -o $bam.$i.42A.bg4
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chr2R:2144349-2386719 --samFlagInclude 16 -o $bam.$i.42A.Minus.bg4
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chr2R:2144349-2386719 --samFlagExclude 16 -o $bam.$i.42A.Plus.bg4
      done<bamsgenome
  done

for i in 10 100 1000
  do 
    ls *.$i.*bg4 > bg4.$i.list
    while read bg4
      do 
        awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' $bg4 > signal.bed
        bedops --chop $i signal.bed | bedmap --echo --echo-map-score - signal.bed | sed -e 's/|/\t/g' > $bg4.chopped.bg4
      done<bg4.$i.list
  done
```
## Inspect ping-pong signature in transposable elements
```
while read bam
  do 
    samtools view -h $bam  chr2R:2144349-2386719 > $bam.sam && 
    python2 signature.py $bam.sam 23 29 1 29 $bam.cluster_42AB.pingpong & 
  done<<<$(ls *.dm3.23_29mer.unique.dup.bam)

while read bam
  do 
    samtools view -h $bam  chrX:21392175-21431907 > $bam.sam &&
    python2 signature.py $bam.sam 23 29 1 29 $bam.cluster_20A.pingpong & 
  done<<<$(ls *.dm3.23_29mer.unique.dup.bam)
```   

## Phasing analysis
A [jupyter notebook](https://github.com/brianpenghe/Luo_2021_piRNA/blob/main/Phasing/220522YichengPhasing.ipynb) is included in this repository to show how phasing analysis can be done.
