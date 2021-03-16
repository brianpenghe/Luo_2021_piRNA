# Reads trimming
## Check sequence quality using FastQC
```
mkdir -p FastQCk6
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

## Hard-trimming to 50bp and filter out the shorter reads
```
#The python2 script can be found here https://github.com/brianpenghe/fastq-scripts/blob/master/trimfastq.py
python trimfastq.py alltrimmedfastq 50 -stdout > allfastq50 
```

# Reads mapping

## Map to the genome
### Use Bowtie to map reads to the genome index
```
bowtie-1.0.1/bowtie genome/bowtie-indexes/dm3 -p 8 -v 2 -k 1 -m 1 -t --sam-nh --best -y --strata -q \
 --sam allfastq50 | samtools-0.1.8/samtools view -bT genome/bowtie-indexes/dm3.fa - | \
 samtools-0.1.16/bin/samtools sort - dm3.50mer.unique.dup
samtools-0.1.16/bin/samtools index dm3.50mer.unique.dup.bam
```
### Remove mitochondria reads
```
samtools-0.1.16/bin/samtools view dm3.50mer.unique.dup.bam | egrep -v chrM | \
samtools-0.1.8/samtools view -bT genome/bowtie-indexes/dm3.fa - -o dm3.50mer.unique.dup.nochrM.bam
samtools-0.1.16/bin/samtools index dm3.50mer.unique.dup.nochrM.bam
```
##### (Optional) Remove duplicate reads in case useful
```
samtools-0.1.16/bin/samtools rmdup -s dm3.50mer.unique.dup.nochrM.bam dm3.50mer.unique.nochrM.bam
samtools-0.1.16/bin/samtools index dm3.50mer.unique.nochrM.bam
```
## Map to vectors
### Use Bowtie to map reads to the vector index
```
bowtie-1.0.1/bowtie genomes/YichengVectors/42AB_UBIG -p 8 -v 2 -k 1 -m 1 -t --sam-nh --best -y --strata -q \
 --sam allfastq50 | samtools-0.1.8/samtools view -bT genomes/YichengVectors/42AB_UBIG.fa - | \
 samtools-0.1.16/bin/samtools sort - dm3.50mer.42AB_UASG.vectoronly.dup
samtools-0.1.16/bin/samtools index dm3.50mer.42AB_UASG.vectoronly.dup.bam
```

# Further processing

## Make piled-up bigWig files for browser track visualization
```
#The python2 script can be found here https://github.com/brianpenghe/bamsam-scripts/blob/master/makewigglefromBAM-NH.py 
python makewigglefromBAM-NH.py --- 50mer.unique.dup.nochrM.bam genome/bowtie-indexes/dm3.chrom.sizes \
dm3.50mer.unique.bg4 -notitle -uniqueBAM -RPM
wigToBigWig -clip dm3.50mer.unique.bg4 genome/bowtie-indexes/dm3.chrom.sizes dm3.50mer.unique.bigWig
```

## Make ChIP-vs-Input enrichment track
```
export PYTHONPATH=$PYTHONPATH:~/deepTools-2.4.2_develop/lib/python2.7/site-packages
deepTools-2.4.2_develop/bin/bamCompare -b1 ChIP.dm3.50mer.unique.dup.nochrM.bam -b2 Input.dm3.50mer.unique.dup.nochrM.bam \
    -of "bigwig" -o ChIP1.Enrichment.bigWig --binSize 10 -bl dm3-blacklist.bed -p 8 
#The blacklist files can be found here: https://github.com/Boyle-Lab/Blacklist/
```

## Calculate read counts per equal-sized bin in vectors
```
ls *vectoronly*.bam | rev | cut -d. -f2- | rev > bams

for i in 10 100 1000
    do 
        while read bam
            do 
              deepTools-2.4.2_develop/bin/bamCoverage \
             -b $bam.bam  -of bedgraph -bs $i -o $bam.$i.bg4 
              deepTools-2.4.2_develop/bin/bamCoverage \
             -b $bam.bam  -of bedgraph -bs $i --samFlagInclude 16 -o $bam.$i.Minus.bg4
              deepTools-2.4.2_develop/bin/bamCoverage \
             -b $bam.bam  -of bedgraph -bs $i --samFlagExclude 16 -o $bam.$i.Plus.bg4
            done<bams
    done

for i in 10 100 1000
    do
        ls *.$i.*bg4 > bg4.$i.list
        while read bg4
            do
                awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' $bg4 > signal.bed
                bedops --chop $i signal.bed | bedmap --echo --echo-map-score - signal.bed \
                       | sed -e 's/|/\t/g' > $bg4.chopped.bg4
            done<bg4.$i.list
    done
```

## Calculate read counts per equal-sized bin in transposons
```
for i in 10 50 100 500 1000
  do
    while read bam
      do
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chr2R:2144349-2386719 -o $bam.$i.42AB.bg4
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chr2R:2144349-2386719 --samFlagInclude 16 -o $bam.$i.42AB.Minus.bg4 
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chr2R:2144349-2386719 --samFlagExclude 16 -o $bam.$i.42AB.Plus.bg4
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chrX:21392175-21431907 -o $bam.$i.20A.bg4
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chrX:21392175-21431907  --samFlagInclude 16 -o $bam.$i.20A.Minus.bg4
        deepTools-2.4.2_develop/bin/bamCoverage -b $bam  -of bedgraph -bs $i --region chrX:21392175-21431907  --samFlagExclude 16 -o $bam.$i.20A.Plus.bg4

        awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' $bam.$i.42AB.Minus.bg4 > signal.bed; bedops --chop $i signal.bed | bedmap --echo --echo-map-score - signal.bed | sed -e 's/|/\t/g' > $bam.$i.42AB.Minus.bg4chopped.bg4
        awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' $bam.$i.42AB.Plus.bg4 > signal.bed; bedops --chop $i signal.bed | bedmap --echo --echo-map-score - signal.bed | sed -e 's/|/\t/g' > $bam.$i.42AB.Plus.bg4chopped.bg4
        awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' $bam.$i.42AB.bg4 > signal.bed; bedops --chop $i signal.bed | bedmap --echo --echo-map-score - signal.bed | sed -e 's/|/\t/g' > $bam.$i.42AB.bg4chopped.bg4
        awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' $bam.$i.20A.Minus.bg4 > signal.bed; bedops --chop $i signal.bed | bedmap --echo --echo-map-score - signal.bed | sed -e 's/|/\t/g' > $bam.$i.20A.Minus.bg4chopped.bg4
        awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' $bam.$i.20A.Plus.bg4 > signal.bed; bedops --chop $i signal.bed | bedmap --echo --echo-map-score - signal.bed | sed -e 's/|/\t/g' > $bam.$i.20A.Plus.bg4chopped.bg4
        awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' $bam.$i.20A.bg4 > signal.bed; bedops --chop $i signal.bed | bedmap --echo --echo-map-score - signal.bed | sed -e 's/|/\t/g' > $bam.$i.20A.bg4chopped.bg4

      done<<<$(ls *dm3.50mer.unique.dup.bam)
  done
```

