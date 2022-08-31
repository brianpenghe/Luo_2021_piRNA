# Reads trimming

## Check sequence quality using FastQC
```
mkdir -p FastQCk6
FastQC-0.11.3/fastqc all.fastq -o FastQCk6 -k 6
```

## Reads trimming
If you know the sequences of some prominent adapters from the previous step, you can use cutadapt to trim your sequences first.

Then, trim the reads to 50bp:

`python2 trimfastq.py all.fastq 50 stdout > allfastq50`

[trimfastq.py](https://github.com/brianpenghe/fastq-scripts/blob/master/trimfastq.py) is written by Georgi Marinov.

### Check sequence quality again using FastQC

`FastQC-0.11.3/fastqc allfastq50 -o FastQCk6 -k 6`

# Map reads
## Ribosomal reads removal
Next, reads are mapped against ribosomal RNA sequence references. The unmapped reads are saved.
```
bowtie dmel_rRNA_unit -p 8 -v 2 -k 1 --best -t --sam-nh -q \
    allfastq50 --un Unmapped50.fastq allfastq.rRNA.mapped50.map
```
For more about output files, please read the manual of [bowtie](http://bowtie-bio.sourceforge.net/manual.shtml).

## Map to transcriptome
We use STAR here for alignment, considering splicing.
```
STAR --genomeDir dm6/ --readFilesIn allfastqrRNAUnmapped50.fastq --runThreadN 8 \
     --genomeLoad NoSharedMemory --outFilterMultimapNmax 1 --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.0 --alignIntronMin 20 \
     --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
     --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate \
     --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
     --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \
     --sjdbScore 1 --limitBAMsortRAM 30000000000 --outFileNamePrefix dm6.50mer \
     --outSAMstrandField intronMotif
```
For more about output files, please read the manual of [STAR](https://github.com/alexdobin/STAR).

Then sort and index the alignment files
```
samtools sort dm6.50merAligned.toTranscriptome.out.bam dm6.50merAligned.toTranscriptome.out.sorted
samtools index dm6.50merAligned.toTranscriptome.out.sorted.bam
```
## Map reads on vectors
```
bowtie UBIG -p 8 --chunkmbs 1024 -v 0 -a -m 1 -t --sam-nh --best --strata -q --sam \
    allfastqrRNAUnmapped50.fastq -k 1 --al dm6.50mer.UBIG.vectoronly.fastq | \
    samtools view -F 4 -bT UBIG.fa - | samtools sort - dm6.50mer.UBIG.vectoronly.dup
```

Then index the alignment files

`samtools index dm6.50mer.UBIG.vectoronly.dup.bam`


# Further process
## Generate bedgraph coverage files
```
STAR --runMode inputAlignmentsFromBAM --inputBAMfile dm6.50merAligned.sortedByCoord.out.bam \
     --outWigType bedGraph --outWigStrand Unstranded --outFileNamePrefix \
     dm6.50merAligned.sortedByCoord --outWigReferencesPrefix chr
```
## Generate read coverage profiles
```
python /woldlab/castor/home/sau/code/gene_coverage_wig_gtf.py Drosophila_melanogaster.BDGP6.32.106.fixed.gtf \
     dm6.50merAligned.sortedByCoordSignal.Unique.str1.out.bgc 1000 dm6.50merAligned.sortedByCoord.coverage
```
[gene_coverage_wig_gtf](https://github.com/brianpenghe/gtfgff-scripts/blob/master/gene_coverage_wig_gtf.py) was written by Georgi Marinov and Sean Upchurch.

## Make bigWig files

```
wigToBigWig -clip dm6.50merAligned.sortedByCoordSignal.Unique.str1.out.bg dm6.chrom.sizes \
    dm6.50merAligned.sortedByCoordSignal.Unique.str1.out.bg.bigWig
wigToBigWig -clip dm6.50merAligned.sortedByCoordSignal.UniqueMultiple.str1.out.bg dm6.chrom.sizes \
    dm6.50merAligned.sortedByCoordSignal.UniqueMultiple.str1.out.bg.bigWig
```

## Transcript quantification
We use Rsem to quantify the transcript abundance.
```
rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed 12345 -p 8 \
    --no-bam-output --ci-memory 30000  --temporary-folder temp \
    dm6.50merAligned.toTranscriptome.out.sorted.bam dm6/rsem \
    dm6.50merAligned.toTranscriptome.out.sorted.rsem
```
For more about output files, please read the manual of [Rsem](https://github.com/deweylab/RSEM).

## Calculate read counts per equal-sized bins in vectors
```
Coordinates_20A='chrX:21519548:21560880'
Coordinates_42AB='chr2R:6256844:6499214'
Coordinates_Flam='chrX:21632639:21883809'

ls *vectoronly*.bam | rev | cut -d. -f2- | rev > bams
for i in 10 100 1000
  do
    while read bam
      do
        bamCoverage -b $bam.bam \
-of bedgraph -bs $i -o $bam.$i.bg4
        bamCoverage -b $bam.bam \
-of bedgraph -bs $i --samFlagInclude 16 -o $bam.$i.Minus.bg4
        bamCoverage -b $bam.bam \
-of bedgraph -bs $i --samFlagExclude 16 -o $bam.$i.Plus.bg4
      done<bams
  done

#bin counts for genome region mappings
ls */*merAligned.sortedByCoord.out.bam > bamsgenome
for i in 10 100 1000
  do
    while read bam
      do
        /woldlab/castor/proj/programs/samtools-0.1.16/bin/samtools index $bam
        bamCoverage -b $bam \
-of bedgraph -bs $i --region $Coordinates_20A -o $bam.$i.20A.bg4
        bamCoverage -b $bam \
-of bedgraph -bs $i --region $Coordinates_20A --samFlagInclude 16 -o $bam.$i.20A.Minus.bg4
        bamCoverage -b $bam \
-of bedgraph -bs $i --region $Coordinates_20A --samFlagExclude 16 -o $bam.$i.20A.Plus.bg4
        bamCoverage -b $bam \
-of bedgraph -bs $i --region $Coordinates_42AB -o $bam.$i.42A.bg4
        bamCoverage -b $bam \
-of bedgraph -bs $i --region $Coordinates_42AB --samFlagInclude 16 -o $bam.$i.42A.Minus.bg4
        bamCoverage -b $bam \
-of bedgraph -bs $i --region $Coordinates_42AB --samFlagExclude 16 -o $bam.$i.42A.Plus.bg4
        bamCoverage -b $bam \
-of bedgraph -bs $i --region $Coordinates_Flam -o $bam.$i.flamenco.bg4
        bamCoverage -b $bam \
-of bedgraph -bs $i --region $Coordinates_Flam --samFlagInclude 16 -o $bam.$i.flamenco.Minus.bg4
        bamCoverage -b $bam \
-of bedgraph -bs $i --region $Coordinates_Flam --samFlagExclude 16 -o $bam.$i.flamenco.Plus.bg4
      done<bamsgenome
  done

#chop to retrieve empty bins
for i in 10 100 1000
  do ls {*.$i.*bg4,*k6/*.$i.*bg4} > bg4.$i.list
    while read bg4
      do
        awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' $bg4 > signal.bed
        bedops --chop $i signal.bed | bedmap --echo --echo-map-score - signal.bed | \
sed -e 's/|/\t/g' > $bg4.chopped.bg4
      done<bg4.$i.list
  done
```

