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

## Map to transcriptome
We use STAR here for alignment, considering splicing.
```
STAR --genomeDir dm6/ --readFilesIn Unmapped50.fastq --runThreadN 8 \
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

Then sort and index the alignment files
```
samtools sort dm6.50merAligned.toTranscriptome.out.bam dm6.50merAligned.toTranscriptome.out.sorted
samtools index dm6.50merAligned.toTranscriptome.out.sorted.bam
```