#!/bin/bash

# Check if the number of arguments is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <sample_number>"
    exit 1
fi

# Set the sample number as a variable
SAMPLE_NUMBER=$1


# Step 1: Trimming with Trimmomatic
trimmomatic SE -phred33 sample${SAMPLE_NUMBER}.fastq.gz sample_trimmed${SAMPLE_NUMBER}.fastq.gz \
ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 2: Alignment with HISAT2
hisat2 -x GRCh38_index -U sample_trimmed${SAMPLE_NUMBER}.fastq.gz --phred33 -S sample_aligned${SAMPLE_NUMBER}.sam

# Step 3: Convert SAM to BAM, sort, and index
samtools view -S -b sample_aligned${SAMPLE_NUMBER}.sam > sample_aligned${SAMPLE_NUMBER}.bam
samtools sort sample_aligned${SAMPLE_NUMBER}.bam -o sample_sorted${SAMPLE_NUMBER}.bam
samtools index sample_sorted${SAMPLE_NUMBER}.bam

# Step 4: Remove the .sam file to save space
rm sample_aligned${SAMPLE_NUMBER}.sam

# Step 5: Run StringTie to assemble transcripts
stringtie sample_sorted${SAMPLE_NUMBER}.bam -G GRCh38_genomic.gtf -e -B \
-o sample${SAMPLE_NUMBER}_assembled.gtf -A sample${SAMPLE_NUMBER}_abundance.txt
