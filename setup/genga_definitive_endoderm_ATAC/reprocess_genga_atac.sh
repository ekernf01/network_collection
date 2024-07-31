#!/bin/bash

# Necessary evil; see https://github.com/ncbi/sra-tools/issues/77    
vdb-config --interactive 
cd from_to_docker
mv ../TruSeq3-PE-2.fa  . # Sorry, I could not get the image to build with this in the right directory.

# hg38 bowtie index
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip
conda install -c bioconda trimmomatic=0.39=hdfd78af_2

for id in  SRR8633116     SRR8633117 SRR8633130 SRR8633131; do
    prefetch $id
    fastq-dump --gzip "${id}/${id}.sra"
    reads_file="${id}.fastq.gz"
    trimmed="${id}_trimmed.fastq.gz"
    aligned="${id}.bam"
    dedup="${id}_dedup.bam"
    macs2_output="${id}.macs2_output"
    trimmomatic SE -threads 15 -summary "${id}_trimsummary.txt" $reads_file $trimmed ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 && rm $reads_file
    bowtie2 -x GRCh38_noalt_as/GRCh38_noalt_as -U $trimmed --threads 15 | samtools view -bh -o $aligned && rm $trimmed
    samtools rmdup $aligned $dedup && rm $aligned
    # There are differences between ATAC and ChIP. 
    # No mock IP or KO controls, so --nolambda is recommended
    # No weird shifting: insertion happens directly in open regions. So --nomodel is recommended. 
    # Shift 100 extsize 200 will smooth out the insertion density like a KDE with a radius-100 uniform kernel. It keeps the original insertion site at the center.
    # Fragments span inaccessible regions between two insertion events. We want the accessible regions. So, for peak calling, **ignore** the paired-endedness of the reads.
    macs2 callpeak -t $dedup -n $macs2_output -g hs --nolambda --nomodel --shift -100 --extsize 200 
done

bedtools intersect -a SRR8633116.macs2_output_peaks.narrowPeak -b SRR8633117.macs2_output_peaks.narrowPeak > peaks_intersected_psc.bed
bedtools intersect -a SRR8633130.macs2_output_peaks.narrowPeak -b SRR8633131.macs2_output_peaks.narrowPeak > peaks_intersected_end.bed
# motif analysis
python ../generate_network.py peaks_intersected_psc.bed psc.parquet
python ../generate_network.py peaks_intersected_end.bed endoderm.parquet
