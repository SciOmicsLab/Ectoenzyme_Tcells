#!/bin/bash

set -e
set -u
set -o pipefail

###RNA-Seq Analysis Pipeline###

#add code to check for programs
#fastqc
#trim_galore
#cut_adapt
#pigz
#hisat2
#samtools
#featurecounts
#multiqc

#caffeinate in background
caffeinate &

#write text file with sample names
#find $FASTQ_FOLDER -name "*.fastq*"> sample_names.txt

#Set folder variables
ROOT_FOLDER=/Users/david/Bioinformatics/Teee_RNA_Seq
FASTQ_FOLDER=/Users/david/Bioinformatics/Teee_RNA_Seq/fastq_files
REFERENCE_GENOME=/Users/david/Bioinformatics/ref_genomes/hg38_ensmbl

#Set sample name array
SAMPLE_NAMES=(Sample_1_CD38_neg_CD39_neg \
Sample_1_CD38_neg_CD39_pos \
Sample_1_CD38_pos_CD39_neg \
Sample_1_CD38_pos_CD39_pos \
Sample_2_CD38_neg_CD39_neg \
Sample_2_CD38_neg_CD39_pos \
Sample_2_CD38_pos_CD39_neg \
Sample_2_CD38_pos_CD39_pos \
Sample_3_CD38_neg_CD39_neg \
Sample_3_CD38_neg_CD39_pos \
Sample_3_CD38_pos_CD39_neg \
Sample_3_CD38_pos_CD39_pos)

#Set hyperparameter variables
THREADS=22

#IF NEEDED: index reference genome  
    #download reference genome
        #curl ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -o GRCh38.dna.primary_assembly.fa.gz
        #curl ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz  -o GRCh38.84.gtf.gz
        #mv GRCh38.dna.primary_assembly.fa.gz ${REFERENCE_GENOME}
        #mv GRCh38.84.gtf.gz ${REFERENCE_GENOME}
        #gunzip GRCh38.dna.primary_assembly.fa.gz
        #gunzip GRCh38.84.gtf.gz
        #hisat2_extract_splice_sites.py GRCh38.84.gtf > genome.ss
        #hisat2_extract_exons.py GRCh38.84.gtf > genome.exon
    #index reference genome
        #hisat2-build \
            #-p $THREADS \
            #--exon ${REFERENCE_GENOME}/genome.exon \
            #--ss ${REFERENCE_GENOME}/genome.ss \
            #${REFERENCE_GENOME}/GRCh38.dna.primary_assembly.fa \
            #hg38_ensmbl

#Create a log file
touch ${ROOT_FOLDER}/log.txt
echo "Pipeline started on $(date +%y-%m-%d\ %H:%M:%S)" "\n" >> log.txt
echo Root folder location: $ROOT_FOLDER >> log.txt
echo Fastq folder location: $FASTQ_FOLDER >> log.txt
echo Reference genome location: $REFERENCE_GENOME "\n" >> log.txt
echo Sample Names: $SAMPLE_NAMES "\n" >> log.txt

#assess quality with fastqc
mkdir $ROOT_FOLDER/fastqc_results
echo "Performing quality assessment using fastqc"
echo "Started fastqc at $(date +%y-%m-%d\ %H:%M:%S)" >> $ROOT_FOLDER/log.txt
fastqc -o ${ROOT_FOLDER}/fastqc_results -t ${THREADS} ${FASTQ_FOLDER}/*.fastq*
echo "Done performing quality assessment"
echo "Finished fastqc at $(date +%y-%m-%d\ %H:%M:%S)" >> $ROOT_FOLDER/log.txt
echo "fastqc results located in ${ROOT_FOLDER}/fastqc_results" "\n" >> $ROOT_FOLDER/log.txt

#remove adapter sequences and reads with <20 PHRED score
mkdir $ROOT_FOLDER/trimmed_files
mkdir $ROOT_FOLDER/trimmed_files/reports
echo "Removing adapter sequences using trim_galore"
echo "Started trim galore at $(date +%y-%m-%d\ %H:%M:%S)" >> $ROOT_FOLDER/log.txt
for sample in $SAMPLE_NAMES
	do
	    trim_galore \
        --cores 8 \
        --paired ${FASTQ_FOLDER}/${sample}_R1.fastq* \
        ${FASTQ_FOLDER}/${sample}_R2.fastq* \
        -o ${ROOT_FOLDER}/trimmed_files
	done
mv ${ROOT_FOLDER}/trimmed_files/*.txt ${ROOT_FOLDER}/trimmed_files/reports
echo "Done removing adapter sequences"
echo "Finished trimming at $(date +%y-%m-%d\ %H:%M:%S)" >> $ROOT_FOLDER/log.txt
echo "Trimmed files results located in ${ROOT_FOLDER}/trimmed_files" >> $ROOT_FOLDER/log.txt
echo "Reports located in ${ROOT_FOLDER}/trimmed_files/reports" "\n" >> $ROOT_FOLDER/log.txt

#map reads
##need to output alignment rate stats##
mkdir $ROOT_FOLDER/bam_files
echo "Started mapping reads at $(date)" >> ${ROOT_FOLDER}/log.txt
for sample in $SAMPLE_NAMES
	do
        echo "Starting mapping of ${sample} at $(date)"
        hisat2 \
            -p $THREADS \
            -x $REFERENCE_GENOME/hg38_ensmbl \
            -1 ${ROOT_FOLDER}/trimmed_files/${sample}_R1* \
            -2 ${ROOT_FOLDER}/trimmed_files/${sample}_R2* \
            -S ${ROOT_FOLDER}/bam_files/${sample}.sam 
        echo "Finished mapping of ${sample} at $(date)"
        #convert sam files to bam files
        echo "Started conversion of ${sample}.sam to bam $(date)"
        samtools view -b \
            -@ $THREADS \
            ${ROOT_FOLDER}/bam_files/${sample}.sam > ${ROOT_FOLDER}/bam_files/${sample}.bam 
        rm ${ROOT_FOLDER}/bam_files/${sample}.sam
        echo "Finished conversion of ${sample}.sam to bam $(date)"
        #echo total time for sample
 	done
echo "Finished mapping reads at $(date)" >> ${ROOT_FOLDER}/log.txt

#sort and index reads
echo "Started sorting and indexing reads at $(date)" >> ${ROOT_FOLDER}/log.txt
for sample in $SAMPLE_NAMES
    do
        #sort
        echo "Starting sorting ${sample} at $(date)"
        samtools sort \
            -@ $THREADS \
            ${ROOT_FOLDER}/bam_files/${sample}.bam \
            -o ${ROOT_FOLDER}/bam_files/${sample}_sorted.bam
        echo "Finished sorting ${sample} at $(date)"
        #index
        echo "Started indexing ${sample} at $(date)"
        samtools index \
            -@ $THREADS \
            ${ROOT_FOLDER}/bam_files/${sample}_sorted.bam
        echo "Finished indexing ${sample} at $(date)"
        rm ${ROOT_FOLDER}/bam_files/${sample}.bam
    done
echo "Finished sorting and indexing reads at $(date)" >> ${ROOT_FOLDER}/log.txt

#generate count table
echo "Started generating read count table at $(date)" 
echo "Started generating read count table at $(date)" >> ${ROOT_FOLDER}/log.txt
featureCounts \
    -T $THREADS \
    -p \
    -s 0 \
    -a ${REFERENCE_GENOME}/GRCh38.84.gtf \
    -o ${ROOT_FOLDER}/counts.out \
    ${ROOT_FOLDER}/bam_files/*_sorted.bam
echo "Finished generating read count table at $(date)" >> ${ROOT_FOLDER}/log.txt
echo "Finished generating read count table at $(date)" 

#MultiQC Report
multiqc ${ROOT_FOLDER}

#End logging
echo "Pipeline finished on $(date +%y-%m-%d\ %H:%M:%S)" "\n" >> log.txt