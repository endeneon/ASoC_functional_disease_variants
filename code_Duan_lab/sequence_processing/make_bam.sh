#!/bin/sh

# Siwei 3 Apr 2018

# Siwei 13 Sept 2018

for EACHFILE_1P in *1P.fastq.gz
do
	####init
	date
	rm *.bam
	rm *.list
	rm *.table
	rm *.bai

	####variables
        echo $EACHFILE_1P
	echo ${EACHFILE_1P/%1P.fastq.gz/2P.fastq.gz}  #variable for 2P
        echo ${EACHFILE_1P/%1P.fastq.gz/1U.fastq.gz}  #variable for 1U
	echo ${EACHFILE_1P/%1P.fastq.gz/2U.fastq.gz}  #variable for 2U

	####alignments

        ####align for unpaired reads
	date >> stat.txt
	echo "Unpaired"
	echo ${EACHFILE_1P/%1P.fastq.gz/U} >> stat.txt
	(bowtie2 -p 23 -X 2000 --mm --qc-filter --met 1 -t --sensitive --no-mixed --no-discordant -x ~/2TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta -U ${EACHFILE_1P/%1P.fastq.gz/1U.fastq.gz} -U ${EACHFILE_1P/%1P.fastq.gz/2U.fastq.gz} | samtools view -b -u -h -@ 23 > Unpaired.bam) 2>>stat.txt

	####align for paired reads
	echo "Paired"
	echo ${EACHFILE_1P/%1P.fastq.gz/P} >> stat.txt
	(bowtie2 -p 23 -X 2000 --mm --qc-filter --met 1 -t --sensitive --no-mixed --no-discordant -x ~/2TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta -1 $EACHFILE_1P -2 ${EACHFILE_1P/%1P.fastq.gz/2P.fastq.gz} | samtools view -b -u -h -@ 23 > Paired.bam) 2>>stat.txt
	ls -lh
	
	####merge reads and sort
	samtools merge -l 0 -@ 23 -f merged_unsorted.bam Paired.bam Unpaired.bam
	samtools sort -l 9 -m 2G -@ 23 -o merged_sorted.bam merged_unsorted.bam
	samtools index -@ 23 merged_sorted.bam
	ls -lh

	####dedup
	~/2TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/2TB/Tools/picard291.jar MarkDuplicates INPUT=merged_sorted.bam OUTPUT=merged_sorted_dedup.bam METRICS_FILE=metrics.txt REMOVE_DUPLICATES=True ASSUME_SORTED=True
	samtools index merged_sorted_dedup.bam

	####add read group
	~/2TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/2TB/Tools/picard291.jar AddOrReplaceReadGroups \
		INPUT=merged_sorted_dedup.bam \
		OUTPUT=merged_sorted_dedup_RG.bam \
		RGID=${EACHFILE_1P/%_R1_001_1P.fastq.gz/} \
		RGLB=lib1 \
		RGPL=illumina \
		RGPU=unit1 \
		RGSM=${EACHFILE_1P/%_R1_001_1P.fastq.gz/}

	samtools index merged_sorted_dedup_RG.bam

	####Re-align reads around known indels
	~/2TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-R ~/2TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		-I merged_sorted_dedup_RG.bam \
		-known ~/2TB/Databases/hg38/INDEX/Mills_and_1000G_gold_standard.indels.hg38.vcf \
		-o realignment_targets.list

	~/2TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R ~/2TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		-I merged_sorted_dedup_RG.bam \
		-known ~/2TB/Databases/hg38/INDEX/Mills_and_1000G_gold_standard.indels.hg38.vcf \
		-targetIntervals realignment_targets.list \
		-o merged_sorted_dedup_RG_realigned.bam

	samtools index -@ 23 merged_sorted_dedup_RG_realigned.bam

	####Base score recalibration
	~/2TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-T BaseRecalibrator \
		-R ~/2TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		-I merged_sorted_dedup_RG_realigned.bam \
		-knownSites ~/2TB/Databases/hg38/INDEX/Mills_and_1000G_gold_standard.indels.hg38.vcf \
		-knownSites ~/2TB/Databases/hg38/dbsnp_146.hg38.vcf \
		-o recal_data.table 

	~/2TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-T PrintReads \
		-R ~/2TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		-I merged_sorted_dedup_RG_realigned.bam \
		-BQSR recal_data.table \
		-o pre_bams_2/${EACHFILE_1P/%_R1_001_1P.fastq.gz/}_new.bam 
	
		samtools index pre_bams_2/${EACHFILE_1P/%_R1_001_1P.fastq.gz/}_new.bam

done

