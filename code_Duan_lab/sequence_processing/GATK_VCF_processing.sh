#!/bin/bash

# 26 Sept 2017 Siwei
# modified 26 Oct 2017 Siwei
# modified 1 Nov 2017 Siwei
# modified 14 Apr 2018 Siwei
vcf_suffix='_999.vcf'

date > log.txt

for EACHFILE in pre_bams/*.bam
do
echo $EACHFILE
echo $EACHFILE >> log.txt
date >> log.txt
#############call variants##############
	
	samtools index -@ 23 $EACHFILE
	
	~/2TB/Tools/oracle_JDK/bin/java -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R ~/2TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		--dbsnp ~/2TB/Databases/hg38/dbsnp_146.hg38.vcf \
		-gt_mode DISCOVERY \
		-nct 23 \
		-stand_call_conf 30 \
		-I $EACHFILE \
		-o raw_variants.vcf

############variants filtering##########
############SNP#########################
############Remove previous calibration files######

############SNP Recalibration###########


	~/2TB/Tools/oracle_JDK/bin/java -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-T VariantRecalibrator \
		-R ~/2TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		--input raw_variants.vcf \
		-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ~/2TB/Databases/hg38/hapmap_3.3.hg38.vcf \
		-resource:omni,known=false,training=true,truth=true,prior=12.0 ~/2TB/Databases/hg38/1000G_omni2.5.hg38.vcf \
		-resource:1000G,known=false,training=true,truth=false,prior=10.0 ~/2TB/Databases/hg38/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ~/2TB/Databases/hg38/dbsnp_146.hg38.vcf \
		-an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -mG 4 -recalFile raw_recalibrate_SNP.recal -tranchesFile raw_recalibrate_SNP.tranches -rscriptFile raw_recalibrate_SNP_plots.R  

###########apply SNP recalib###########

	~/2TB/Tools/oracle_JDK/bin/java -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
                -T ApplyRecalibration \
		-R ~/2TB/Databases/hg38/Homo_sapiens_assembly38.fasta \
		--input raw_variants.vcf \
		-mode SNP \
		--ts_filter_level 99.9 \
		-recalFile raw_recalibrate_SNP.recal \
		-tranchesFile raw_recalibrate_SNP.tranches \
		-o recalibrated_snps_raw_indels.vcf

############INDEL#####################
############Remove previous calibration file####


############INDEL Recalibration#######

	~/2TB/Tools/oracle_JDK/bin/java -Xmx50G -Xms10G -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
                -T VariantRecalibrator \
                -R ~/2TB/Databases/hg38/Homo_sapiens_assembly38.fasta \
                --input recalibrated_snps_raw_indels.vcf \
		-resource:mills,known=true,training=true,truth=true,prior=12.0 ~/2TB/Databases/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf \
		-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode INDEL -tranche 90.0 -tranche 99.0 -tranche 99.5 -recalFile raw_recalibrate_INDEL.recal -tranchesFile raw_recalibrate_INDEL.tranches -rscriptFile raw_recalibrate_INDEL_plots.R

##########apply INDEL recalib########

        ~/2TB/Tools/oracle_JDK/bin/java -Xmx50G -Xms10G -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
                -T ApplyRecalibration \
                -R ~/2TB/Databases/hg38/Homo_sapiens_assembly38.fasta \
                --input recalibrated_snps_raw_indels.vcf \
		-mode INDEL \
		--ts_filter_level 90.0 \
		-recalFile raw_recalibrate_INDEL.recal \
		-tranchesFile raw_recalibrate_INDEL.tranches \
		-o $EACHFILE$vcf_suffix
done

