#!/bin/bash

# 26 Sept 2017 Siwei
# modified 26 Oct 2017 Siwei
# modified 1 Nov 2017 Siwei
# modified 9 Nov 2017 Siwei
# modified 9 Mar 2018
# vcf_suffix='_995.vcf'

date > log.txt

for EACHFILE in pre_bams/*.bam
do
echo $EACHFILE
echo $EACHFILE >> log.txt
echo $date >> log.txt
#############call variants##############
	
	# samtools index -@ 23 $EACHFILE
	~/2TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R ~/2TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		--dbsnp ~/2TB/Databases/hg38/INDEX/dbsnp_146.hg38.vcf \
		-gt_mode DISCOVERY \
		-nct 23 \
		-stand_call_conf 30 \
		-I $EACHFILE \
		-o ${EACHFILE/%.bam/_nofilter.vcf}
done

