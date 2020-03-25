#!/bin/sh

suffix="_SNP.vcf"

for EACHFILE in *.vcf
do
	~/2TB/Tools/oracle_JDK/bin/java -jar ~/2TB/Tools/GATK37/GenomeAnalysisTK.jar \
		-R ~/2TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta \
		-T SelectVariants \
		-o output/$EACHFILE$suffix \
		-env -selectType SNP \
		-restrictAllelesTo BIALLELIC \
		-V $EACHFILE
done

