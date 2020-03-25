#!/bin/sh

#Siwei 3 Apr 2018
# Trimmomatic ILLUMINACLIP:NexteraPE-PE.fa:2:30:7 SLIDINGWINDOW:3:18 MINLENGTH:26
# 2 seed mismatch
# Palindrome clip threshold 30
# Simple clip threshold 7
# Sliding window 3 bp avg quality 18
# Minimum read length 26 bp

for EACHFILE in *R1_001.fastq.gz
do
	echo "$EACHFILE"
	~/2TB/Tools/oracle_JDK/bin/java -Xmx60g -jar ~/2TB/Tools/Trimmomatic-0.36/trimmomatic-0.36.jar \
		PE -threads 23 -basein $EACHFILE -baseout Trimmo/$EACHFILE \
		ILLUMINACLIP:/home/cpg/2TB/Tools/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:7 SLIDINGWINDOW:3:18 MINLEN:26
done

