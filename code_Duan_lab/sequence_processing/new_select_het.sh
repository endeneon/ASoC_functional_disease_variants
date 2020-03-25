#!/bin/sh

suffix="_het.vcf"

for EACHFILE in *.vcf
do
	echo $EACHFILE
#	cat $EACHFILE | grep "^#" > header.txt
#	cat $EACHFILE | grep -v "1/1" > body.txt
#	cat header.txt body.txt > output/$EACHFILE$suffix
	cat $EACHFILE | grep -v "1/1" > output/$EACHFILE$suffix
	wc -l output/$EACHFILE$suffix
done

