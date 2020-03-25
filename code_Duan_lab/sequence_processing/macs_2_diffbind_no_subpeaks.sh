#!/bin/bash

# Siwei 5 Jun 2018
# Siwei 31 Oct 2018


for EACHFILE in *.bam
do
	echo $EACHFILE
	echo ${EACHFILE/%_WASPed.bam/}
	macs2 callpeak -t $EACHFILE \
		-n ${EACHFILE/%_WASPed.bam/} \
		--outdir MACS2_output_no_subpeaks/ \
		--nomodel \
		-f BAMPE \
		-g 3.2e9 \
		--keep-dup all \
		--verbose 3
done

