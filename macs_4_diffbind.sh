#!/bin/sh

# Siwei 5 Jun 2018

for EACHFILE in *.bam
do
	echo $EACHFILE
	echo ${EACHFILE/%_WASPed.bam/}
	macs2 callpeak -t $EACHFILE -n ${EACHFILE/%_WASPed.bam/} --outdir DiffBind_output/ -f BAMPE -g 3e9 --nomodel --keep-dup all --call-summits -B --verbose 3
done

