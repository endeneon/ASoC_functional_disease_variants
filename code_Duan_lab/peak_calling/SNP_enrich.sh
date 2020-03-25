#!/bin/sh

# Siwei 24 May 2018
# Siwei 21 Nov 2018

for EACHFILE in *.narrowPeak
do
	echo $EACHFILE
	echo ${EACHFILE/%_peaks.narrowPeak/}
	mkdir ${EACHFILE/%_peaks.narrowPeak/}
	findMotifsGenome.pl $EACHFILE hg38 ${EACHFILE/%_peaks.narrowPeak/} -size 100 -mask -cache 30000
done

