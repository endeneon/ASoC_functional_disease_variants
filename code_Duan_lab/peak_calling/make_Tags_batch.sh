#!/bin/sh

####Siwei 5 Dec 2017

####read each bam file and make a single tag file each

for EACHFILE in *.bam
do
	date
	echo $EACHFILE
	makeTagDirectory tags/$EACHFILE -single $EACHFILE
done

