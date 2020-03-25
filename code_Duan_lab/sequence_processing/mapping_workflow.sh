#!/bin/sh

# Code for cpgmac

# Set these environment vars to point to
# your local installation of WASP
WASP=$HOME/2TB/Tools/WASP-master
DATA_DIR=~/SSD/WASP_temp
OUTPUT_FINAL=WASPed_BAMs

# These environment vars point to the reference genome and bowtie2.
# in the examples below, the reference genome is assumed
# to be indexed for use with bowtie2
INDEX=$HOME/2TB/Databases/hg38/INDEX/Homo_sapiens_assembly38.fasta

# Start the cycle
for EACHFILE in *.bam
do

# pre-cleanup
	rm -r $DATA_DIR/

# initialise
	mkdir $DATA_DIR/
	echo $EACHFILE
	echo $DATA_DIR/${EACHFILE/%.bam/.fq.gz}
	date


	ls $DATA_DIR/
# Pull out reads that need to be remapped to check for bias
# Use the -p option for paired-end reads.

	python ~/2TB/Tools/WASP-master/mapping/find_intersecting_snps.py --is_paired_end --is_sorted --output_dir $DATA_DIR/ --snp_dir ~/2TB/Databases/hg38/INDEX/SNP_chr_het_14Apr18 $EACHFILE

# Remap the reads, using same the program and options as before.
# NOTE: If you use an option in the first mapping step that modifies the
# reads (e.g. the -5 read trimming option to bowtie2) you should omit this
# option during the second mapping step here (otherwise the reads will be modified
# twice)!         
	echo "map for paired ends"      
        bowtie2 -p 23 -X 2000 --mm --qc-filter --met 1 -t \
                --sensitive --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap \
                -x $INDEX -1 $DATA_DIR/${EACHFILE/%.bam/.remap.fq1.gz} -2 $DATA_DIR/${EACHFILE/%.bam/.remap.fq2.gz} \
                | samtools view -b -u -h -@ 23 | samtools sort -l 9 -@ 23 -m 2G -o $DATA_DIR/${EACHFILE/%.bam/.remap_paired.bam}

# remap to same position
        python $WASP/mapping/filter_remapped_reads.py \
                $DATA_DIR/${EACHFILE/%.bam/.to.remap.bam} \
                $DATA_DIR/${EACHFILE/%.bam/.remap_paired.bam} \
                $DATA_DIR/${EACHFILE/%.bam/.remap_paired.keep.bam}
                ls -lh $DATA_DIR/

################### map and filter for unpaired ones
        
        echo "map for unpaired ends"
        bowtie2 -p 23 -X 2000 --mm --qc-filter --met 1 -t \
                --sensitive --no-mixed --no-discordant --no-dovetail --no-contain --no-overlap \
                -x $INDEX -U $DATA_DIR/${EACHFILE/%.bam/.remap.single.fq.gz} \
                | samtools view -b -u -h -@ 23 | samtools sort -l 9 -@ 23 -m 2G -o $DATA_DIR/${EACHFILE/%.bam/.remap_unpaired.bam}

# remap to same position
	echo "remap"
        python $WASP/mapping/filter_remapped_reads.py \
                $DATA_DIR/${EACHFILE/%.bam/.to.remap.bam} \
                $DATA_DIR/${EACHFILE/%.bam/.remap_unpaired.bam} \
                $DATA_DIR/${EACHFILE/%.bam/.remap_unpaired.keep.bam}
        ls -lh $DATA_DIR/

################################## merge, sort, and dedup

# Create a merged BAM containing [1] reads (PE+SE) that did
# not need remapping [2] filtered remapped reads
	echo "merge bam"
       samtools merge -l 9 -@ 23 -f \
                $DATA_DIR/sim_pe_reads.keep.merged.bam \
                $DATA_DIR/${EACHFILE/%.bam/.remap_paired.keep.bam} \
                $DATA_DIR/${EACHFILE/%.bam/.remap_unpaired.keep.bam} \
                $DATA_DIR/${EACHFILE/%.bam/.keep.bam}

################### Sort and index the bam file
	echo "sort bam"
        samtools sort -l 9 -m 2G -@ 23 \
                -o $DATA_DIR/sim_pe_reads.keep.merged.sorted.bam \
                $DATA_DIR/sim_pe_reads.keep.merged.bam 
        samtools index -@ 23 $DATA_DIR/sim_pe_reads.keep.merged.sorted.bam
	ls -lh $DATA_DIR/

################### Dedup
        ~/2TB/Tools/oracle_JDK/bin/java -Xmx40g -jar ~/2TB/Tools/picard291.jar MarkDuplicates \
		I=$DATA_DIR/sim_pe_reads.keep.merged.sorted.bam \
		O=$DATA_DIR/merged_sorted_dedup.bam M=metrics.txt \
                REMOVE_DUPLICATES=True
	ls -lh $DATA_DIR/

################## Re-add readgroup     
        ~/2TB/Tools/oracle_JDK/bin/java -Xmx40g -jar ~/2TB/Tools/picard291.jar AddOrReplaceReadGroups \
                        INPUT=$DATA_DIR/merged_sorted_dedup.bam \
                        OUTPUT=$OUTPUT_FINAL/${EACHFILE/%.bam/_WASPed.bam} \
                        RGID=${EACHFILE/%.bam/} \
                        RGLB=lib1 \
                        RGPL=illumina \
                        RGPU=unit1 \
                        RGSM=${EACHFILE/%.bam/}

	ls -lh $DATA_DIR/
        ls -lh $OUTPUT_FINAL/

################################## clean up the temp folder

#        rm -r $DATA_DIR/


done
