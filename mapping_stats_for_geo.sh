#!/bin/bash

WORK_DIR=$(dirname $PWD)
BAM_FILES="${WORK_DIR}/mapping2/*.sorted.bam"
READ1_FILES="${WORK_DIR}/data/*R1*.gz"

printf "file name 1\tfile name 2\taverage insert size\tstandard deviation\n" >| mapping_stats_for_geo.txt


for READ1 in ${READ1_FILES} 
	do
		prefix=$(echo $READ1| cut -d'_' -f 1)
		prefix=${prefix##*/}
		READ1=${READ1##*/}
		READ2=$(echo ${READ1/R1/R2})
		BAM_FILE="${prefix}.sorted.bam"
		
		$(java -jar /opt/feifei/bin/picard.jar CollectInsertSizeMetrics I=${BAM_FILE} O=${prefix}_insert_size_metrics.txt H=${prefix}_insert_size_histogram.pdf M=0.5)
		insert_size=$(head -8 ${prefix}_insert_size_metrics.txt | tail -1 |cut  -f6)
		sd=$(head -8 ${prefix}_insert_size_metrics.txt | tail -1 |cut  -f7)
		printf "$READ1\t$READ2\t" >>mapping_stats_for_geo.txt
		$(export LC_NUMERIC="en_US.UTF-8")
		printf "%.0f\t" ${insert_size} >>mapping_stats_for_geo.txt
		printf "%.0f\n" ${sd} >>mapping_stats_for_geo.txt

done

