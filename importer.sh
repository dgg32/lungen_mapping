#!/bin/bash

script_folder=$(pwd)

inputfolder=$1
negativefolder=$2

projectname=$(basename "${inputfolder}")
negativeproject=$(basename "${negativefolder}")

cpu=28

R1s=$(find $inputfolder -maxdepth 1 -regex '\S+fastq'  -not -name "*trimmed*" | grep "R1")

R1_list=""
R2_list=""
total_list=""

for R1 in $R1s;
do
	R2=${R1/R1/"R2"}
	#printf "$R1\n$R2\n"
	R1trim=${R1/.fastq/"_trimmed.fastq"}
	R2trim=${R2/.fastq/"_trimmed.fastq"}
	single=${R1/.fastq/"_single_trimmed.fastq"}
	single=${single/R1/""}

	R1_no_carp=${R1/.fastq/"_trimmed_no_carp.fastq"}
	R2_no_carp=${R2/.fastq/"_trimmed_no_carp.fastq"}

	# R1_no_conta=${R1/.fastq/"_trimmed_no_carp_no_conta.fastq"}
	# R2_no_conta=${R2/.fastq/"_trimmed_no_carp_no_conta.fastq"}

	# mapping_bam=${R1/.fastq/"_trimmed_no_carp.bam"}
	# bothReadsUnmapped=${R1/.fastq/"_trimmed_no_carp_bothReadsUnmapped.bam"}
	# sort_bothReadsUnmapped=${R1/.fastq/"_trimmed_no_carp_bothReadsUnmapped_sorted.bam"}

	if [ ! -f $R1trim ];
	then
		sickle pe -f $R1 -r $R2 -o $R1trim -p $R2trim -s $single --quiet -t sanger
	fi


	if [ ! -f $R1_no_carp ];
	then
		python filter_fasta_multiprocess.py './taboo.txt' $R1trim $cpu
	fi


	R1_list="${R1_list},$R1_no_carp"
	R2_list="${R2_list},$R2_no_carp"


done

R1_list="${R1_list:1}"
R2_list="${R2_list:1}"


megahitout="${negativefolder}/${negativeproject}_mix_megahit"




bamoutfolder="${negativefolder}/${projectname}_positive_mix_bam"





if [ ! -d $bamoutfolder ];
then
mkdir $bamoutfolder

bowtie2 -x $megahitout/final.contigs -1 $R1_list  -2 $R2_list  | samtools view -bS -o $bamoutfolder/mapping.bam

samtools sort -o $bamoutfolder/mapping_sort.bam $bamoutfolder/mapping.bam  && samtools index $bamoutfolder/mapping_sort.bam

fi



