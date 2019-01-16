#!/bin/bash

home_folder=/BioII/lulab_b/songyabing/motif_analysis/motif_analysis
# input:get genome info from R
genome_info_folder=/BioII/lulab_b/songyabing/motif_analysis/motif_analysis/genome_info
# input:get interested genes form upstream analysis
interested_gene_folder=/BioII/lulab_b/songyabing/motif_analysis/motif_analysis/interested_genes
# input:gtf file used 
gtf_file=/BioII/lulab_b/songyabing/genome/gencode.v27.annotation.gtf
# input:genome file
genome=/BioII/lulab_b/songyabing/genome/GRCh38.p10.genome.fa

#bin:scripts used
scripts=/BioII/lulab_b/songyabing/motif_analysis/motif_analysis/scripts

# output:get interested genome info 
interested_info_folder=/BioII/lulab_b/songyabing/motif_analysis/motif_analysis/interested_genome_info

for j in `ls $interested_gene_folder`
do 
	sample_folder=$interested_gene_folder/$j
	interested_sample_folder=$interested_info_folder/$j
	mkdir -p $interested_sample_folder
  for i in `ls $sample_folder`
  do 
    # file not empty
		if [ -s $sample_folder/$i ]
		then 
			echo -e "==start analysis of file $i==" 
			echo -e "\n\n==$i|step1:get genome info==\ndone before"
			
			echo -e "\n\n==$i|step2:get interested genome info=="
			sample_prefix=${i%.*}
			echo -e "\n$i:get interested three prime info to interested_sample_folder"
		  sort -t $'\t' -k 2 $genome_info_folder/utr3p.info|join -o 1.3 2.1 1.2 1.9 1.4 1.5 1.6 1.7 1.8 1.10 -t $'\t' -1 2 -2 2 - \
		    <(cut -f 1 $sample_folder/$i |sort |join -t $'\t' -1 1 -2 1 - \
		    <(sort -t $'\t' -k 1 <(grep -o -P -e "gene_id.*; transcript_id.*?;" $gtf_file |\
		    sort |uniq|sed -e 's/gene_id "//' -e 's/"; transcript_id "/\t/' -e 's/";//'))|sort -t $'\t' -k 2  ) \
		    |sort -t $'\t' -k 1 >$interested_sample_folder/${sample_prefix}_3_UTR.info
		  
		  echo -e "\n$i:get interested promoter info to interested_sample_folder"
		  sort -t $'\t' -k 7 $genome_info_folder/promoter.info|join -o 1.1 2.1 1.7 1.2  1.3 1.4 1.5 1.6 -t $'\t' -1 7 -2 2 -  \
		    <(cut -f 1 $sample_folder/$i |sort |join -t $'\t' -1 1 -2 1 - \
		    <(sort -t $'\t' -k 1 <(grep -o -P -e "gene_id.*; transcript_id.*?;" $gtf_file \
		    |sort |uniq|sed -e 's/gene_id "//' -e 's/"; transcript_id "/\t/' -e 's/";//' ))|\
		    sort -t $'\t' -k 2  ) |sort -t $'\t' -k 1 >$interested_sample_folder/${sample_prefix}_promoter.info
		  
		  
		  echo -e "\n\n==$i|step3:convert to bed format=="
		  echo -e "\n$i:convert interested three prime info to bed format"
		  cat $interested_sample_folder/${sample_prefix}_3_UTR.info \
		  	| awk '{print $1 "\t" $5-1 "\t" $6 "\t" $3 "\t" $2 "\t" $8}' | sort -u  > $interested_sample_folder/${sample_prefix}_3_UTR.bed
		  echo -e "\n$i:convert interested promoter info to bed format"
		  cat $interested_sample_folder/${sample_prefix}_promoter.info\
		  	| awk '{print $1 "\t" $4-1 "\t" $5 "\t" $3 "\t" $2 "\t" $7}' | sort -u  > $interested_sample_folder/${sample_prefix}_promoter.bed
		  	
		  	
		  echo -e "\n\n==$i|step4:get genome sequence=="
		  echo -e "\n$i:get interested 3_UTR related genome sequence"
		  bedtools getfasta -s -name -fi $genome -bed $interested_sample_folder/${sample_prefix}_3_UTR.bed \
		  	-fo $interested_sample_folder/${sample_prefix}_3_UTR.fa
		  echo -e "\n$i:concatenate sequences of the same 3¡¯ UTR"
		  Rscript $scripts/concat.r $interested_sample_folder/${sample_prefix}_3_UTR.fa
		  
		  echo -e "\n$i:get interested promoter related genome sequence"
		  bedtools getfasta -s -name -fi $genome -bed $interested_sample_folder/${sample_prefix}_promoter.bed \
		  	-fo $interested_sample_folder/${sample_prefix}_promoter.fa
		  echo -e "\n$i:concatenate sequences of the same promoter"
		  Rscript $scripts/concat.r $interested_sample_folder/${sample_prefix}_promoter.fa
		  
		  
		  echo -e "\n\n==$i|step5:get bg sequence=="
		  echo -e "\n$i:get random 3_UTR sequence"
		  # exits some problems :maybe overlap(before concat UTR sequence)
		  #Rscript $scripts/bg_down.r $interested_sample_folder/${sample_prefix}_3_UTR.bed  $interested_sample_folder/${sample_prefix}_3_UTR_bg.bed 
		  #bedtools getfasta -s -name -fi $genome -bed $interested_sample_folder/${sample_prefix}_3_UTR_bg.bed  \
		  #	-fo $interested_sample_folder/${sample_prefix}_3_UTR_bg.fa	
		  fasta-shuffle-letters $interested_sample_folder/${sample_prefix}_3_UTR.fa $interested_sample_folder/${sample_prefix}_3_UTR.control
		  echo -e "\n$i:get random promoter sequence"
		  fasta-shuffle-letters $interested_sample_folder/${sample_prefix}_promoter.fa $interested_sample_folder/${sample_prefix}_promoter.control
		  
		  
  done
done