#!/bin/bash

# input:motif db
motif_db=/BioII/lulab_b/songyabing/motif_analysis/motif_db/motif_databases
# input:get interested genome info (output folder of motif.sh)
interested_info_folder=/BioII/lulab_b/songyabing/motif_analysis/motif_analysis/interested_genome_info

#bin:scripts used
scripts=/BioII/lulab_b/songyabing/motif_analysis/motif_analysis/scripts

# output:result 
results_folder=/BioII/lulab_b/songyabing/motif_analysis/motif_analysis/results

# to improve running efficiency 
for i in `ls $interested_info_folder`
do  
	sample_de_novo_folder=$results_folder/meme_de_novo/$i
	mkdir -p $sample_de_novo_folder
	sample_enrich_motif=$results_folder/ame_enrich/$i
	mkdir -p $sample_enrich_motif
  prefix_list=$(for j in `ls $interested_info_folder/$i`;do name=`echo $j|sed 's/_promoter.*//; s/_3_UTR.*//'`; echo $name;done |uniq)
	for z in `echo $prefix_list`
	do
		meme -dna -maxsize 1000000000 -minw 4 -maxw 12 -oc $sample_de_novo_folder/${z}_UTR_de_novo -mod zoops \
			-nmotifs 5 $interested_info_folder/$i/${z}_3_UTR.fa
		
		meme -dna -maxsize 1000000000 -minw 4 -maxw 12 -oc $sample_de_novo_folder/${z}_promoter_de_novo -mod zoops \
			-nmotifs 5 $interested_info_folder/$i/${z}_promoter.fa
	  
	  ame --control $interested_info_folder/$i/${z}_3_UTR.control --oc $sample_enrich_motif/${z}_UTR_enrich $interested_info_folder/$i/${z}_3_UTR.fa \
	  	$motif_db/CISBP-RNA/Homo_sapiens.meme $motif_db/RNA/Ray2013_rbp_Homo_sapiens.meme $sample_de_novo_folder/${z}_UTR_de_novo/meme.txt
	  
	  ame --control $interested_info_folder/$i/${z}_promoter.control --oc $sample_enrich_motif/${z}_promoter_enrich $interested_info_folder/$i/${z}_promoter.fa \
	  	$motif_db/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme \
	  	$motif_db/CIS-BP/Homo_sapiens.meme $motif_db/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme \
	  	$sample_de_novo_folder/${z}_promoter_de_novo/meme.txt
	done
done
