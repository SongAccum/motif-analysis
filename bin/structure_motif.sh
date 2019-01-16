#!/bin/bash

# software folder
beam=/BioII/lulab_b/songyabing/motif_analysis/software/BEAM/beam-2.0

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
	sample_structure_folder=$results_folder/structure/$i
	mkdir -p $sample_structure_folder
	
  prefix_list=$(for j in `ls $interested_info_folder/$i`;do name=`echo $j|sed 's/_promoter.*//; s/_3_UTR.*//'`; echo $name;done |uniq)
	for z in `echo $prefix_list`
	do
		# get fastB format
		for m in "promoter" "3_UTR"
		do
			folder_prefix=$(echo $m |sed 's/3_//')
			diff_type_folder=$sample_structure_folder/${z}_${folder_prefix}_structure
			mkdir -p $diff_type_folder
			# RNAfold produce ps files delete it after RNAfold run finish
			#echo "$interested_info_folder/$i/${z}_${m}.fa"
			RNAfold  < $interested_info_folder/$i/${z}_${m}.fa >$diff_type_folder/${z}_${m}_dot.fa
			RNAfold  < $interested_info_folder/$i/${z}_${m}.control >$diff_type_folder/${z}_${m}_dot.control
			# fastB (fastBEAR)
			awk '/^>/ {print; getline; print; getline; print $1}' $diff_type_folder/${z}_${m}_dot.fa >$diff_type_folder/${z}_${m}_dot_to_encode.fa
			java -jar $beam/BearEncoder.new.jar $diff_type_folder/${z}_${m}_dot_to_encode.fa  $diff_type_folder/${z}_${m}_BEAMready.fa
			
			awk '/^>/ {print; getline; print; getline; print $1}' $diff_type_folder/${z}_${m}_dot.control >$diff_type_folder/${z}_${m}_dot_to_encode.control
			java -jar $beam/BearEncoder.new.jar $diff_type_folder/${z}_${m}_dot_to_encode.control  $diff_type_folder/${z}_${m}_BEAMready.control
			
			# remove all file with ps extention
			rm $beam/*.ps
		done
		
		# run structure motif 
		# get into the folder to output in whole folder
		cd $beam
		
		for t in "promoter" "3_UTR"
		do
			folder_prefix=$(echo $t |sed 's/3_//')
			diff_type_folder=$sample_structure_folder/${z}_${folder_prefix}_structure
			mkdir -p $diff_type_folder
			
			#change file name(replace dot with _)
			new_name=$(echo "$z"|sed 's/\./_/g')
			mv $diff_type_folder/${z}_${t}_BEAMready.fa  $diff_type_folder/${new_name}_${t}_BEAMready.fa
			mv $diff_type_folder/${z}_${t}_BEAMready.control  $diff_type_folder/${new_name}_${t}_BEAMready.control
			
			# change python version used
			PATH=/apps/anaconda2/bin:$PATH
			echo "$(python --version)"
			java -jar $beam/BEAM_release1.6.1.jar -f $diff_type_folder/${new_name}_${t}_BEAMready.fa \
				-g $diff_type_folder/${new_name}_${t}_BEAMready.control -w 10 -W 40 -M 3
			echo "$(python --version)"
			
			# plot weblogo
			logo_motif_folder=$beam/risultati/${new_name}_${t}_BEAMready/webLogoOut/motifs
			for q in `ls $logo_motif_folder`
			do
				weblogo -a 'ZAQXSWCDEVFRBGTNHY' -f $logo_motif_folder/$q -D fasta --resolution 600 \
					-o $logo_motif_folder/${q}.jpeg -F jpeg --composition="none" \
					-C red ZAQ 'Stem' -C blue XSW 'Loop' -C forestgreen CDE 'InternalLoop' \
					-C orange VFR 'StemBranch' -C DarkOrange B 'Bulge' \
					-C lime G 'BulgeBranch' -C purple T 'Branching' \
					-C limegreen NHY 'InternalLoopBranch'
			done
			
			# plot 2D structure
			benchmark_folder=$beam/risultati/${new_name}_${t}_BEAMready/benchmark/motifs
			for w in `ls $benchmark_folder|grep 'search.txt$'`
			do
				#motif_search_file=$benchmark_folder/$w
				python $scripts/visual_2D.py  $benchmark_folder/$w $benchmark_folder
			done
			
		done
		
	done
done
