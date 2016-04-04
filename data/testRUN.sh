#!/bin/bash

#Step1
echo "Step1 mapping short reads (mature miRNAs) against long reads (pre-miRNAs) " 
perl ..//unique_reads_for_mapping.pl ./long.fa 20 1000 1
perl ..//unique_reads_for_mapping.pl ./short.fq 20 1000
../soap.short -a ./short.fq.uniq -d ./prec_mmu17.fa -o short_map2_miRBase -s 6 -v 2 -n 1 -r 1 -p 8
touch long_map2_miRBase
../soap.short -a ./short.fq.uniq -d ./long.fa.uniq -o short_map2_long -s 6 -v 0 -n 1 -r 2 -p 8
perl ..//S1_reversed.pl short_map2_long ./long.fa.uniq S1_reversed
perl ..//S1_long_sharpness.pl S1_reversed S1_reversed
echo "Step1 mapping short reads against long reads Finished" 


#Step2
echo "Step2 filtering potential pre-miRNAs Started" 
perl ..//S1_mature_clustering.pl S1_reversed.2.maturefam5p S1_reversed.2.maturefam5p.s
perl ..//S1_mature_clustering.pl S1_reversed.2.maturefam3p S1_reversed.2.maturefam3p.s
perl ..//S1_long_selecting.pl S1_reversed.2.cf S1_reversed.2.maturefam5p.s S1_reversed.2.maturefam3p.s S1_reversed.2.cfs
perl ..//S1_long_clustering.pl S1_reversed.2.cfs ./long.fa.uniq S1_reversed.2.cfsc
sort -u S1_reversed.2.cfsc | perl -ne 'chomp; @a=split("\t",$_); if(($a[16]+$a[29]) >= 0.75){print join("\t",@a),"\n";}' > S1_reversed.2.cfsc.S75
cut -f1 S1_reversed.2.cfsc.S75 > S1_reversed.2.cfsc.S75.id
perl ..//get_selected_from_pool_fa.pl S1_reversed.2.cfsc.S75.id ./long.fa.uniq S1_reversed.2.cfsc.S75.fa
../blat -minIdentity=90 -stepSize=5 -tileSize=8 -minScore=40 -noHead S1_reversed.2.cfsc.S75.fa S1_reversed.2.cfsc.S75.fa S1_reversed.2.cfsc.S75.psl
perl ..//S3_condense.pl S1_reversed.2.cfsc.S75 S1_reversed.2.cfsc.S75.psl S1_reversed.2.cfsc.S75.psl.parsed
cut -f1 S1_reversed.2.cfsc.S75.psl.parsed | sort -u > S1_reversed.2.cfsc.S75.psl.parsed.id
perl ..//get_selected_from_pool_fa.pl S1_reversed.2.cfsc.S75.psl.parsed.id ./long.fa.uniq S1_reversed.2.cfsc.S75.psl.parsed.fa
cat S1_reversed.2.cfsc.S75.psl.parsed.id | sort -u > S1_all_long.id
perl ..//get_selected_from_pool_fa.pl S1_all_long.id ./long.fa.uniq S1_all_long.fa
echo "Step2 filtering potential pre-miRNAs Finished" 


#Step3
echo "Step3 estimate hairpin structure Started" 
../RNAfold -p -d2 -noLP < S1_all_long.fa > S1_all_long.str
grep -A 2 ">" S1_all_long.str | grep -v '^-' > S1_all_long.str1
echo "Step3 estimate hairpin structure Finished" 


#Step4
echo "Step4 estimate hairpin stability Started" 
../randfold -d S1_all_long.fa 199 > S1_all_long.rand
echo "Step4 estimate hairpin stability Finished" 


#Step5
echo "Step5 Scoring candidate pre-miRNAs Started" 
perl ..//S4_scoring.pl S1_all_long.rand S1_all_long.str long_map2_miRBase S1_reversed.2.cfsc.S75.psl.parsed S4_scoring 
perl -ne 'chomp; my @a=split("\t",$_); if(scalar(@a) > 10){ print join("\t",@a),"\n";}' S4_scoring > S4_score
perl ..//S5_validating_with_miRBase.pl S1_reversed S4_scoring S4_score short_map2_miRBase S5_full mmu
perl -ne 'chomp; my @a=split("\t",$_); if(scalar(@a) > 40){ print join("\t",@a),"\n";}' S5_full > S5_score
awk '$32 < 0.2 && $45 > 0.6 ' S5_score > S5_score_fg 
awk '$32 >= 0.2 || $45 <= 0.6 ' S5_score > S5_score_bg 
perl ..//S6_Bayesian_scoring.pl ./prec_mmu17_score S5_score_bg S5_score_fg S6_score
awk '$7 > 0.95 && $12 >=5'  S6_score | cut -f9 > Candidate_premiRNAs.id
perl ..//S0_print.pl S1_reversed ./long.fa.uniq Candidate_premiRNAs.id S1_all_long.str1 short_map2_miRBase long_map2_miRBase Candidate_premiRNA.display.txt
echo "Step5 Scoring candidate pre-miRNAs Finished" 


