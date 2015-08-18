=== 0. Note ===
1. miRGrep stand for "miRNA Genome Reference Free Prediction".
2. It is a experimental/computational pipeline that utilize RNA-Seq technology to discover miRNAs and their precursors, INDEPENDENT of genomic reference sequences.
3. For non-commercial purposes the code is released under the General Public License (version 3). See the file LICENSE for more detail. If you are interested in commercial use of this software contact the author.
4. Please kindly cite doi: 10.1093/nar/gkt072 if you use this tool.
5. Contact: arthur.yxt@gmail.com
=== 0. Note ===




=== 1. Tutorial ===
# make pipeline using the toy data provided in the sub-folder ./data/
cd ./data/
perl ../miRGrep_make.pl SPEC.txt testRUN.sh

# run pipeline
bash ./testRUN.sh

# eye-inspect results
Candidate_premiRNA.display.txt
and
Candidate_premiRNA.display.txt.MIR (if there is high similarity between long_reads and provided hairpin sequences)

Here you have one block for each candidate pre-miRNAs. Within each block, there are:
1. fasta header of the long read
2. sequence of the long read
3. predicted secondary structure of the long read
4. multiple lines of short reads that aligned to the long read, with the occurance of each short read showing to the right

For example:
Truseq_HWI-ST540_124_8_1209_5476_21869__2
AGCTGGTTCCTGCCACGGTGATAATCCAATCTATTGTGACTATCACCGGGTGTGAATCAGTTCGG
(((((((((.((((.((((((((.(((((....))).)).)))))))))))).)))))))))... (-29.70)
AGCTGGTTCCTGCCACGGTGATA                                                 80
AGCTGGTTCCTGCCACGGTGATAA                                                4
                       ATCCAATCTATTGTGACTA                              6
                       ATCCAATCTATTGTGACTAT                             1
                                          TCACCGGGTGTGAATCAGTTCG        84633
                                          TCACCGGGTGTGAATCAGTTCGG       222

Pipeline tested under : Linux 2.6.38-15-server #66~lucid1-Ubuntu SMP
=== 1. Tutorial ===






=== 2. parameters of miRGrep ===
# Note that "folder-to-miRGrep" is where you put miRGrep, such as : /data/analysis/miRGrep/
Below are mandatory parameters, one has to make sure the paths are correct
perl -e 'print join("\t","miRGrep","/folder-to-miRGrep/"),"\n";'                             >  SPEC.txt
perl -e 'print join("\t","soap1","/folder-to-miRGrep/soap.short"),"\n";'                     >> SPEC.txt
perl -e 'print join("\t","BLAT","/folder-to-miRGrep/blat"),"\n";'                            >> SPEC.txt
perl -e 'print join("\t","RNAfold","/folder-to-miRGrep/RNAfold"),"\n";'                      >> SPEC.txt
perl -e 'print join("\t","randfold","/folder-to-miRGrep/randfold"),"\n";'                    >> SPEC.txt
perl -e 'print join("\t","short","/folder-to-miRGrep/data/short.fq"),"\n";'                   >> SPEC.txt
perl -e 'print join("\t","long","/folder-to-miRGrep/data/long.fa"),"\n";'                     >> SPEC.txt
perl -e 'print join("\t","knownMIRscore","/folder-to-miRGrep/data/prec_mmu17_score"),"\n";'   >> SPEC.txt
perl -e 'print join("\t","miRBase","/folder-to-miRGrep/data/prec_mmu17.fa"),"\n";'            >> SPEC.txt
perl -e 'print join("\t","mainRefSpecies","mmu"),"\n";'                                       >> SPEC.txt

Below are optional and default parameters, tweak them if needed.
perl -e 'print join("\t","Threads","8"),"\n";'                    >> SPEC.txt
perl -e 'print join("\t","Sharpness","0.75"),"\n";'               >> SPEC.txt
perl -e 'print join("\t","MFEperbase","0.6"),"\n";'               >> SPEC.txt
perl -e 'print join("\t","randfoldpval","0.2"),"\n";'             >> SPEC.txt
perl -e 'print join("\t","likelihood","0.95"),"\n";'              >> SPEC.txt
perl -e 'print join("\t","minreadcount","5"),"\n";'               >> SPEC.txt
perl -e 'print join("\t","randfoldtrails","199"),"\n";'           >> SPEC.txt
=== 2. parameters of miRGrep ===






=== 3. prepare the annotations for miRGrep ===
# in ./data/ folder
# prepare sample long_reads
cut -f1 ../S5_score > S5_score.id
perl ~/commoncode/get_selected_from_pool_fa.pl long.id ../long.fa long.fa

# prepare sample short_reads
perl ~/commoncode/get_selected_from_pool_singleline.pl long.id ../S1_reversed long_short
cut -f5 long_short | sed -e 's/\#/\n/g' | cut -d\| -f1 | sort -u > short.id
perl ~/commoncode/get_selected_from_pool_fq.pl short.id ../short.fq short.fq

# prepare miRBase score
# take mouse miRNA as example: download <hairpin.fa> and <mature.fa> from miRBase and change the "U" to "T".
### IMPORTANT ###
Make sure that <hairpin.fa> and <mature.fa> are strictly in the dual-liner format as:
"
>fasta_header1
ATCGATGC
>fasta_header2
TGCATGCA
"

grep -A 1 'mmu' hairpin.fa | grep -v "^-" > prec_mmu17.fa
RNAfold -p -d2 -noLP < prec_mmu17.fa > prec_mmu17.fa.str
rm -f *ps
grep -A 2 ">" prec_mmu17.fa.str | grep -v '^-' > prec_mmu17.fa.str1

randfold -d prec_mmu17.fa 199 > prec_mmu17.fa.rand

grep -A 1 'mmu' mature.fa | grep -v "^-" > mature_mmu17.fa
soap.short -a mature_mmu17.fa -d prec_mmu17.fa -o miRBAse_short_long -s 6 -v 0 -n 1 -r 2
perl S1_reversed.pl miRBAse_short_long prec_mmu17.fa S1_miRBase_reversed
perl S1_long_sharpness.pl S1_miRBase_reversed S1_miRBase_reversed
perl S4_scoring.pl prec_mmu17.fa.rand prec_mmu17.fa.str1 long_miRBase S1_miRBase_reversed.2 prec_mmu17_scoring
perl -ne 'chomp; my @a=split("\t",$_); if(scalar(@a) > 20){ print join("\t",@a),"\n"; }' pre_mmu17_scoring > prec_mmu17_score

Now these two files <pre_mmu17.fa> and <prec_mmu17_score> should be written into the SPEC file to generate the miRGrep pipeline
=== 3. prepare the annotations for miRGrep ===










=== 4. Tips ===
The runtime of randfold would be prohibitively long, if given over thousands of long_reads to do through.
To speed up, one can run Step1-Step3 as pipelined, isolate Step4 and run in parallel (if SEG is supported), and finally run Step5.
In stead of running this :
"
randfold -d S1_all_long.fa 199 > S1_all_long.rand
"
run the following instead:
"
split -d -a 3 -l 2000 S1_all_long.fa SEPE
ls -htrml SEPE*  | perl -ne 'chomp;  my @a=split(" ",$_); print $a[-1],"\n"' > list.sepe
perl batch_randfold.pl list.sepe /folder-to-miRGrep/randfold 199
cat SEPE*.rand > S1_all_long.rand
"
=== 4. Tips ===
