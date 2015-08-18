#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"config_file\"   \"output_sh\"  " if (@ARGV < 2);
my $filein=$ARGV[0];
my $fileout=$ARGV[1];
my %SPEC;
open(IN,$filein) or die "Cannot open config_file $filein";
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
    if ((scalar(@a) < 2) or ($a[0] eq "")) { next; }
    $SPEC{$a[0]}=$a[1];
}
close IN;
my $thread=8;
my $sharpness=0.75;
my $mfepb=0.6;
my $randfoldpval=0.2;
my $likelihood=0.95;
my $minreadcount=5;
my $randfoldtrails=199;
my $species="mmu";
my $maplong="no";
# check if all parameters are set
if (!exists $SPEC{"short"}) { die "short reads fastq must by specified in the config_file $filein";}
if (!exists $SPEC{"long"}) { die "long reads fasta must by specified in the config_file $filein";}
if (!exists $SPEC{"miRGrep"}) { die "miRGrep_folder must by specified in the config_file $filein";}
if (!exists $SPEC{"soap1"}) { die "executable soap.short must by specified in the config_file $filein";}
if (!exists $SPEC{"BLAT"}) { die "executable blat must by specified in the config_file $filein";}
if (!exists $SPEC{"RNAfold"}) { die "executable RNAfold must by specified in the config_file $filein";}
if (!exists $SPEC{"randfold"}) { die "executable randfold must by specified in the config_file $filein";}
if (!exists $SPEC{"knownMIRscore"}) { die "score for known miRNAs must by specified in the config_file $filein";}
if (!exists $SPEC{"miRBase"}) { die "fasta file for reference pre-miRNA hairpins must by specified in the config_file $filein";}
if (exists $SPEC{"Thread"}) { $thread=$SPEC{"Thread"}; }
if (exists $SPEC{"Sharpness"}) { $sharpness=$SPEC{"Sharpness"}; }
if (exists $SPEC{"MFEperbase"}) { $mfepb=$SPEC{"MFEperbase"}; }
if (exists $SPEC{"randfoldpval"}) { $randfoldpval=$SPEC{"randfoldpval"}; }
if (exists $SPEC{"likelihood"}) { $likelihood=$SPEC{"likelihood"}; }
if (exists $SPEC{"minreadcount"}) { $minreadcount=$SPEC{"minreadcount"}; }
if (exists $SPEC{"randfoldtrails"}) { $randfoldtrails=$SPEC{"randfoldtrails"}; }
if (exists $SPEC{"mainRefSpecies"}) { $species=$SPEC{"mainRefSpecies"}; }
if (exists $SPEC{"maplong"}) { $maplong=$SPEC{"maplong"}; }

my $command="";
open(OUT, ">".$fileout) or die "Cannot open output_sh file $fileout";
print OUT "#!/bin/bash\n\n";

print OUT "#Step1\n";
print OUT "echo \"Step1 mapping short reads \(mature miRNAs\) against long reads \(pre-miRNAs\) \" \n";
$command=$SPEC{"soap1"}." -a ".$SPEC{"short"}." -d ".$SPEC{"miRBase"}." -o short_map2_miRBase -s 6 -v 2 -n 1 -r 1 -p ".$thread;
print OUT $command,"\n";
if ($maplong ne "no") {
    $command=$SPEC{"soap1"}." -a ".$SPEC{"long"}." -d ".$SPEC{"miRBase"}." -o long_map2_miRBase -s 6 -v 2 -n 1 -r 1 -p ".$thread;
    print OUT $command,"\n";
}
else {
    $command="touch long_map2_miRBase";
    print OUT $command,"\n";
}
$command=$SPEC{"soap1"}." -a ".$SPEC{"short"}." -d ".$SPEC{"long"}." -o short_map2_long -s 6 -v 0 -n 1 -r 2 -p ".$thread;
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/S1_reversed.pl short_map2_long ".$SPEC{"long"}." S1_reversed";
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/S1_long_sharpness.pl S1_reversed S1_reversed";
print OUT $command,"\n";
print OUT "echo \"Step1 mapping short reads against long reads Finished\" \n\n\n";


print OUT "#Step2\n";
print OUT "echo \"Step2 filtering potential pre-miRNAs Started\" \n";
$command="perl ".$SPEC{"miRGrep"}."/S1_mature_clustering.pl S1_reversed.2.maturefam5p S1_reversed.2.maturefam5p.s";
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/S1_mature_clustering.pl S1_reversed.2.maturefam3p S1_reversed.2.maturefam3p.s";
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/S1_long_selecting.pl S1_reversed.2.cf S1_reversed.2.maturefam5p.s S1_reversed.2.maturefam3p.s S1_reversed.2.cfs";
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/S1_long_clustering.pl S1_reversed.2.cfs ".$SPEC{"long"}." S1_reversed.2.cfsc";
print OUT $command,"\n";
$command="sort -u S1_reversed.2.cfsc \| perl -ne \'chomp; \@a=split(\"\\t\",\$_); if((\$a[16]+\$a[29]) >= ".$sharpness."){print join(\"\\t\",\@a),\"\\n\";}\' \> S1_reversed.2.cfsc.S75";
print OUT $command,"\n";
$command="cut -f1 S1_reversed.2.cfsc.S75 \> S1_reversed.2.cfsc.S75.id";
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/get_selected_from_pool_fa.pl S1_reversed.2.cfsc.S75.id ".$SPEC{"long"}." S1_reversed.2.cfsc.S75.fa";
print OUT $command,"\n";
$command=$SPEC{"BLAT"}." -minIdentity=90 -stepSize=5 -tileSize=8 -minScore=40 -noHead S1_reversed.2.cfsc.S75.fa S1_reversed.2.cfsc.S75.fa S1_reversed.2.cfsc.S75.psl";
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/S3_condense.pl S1_reversed.2.cfsc.S75 S1_reversed.2.cfsc.S75.psl S1_reversed.2.cfsc.S75.psl.parsed";
print OUT $command,"\n";
$command="cut -f1 S1_reversed.2.cfsc.S75.psl.parsed \| sort -u \> S1_reversed.2.cfsc.S75.psl.parsed.id";
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/get_selected_from_pool_fa.pl S1_reversed.2.cfsc.S75.psl.parsed.id ".$SPEC{"long"}." S1_reversed.2.cfsc.S75.psl.parsed.fa";
print OUT $command,"\n";
$command="cat S1_reversed.2.cfsc.S75.psl.parsed.id \| sort -u \> S1_all_long.id";
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/get_selected_from_pool_fa.pl S1_all_long.id ".$SPEC{"long"}." S1_all_long.fa";
print OUT $command,"\n";
print OUT "echo \"Step2 filtering potential pre-miRNAs Finished\" \n\n\n";


print OUT "#Step3\n";
print OUT "echo \"Step3 estimate hairpin structure Started\" \n";
$command=$SPEC{"RNAfold"}." -p -d2 -noLP \< S1_all_long.fa \> S1_all_long.str";
print OUT $command,"\n";
$command="grep -A 2 \"\>\" S1_all_long.str \| grep -v \'\^-\' \> S1_all_long.str1";
print OUT $command,"\n";
print OUT "echo \"Step3 estimate hairpin structure Finished\" \n\n\n";

print OUT "#Step4\n";
print OUT "echo \"Step4 estimate hairpin stability Started\" \n";
$command=$SPEC{"randfold"}." -d S1_all_long.fa ".$randfoldtrails." > S1_all_long.rand";
print OUT $command,"\n";
print OUT "echo \"Step4 estimate hairpin stability Finished\" \n\n\n";


print OUT "#Step5\n";
print OUT "echo \"Step5 Scoring candidate pre-miRNAs Started\" \n";
$command="perl ".$SPEC{"miRGrep"}."/S4_scoring.pl S1_all_long.rand S1_all_long.str long_map2_miRBase S1_reversed.2.cfsc.S75.psl.parsed S4_scoring ";
print OUT $command,"\n";
$command="perl -ne \'chomp; my \@a=split(\"\\t\",\$_); if(scalar(\@a) > 10){ print join\(\"\\t\",\@a\),\"\\n\";}\' S4_scoring \> S4_score";
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/S5_validating_with_miRBase.pl S1_reversed S4_scoring S4_score short_map2_miRBase S5_full ".$species;
print OUT $command,"\n";
$command="perl -ne \'chomp; my \@a=split(\"\\t\",\$_); if(scalar(\@a) > 40){ print join\(\"\\t\",\@a\),\"\\n\";}\' S5_full \> S5_score";
print OUT $command,"\n";
$command="awk \'\$32 \< ".$randfoldpval." \&\& \$45 \> ".$mfepb." \' S5_score \> S5_score_fg ";
print OUT $command,"\n";
$command="awk \'\$32 \>\= ".$randfoldpval." \|\| \$45 \<\= ".$mfepb." \' S5_score \> S5_score_bg ";
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/S6_Bayesian_scoring.pl ".$SPEC{"knownMIRscore"}." S5_score_bg S5_score_fg S6_score";
print OUT $command,"\n";
$command="awk \'\$7 > ".$likelihood." \&\& \$12 \>\=".$minreadcount."\'  S6_score \| cut -f9 \> Candidate_premiRNAs.id";
print OUT $command,"\n";
$command="perl ".$SPEC{"miRGrep"}."/S0_print.pl S1_reversed ".$SPEC{"long"}." Candidate_premiRNAs.id S1_all_long.str1 short_map2_miRBase long_map2_miRBase Candidate_premiRNA.display.txt";
print OUT $command,"\n";
print OUT "echo \"Step5 Scoring candidate pre-miRNAs Finished\" \n\n\n";



close OUT;
