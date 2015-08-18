#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"fasta_to_run\"  \"executable_to_randrold\"   \"number_of_iterations\" " if (@ARGV < 2);
my $file=$ARGV[0];
my $randfold=$ARGV[1];
my $Nr=199;
if (scalar(@ARGV) > 2) {$Nr=$ARGV[2];}
open IN,$file;
while(my $id=<IN>) {
    chomp $id;
    system("qsub -N $id  -cwd -b y '$randfold -d $id $Nr > $id.rand'");
}