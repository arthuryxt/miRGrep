#!/usr/bin/perl -w
use strict;
my $filein=$ARGV[0];    # interesting list
my $DB=$ARGV[1];        # database file, all entry
my $fileout=$ARGV[2];   # resulting interesting entry
open IN,$filein;
my %uniq;
while(<IN>) {
    chomp $_;
    $uniq{$_}=1;
}
close IN;
open IN1,$DB;
open OUT1,">".$fileout;
while(<IN1>) {
    chomp $_;
    s/^>//;
    my $id=$_;
    my $seq=<IN1>;
    if (exists $uniq{$id}) {
        print OUT1 ">".$id,"\n",$seq;
    }
}
