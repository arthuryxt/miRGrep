#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"mapped file\"  \"long reads fastq file\"  \"output name\" " if (@ARGV ne 3);
my %uniq;
my $filename=$ARGV[0];  
my $db=$ARGV[1];       
my %db_Seq;
my %db_Qua;
my $fileout=$ARGV[2];
open OUT,">".$fileout;
open IN1,$db;
my $isfasta=0;
while(<IN1>) {
    chomp;
    if (m/^@/) {
        s/^@//;
        my $id=$_;
        my $seq=<IN1>;
        chomp $seq;
        my $qua=<IN1>;
        $qua=<IN1>;
        chomp $qua;
        $db_Seq{$id}=$seq;
        $db_Qua{$id}=$qua;
        $isfasta=1;
    }
    else {
        s/^>//;
        my $id=$_;
        my $seq=<IN1>;
        chomp $seq;
        $db_Seq{$id}=$seq;
    }
}
close IN1;
print "fastq database read! \n";
open IN,$filename;
while(<IN>) {
    chomp $_;
    my @array=split("\t",$_);
    if (($array[6] ne "+") or (scalar(@array) ne 10) or ($array[1]=~m/\./)) {next}
    my $info=join("\|",$array[0],$array[5],$array[8]);
    if (exists $uniq{$array[7]}) {
        $uniq{$array[7]}=$uniq{$array[7]}."#".$info;
    }
    else {
        $uniq{$array[7]}=$info;
    }
}
close IN;
foreach my $id (keys %uniq) {
    my @tmp=split(/\#/,$uniq{$id});
    my %checkuniq;  # check if one short can be mapped to the same long twice, which is wrong
    my $f=0;
    for(@tmp) {
        my @b=split(/\|/,$_);
        if (exists $checkuniq{$b[0]}) {$f++;}
        else {$checkuniq{$b[0]}=1;}
    }
    if ($f eq 0) {
        if ($isfasta eq 0) { print OUT join("\t",$id,$db_Seq{$id},"+",scalar(@tmp),$uniq{$id}),"\n"; }
        else { print OUT join("\t",$id,$db_Seq{$id},$db_Qua{$id},scalar(@tmp),$uniq{$id}),"\n"; }
    }
}
