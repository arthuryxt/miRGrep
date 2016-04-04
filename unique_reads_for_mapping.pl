#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"fasta-or-fastq\"  \"\(optional\)min_len\"  \"\(optional\)max_len\"  \(optional\)outputfmt auto=0 fa=1  " if (@ARGV < 1);
my $filein=$ARGV[0];
my $min_len=16;
my $max_len=999999;
if (scalar(@ARGV) > 1) { $min_len = $ARGV[1]; }
if ($min_len < 1) { die "min_len should be >= 1" }
if (scalar(@ARGV) > 2) { $max_len = $ARGV[2]; }
if ($max_len < $min_len) { die "max_len should be >= 1 min_len" }
my $outputfmt=0;
if (scalar(@ARGV) > 3) { $outputfmt = $ARGV[3]; }
open IN,$filein;
open OUT,">".$filein.".uniq";
my $count=0;
my %uniq;
my %qua;
my %ider;
my $filetype;
while(my $id=<IN>) {
    chomp $id;
    if ($id=~m/^@/) {$filetype="fq";}
    elsif($id=~m/^>/) {$filetype="fa";}
    if($filetype eq "fq") {
        $id=~s/^@//;
        my @tmp=split("\_\_",$id);
        my $count=1;
        if (scalar(@tmp) > 1) {$count=$tmp[1];}
        my $seq=<IN>;
        chomp $seq;
        my $seqlen=length($seq);
        if (($seqlen < $min_len ) or ($seqlen > $max_len)) { next; }
        if (exists $uniq{$seq}) {
            $uniq{$seq}+=$count;
            $seq=<IN>;
            $seq=<IN>;
        }
        else {
            $uniq{$seq}=$count;
            $ider{$seq}=$tmp[0];
            my $qua=<IN>;
            $qua=<IN>;
            chomp $qua;
            $qua{$seq}=$qua;
        }
    }
    elsif($filetype eq "fa") {
        $id=~s/^>//;
        my @tmp=split("\_\_",$id);
        my $count=1;
        if (scalar(@tmp) > 1) {$count=$tmp[1];}
        my $seq=<IN>;
        chomp $seq;
        my $seqlen=length($seq);
        if (($seqlen < $min_len ) or ($seqlen > $max_len)) { next; }
        if (exists $uniq{$seq}) {
            $uniq{$seq}+=$count;
        }
        else {
            $uniq{$seq}=$count;
            $ider{$seq}=$tmp[0];
        }
    }    
}
if ($outputfmt eq 1) { $filetype="fa"; }
if($filetype eq "fq") {
    foreach my $seq (sort{$uniq{$b} <=> $uniq{$a}} keys %uniq) {
        my $seqlen=length($seq);
        if (($seqlen < $min_len ) or ($seqlen > $max_len)) { next; }
        print OUT "@".$ider{$seq}."__".$uniq{$seq},"\n";
        print OUT $seq,"\n";
        print OUT "+\n";
        print OUT $qua{$seq},"\n";
    }
}
else {
    foreach my $seq (sort{$uniq{$b} <=> $uniq{$a}} keys %uniq) {
        my $seqlen=length($seq);
        if (($seqlen < $min_len ) or ($seqlen > $max_len)) { next; }
        print OUT ">".$ider{$seq}."__".$uniq{$seq},"\n";
        print OUT $seq,"\n";
    }
}
