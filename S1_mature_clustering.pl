#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"*maturefamily\"  \"output name\" " if (@ARGV ne 2);
my $filein=$ARGV[0];    # S1_reversed.2.maturefam5p or 3p
my $fileout=$ARGV[1];   # S1_reversed.2.maturefam5p.s
open IN,$filein;
open OUT,">".$fileout;
my %reversed;
my $count=0;
while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    my @info=split(/\_\_/,$a[0]);
    if (exists $reversed{$a[0]}) {
        for(my $i=1; $i<scalar(@a); $i++) {
            $reversed{$a[$i]}=$reversed{$a[0]};
        }
    }
    else {
        for(my $i=1; $i<scalar(@a); $i++) {
            $reversed{$a[$i]}=$count;
        }
        $count++;
    }
}
print $count,"\n";
my %uniq;
foreach my $id(keys %reversed) {
    if (exists $uniq{$reversed{$id}}) {
        if ($id =~ m/$uniq{$reversed{$id}}/) {}
        else {$uniq{$reversed{$id}}=$uniq{$reversed{$id}}."\t".$id;}
    }
    else {
        $uniq{$reversed{$id}}=$id;
    }
}
foreach my $id(keys %uniq) {
    my @a=split("\t",$uniq{$id});
    my $max=0;
    my $maxid='';
    foreach(@a) {
        my @b=split(/\_\_/,$_);
        if ($b[1] > $max) {
            $max=$b[1];
            $maxid=join("\_\_",@b);
        }
    }
    print OUT $maxid,"\t",$uniq{$id},"\n";
}