#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"*reversed.2.cf\"  \"reversed.2.maturefam5p.s\"  \"reversed.2.maturefam3p.s\"  \"output name\" " if (@ARGV ne 4);
my $filein=$ARGV[0];    # _reversed.2.cf
my $mature5p=$ARGV[1];  # _reversed.2.maturefam5ps
my $mature3p=$ARGV[2];  # _reversed.2.maturefam3ps
my $fileout=$ARGV[3];   # _reversed.2.cfs
open IN,$filein;
open IN2,$mature5p;
open IN3,$mature3p;
open OUT,">".$fileout;
my %mat5p;
my %mat3p;
while(<IN2>) {
    chomp;
    my @a=split("\t",$_);
    for (my $i=1; $i<scalar(@a); $i++) {
        $mat5p{$a[$i]}=$a[0];
    }
}
close IN2;
while(<IN3>) {
    chomp;
    my @a=split("\t",$_);
    for (my $i=1; $i<scalar(@a); $i++) {
        $mat3p{$a[$i]}=$a[0];
    }
}
close IN3;
while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    #if (($a[19] ne 0) and (($a[1] - $a[17] - $a[18]) > 8)){next;}
    my $f=0;
    if ($a[30] ne "NA") {
        # 5p mature exists
        if (exists $mat5p{$a[30]}) { $a[30]=$mat5p{$a[30]}; $f+=1; }
        else {$f=$f - 1;}
    }
    if ($a[31] ne "NA") {
        # 3p mature exists
        if (exists $mat3p{$a[31]}) { $a[31]=$mat3p{$a[31]}; $f+=10; }
        else {$f=$f - 10;}
    }
    if (($f eq 1) or ($f eq 10) or ($f eq 11)) {
        print OUT join("\t",@a),"\n";
    }
}