#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"miRBase precursor score\"  \"background score\"  \"foreground score\"  \"output\"" if (@ARGV < 4);
my $miRBase=$ARGV[0];   # pre_mmu17_score
my $longlist=$ARGV[1];  # S5_score_bg
my $shortlist=$ARGV[2]; # S5_score_fg
my $fileout=$ARGV[3];
open IN, $miRBase;
open IN2, $longlist;
open IN3, $shortlist;
open OUT,">".$fileout;
my %positive;
my %background;
my $count=0;
while(<IN>) {
    chomp;
    $count++;
    my @a=split("\t",$_);
    my $mfe=round(1-60*$a[30]/$a[1]);
    $positive{"mfe"}{$mfe}++;
    my $stability=int($a[31]*100);
    if ($stability eq 0) {$stability = 1;}
    $positive{"stab"}{$stability}++;
    my $p5=$a[32];
    $positive{"p5"}{$p5}++;
    my $overhang=abs($a[34] - 2);
    $positive{"overhang"}{$overhang}++;
    #my $internal_pair=$a[36];
    #$positive{"int_pair"}{$internal_pair}++;
    #my $basepair=$a[39];
    #$positive{"basepair"}{$basepair}++;
    if ($a[13] ne 0) {$positive{"mat_len"}{$a[13]}++;}
    if ($a[26] ne 0) {$positive{"mat_len"}{$a[26]}++;}
    my $pre_length=$a[1];
    $positive{"pre_len"}{$pre_length}++;
}
# smooth hist to pdf
##### basepair
my $tmp=0;
#for (my $i=0; $i <= 25; $i++) {
#    if (exists $positive{"basepair"}{$i}) { $positive{"basepair"}{$i} = $positive{"basepair"}{$i} / $count;}
#    else { $positive{"basepair"}{$i}=0.5/$count; }
#    if($positive{"basepair"}{$i} > $tmp) {$tmp=$positive{"basepair"}{$i};}
#}
#for (my $i=1; $i < 25; $i++) {
#    if ($positive{"basepair"}{$i} ne $tmp) {
#        my $count=1;
#        my $sum5=$positive{"basepair"}{$i};
#        if (exists $positive{"basepair"}{$i-2}) {$sum5+=$positive{"basepair"}{$i-2}; $count++;}
#        if (exists $positive{"basepair"}{$i-1}) {$sum5+=$positive{"basepair"}{$i-1}; $count++;}
#        if (exists $positive{"basepair"}{$i+1}) {$sum5+=$positive{"basepair"}{$i+1}; $count++;}
#        if (exists $positive{"basepair"}{$i+2}) {$sum5+=$positive{"basepair"}{$i+2}; $count++;}
#        $positive{"basepair"}{$i}=$sum5/5;
#    }
#}
##### int_pair
#$tmp=0;
#for (my $i=0; $i <= 20; $i++) {
#    if (exists $positive{"int_pair"}{$i}) { $positive{"int_pair"}{$i} = $positive{"int_pair"}{$i} / $count;}
#    else { $positive{"int_pair"}{$i}=0.5/$count; }
#    if($positive{"int_pair"}{$i} > $tmp) {$tmp=$positive{"int_pair"}{$i};}
#}
#for (my $i=1; $i < 20; $i++) {
#    if ($positive{"int_pair"}{$i} ne $tmp) {
#        my $count=1;
#        my $sum5=$positive{"int_pair"}{$i};
#        if (exists $positive{"int_pair"}{$i-2}) {$sum5+=$positive{"int_pair"}{$i-2}; $count++;}
#        if (exists $positive{"int_pair"}{$i-1}) {$sum5+=$positive{"int_pair"}{$i-1}; $count++;}
#        if (exists $positive{"int_pair"}{$i+1}) {$sum5+=$positive{"int_pair"}{$i+1}; $count++;}
#        if (exists $positive{"int_pair"}{$i+2}) {$sum5+=$positive{"int_pair"}{$i+2}; $count++;}
#        $positive{"int_pair"}{$i}=$sum5/5;
#    }
#}
##### mfe
#$tmp=0;
for (my $i=0; $i <= 70; $i++) {
    if (exists $positive{"mfe"}{$i}) { $positive{"mfe"}{$i} = $positive{"mfe"}{$i} / $count;}
    else { $positive{"mfe"}{$i}=0.5/$count; }
    if($positive{"mfe"}{$i} > $tmp) {$tmp=$positive{"mfe"}{$i};}
}
for (my $i=1; $i < 70; $i++) {
    if ($positive{"mfe"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$positive{"mfe"}{$i};
        if (exists $positive{"mfe"}{$i-2}) {$sum5+=$positive{"mfe"}{$i-2}; $count++;}
        if (exists $positive{"mfe"}{$i-1}) {$sum5+=$positive{"mfe"}{$i-1}; $count++;}
        if (exists $positive{"mfe"}{$i+1}) {$sum5+=$positive{"mfe"}{$i+1}; $count++;}
        if (exists $positive{"mfe"}{$i+2}) {$sum5+=$positive{"mfe"}{$i+2}; $count++;}
        $positive{"mfe"}{$i}=$sum5/5;
    }
}
##### overhang
$tmp=0;
for (my $i=0; $i <= 50; $i++) {
    if (exists $positive{"overhang"}{$i}) { $positive{"overhang"}{$i} = $positive{"overhang"}{$i} / $count;}
    else { $positive{"overhang"}{$i}=0.5/$count; }
    if($positive{"overhang"}{$i} > $tmp) {$tmp=$positive{"overhang"}{$i};}
}
for (my $i=1; $i < 50; $i++) {
    if ($positive{"overhang"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$positive{"overhang"}{$i};
        if (exists $positive{"overhang"}{$i-2}) {$sum5+=$positive{"overhang"}{$i-2}; $count++;}
        if (exists $positive{"overhang"}{$i-1}) {$sum5+=$positive{"overhang"}{$i-1}; $count++;}
        if (exists $positive{"overhang"}{$i+1}) {$sum5+=$positive{"overhang"}{$i+1}; $count++;}
        if (exists $positive{"overhang"}{$i+2}) {$sum5+=$positive{"overhang"}{$i+2}; $count++;}
        $positive{"overhang"}{$i}=$sum5/5;
    }
}
##### p5
$tmp=0;
for (my $i=0; $i <= 50; $i++) {
    if (exists $positive{"p5"}{$i}) { $positive{"p5"}{$i} = $positive{"p5"}{$i} / $count;}
    else { $positive{"p5"}{$i}=0.5/$count; }
    if($positive{"p5"}{$i} > $tmp) {$tmp=$positive{"p5"}{$i};}
}
for (my $i=1; $i < 50; $i++) {
    if ($positive{"p5"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$positive{"p5"}{$i};
        if (exists $positive{"p5"}{$i-2}) {$sum5+=$positive{"p5"}{$i-2}; $count++;}
        if (exists $positive{"p5"}{$i-1}) {$sum5+=$positive{"p5"}{$i-1}; $count++;}
        if (exists $positive{"p5"}{$i+1}) {$sum5+=$positive{"p5"}{$i+1}; $count++;}
        if (exists $positive{"p5"}{$i+2}) {$sum5+=$positive{"p5"}{$i+2}; $count++;}
        $positive{"p5"}{$i}=$sum5/5;
    }
}
##### stab
$tmp=0;
for (my $i=0; $i <= 100; $i++) {
    if (exists $positive{"stab"}{$i}) { $positive{"stab"}{$i} = $positive{"stab"}{$i} / $count;}
    else { $positive{"stab"}{$i}=0.005/$count; }
    if($positive{"stab"}{$i} > $tmp) {$tmp=$positive{"stab"}{$i};}
}
for (my $i=1; $i < 100; $i++) {
    if ($positive{"stab"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$positive{"stab"}{$i};
        if (exists $positive{"stab"}{$i-2}) {$sum5+=$positive{"stab"}{$i-2}; $count++;}
        if (exists $positive{"stab"}{$i-1}) {$sum5+=$positive{"stab"}{$i-1}; $count++;}
        if (exists $positive{"stab"}{$i+1}) {$sum5+=$positive{"stab"}{$i+1}; $count++;}
        if (exists $positive{"stab"}{$i+2}) {$sum5+=$positive{"stab"}{$i+2}; $count++;}
        $positive{"stab"}{$i}=$sum5/5;
    }
}
##### mat_len
$tmp=0;
for (my $i=15; $i <= 36; $i++) {
    if (exists $positive{"mat_len"}{$i}) { $positive{"mat_len"}{$i} = $positive{"mat_len"}{$i} / $count;}
    else { $positive{"mat_len"}{$i}=0.5/$count; }
    if($positive{"mat_len"}{$i} > $tmp) {$tmp=$positive{"mat_len"}{$i};}
}
for (my $i=16; $i < 36; $i++) {
    if ($positive{"mat_len"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$positive{"mat_len"}{$i};
        if (exists $positive{"mat_len"}{$i-2}) {$sum5+=$positive{"mat_len"}{$i-2}; $count++;}
        if (exists $positive{"mat_len"}{$i-1}) {$sum5+=$positive{"mat_len"}{$i-1}; $count++;}
        if (exists $positive{"mat_len"}{$i+1}) {$sum5+=$positive{"mat_len"}{$i+1}; $count++;}
        if (exists $positive{"mat_len"}{$i+2}) {$sum5+=$positive{"mat_len"}{$i+2}; $count++;}
        $positive{"mat_len"}{$i}=$sum5/5;
    }
}
##### pre_len
$tmp=0;
for (my $i=40; $i <= 100; $i++) {
    if (exists $positive{"pre_len"}{$i}) { $positive{"pre_len"}{$i} = $positive{"pre_len"}{$i} / $count;}
    else { $positive{"pre_len"}{$i}=0.5/$count; }
    if($positive{"pre_len"}{$i} > $tmp) {$tmp=$positive{"pre_len"}{$i};}
}
for (my $i=41; $i < 100; $i++) {
    if ($positive{"pre_len"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$positive{"pre_len"}{$i};
        if (exists $positive{"pre_len"}{$i-2}) {$sum5+=$positive{"pre_len"}{$i-2}; $count++;}
        if (exists $positive{"pre_len"}{$i-1}) {$sum5+=$positive{"pre_len"}{$i-1}; $count++;}
        if (exists $positive{"pre_len"}{$i+1}) {$sum5+=$positive{"pre_len"}{$i+1}; $count++;}
        if (exists $positive{"pre_len"}{$i+2}) {$sum5+=$positive{"pre_len"}{$i+2}; $count++;}
        $positive{"pre_len"}{$i}=$sum5/5;
    }
}

#for (my $i=0; $i <= 100; $i++) { print $i,"\t",$positive{"stab"}{$i/100},"\n";}
## turn pdf to cdf
#foreach my $type (keys %positive) {
#    my $cdf=0;
#    #foreach my $id (sort{$a <=> $b} keys %{$positive{$type}}) {
#    #    $positive{$type}{$id}=$positive{$type}{$id}+$cdf;
#    #    $cdf=$positive{$type}{$id};
#    #}
#    foreach my $id (sort{$a <=> $b} keys %{$positive{$type}}) {
#        print join("\t","pos",$type,$id,$positive{$type}{$id}),"\n";
#    }
#}


$count=0;
while(<IN2>) {
    chomp;
    $count++;
    my @a=split("\t",$_);
    my $mfe=round(1-60*$a[30]/$a[1]);
    $background{"mfe"}{$mfe}++;
    my $stability=int($a[31]*100);
    if ($stability eq 0) {$stability = 1;}
    $background{"stab"}{$stability}++;
    my $p5=$a[32];
    $background{"p5"}{$p5}++;
    my $overhang=abs($a[34] - 2);
    $background{"overhang"}{$overhang}++;
    #my $internal_pair=$a[36];
    #$background{"int_pair"}{$internal_pair}++;
    #my $basepair=$a[39];
    #$background{"basepair"}{$basepair}++;
    if ($a[13] ne 0) {$background{"mat_len"}{$a[13]}++;}
    if ($a[26] ne 0) {$background{"mat_len"}{$a[26]}++;}
    my $pre_length=$a[1];
    $background{"pre_len"}{$pre_length}++;
}
# smooth hist to pdf
##### basepair
#$tmp=0;
#for (my $i=0; $i <= 25; $i++) {
#    if (exists $background{"basepair"}{$i}) { $background{"basepair"}{$i} = $background{"basepair"}{$i} / $count;}
#    else { $background{"basepair"}{$i}=0.5/$count; }
#    if($background{"basepair"}{$i} > $tmp) {$tmp=$background{"basepair"}{$i};}
#}
#for (my $i=1; $i < 25; $i++) {
#    if ($background{"basepair"}{$i} ne $tmp) {
#        my $count=1;
#        my $sum5=$background{"basepair"}{$i};
#        if (exists $background{"basepair"}{$i-2}) {$sum5+=$background{"basepair"}{$i-2}; $count++;}
#        if (exists $background{"basepair"}{$i-1}) {$sum5+=$background{"basepair"}{$i-1}; $count++;}
#        if (exists $background{"basepair"}{$i+1}) {$sum5+=$background{"basepair"}{$i+1}; $count++;}
#        if (exists $background{"basepair"}{$i+2}) {$sum5+=$background{"basepair"}{$i+2}; $count++;}
#        $background{"basepair"}{$i}=$sum5/5;
#    }
#}
##### int_pair
#$tmp=0;
#for (my $i=0; $i <= 20; $i++) {
#    if (exists $background{"int_pair"}{$i}) { $background{"int_pair"}{$i} = $background{"int_pair"}{$i} / $count;}
#    else { $background{"int_pair"}{$i}=0.5/$count; }
#    if($background{"int_pair"}{$i} > $tmp) {$tmp=$background{"int_pair"}{$i};}
#}
#for (my $i=1; $i < 20; $i++) {
#    if ($background{"int_pair"}{$i} ne $tmp) {
#        my $count=1;
#        my $sum5=$background{"int_pair"}{$i};
#        if (exists $background{"int_pair"}{$i-2}) {$sum5+=$background{"int_pair"}{$i-2}; $count++;}
#        if (exists $background{"int_pair"}{$i-1}) {$sum5+=$background{"int_pair"}{$i-1}; $count++;}
#        if (exists $background{"int_pair"}{$i+1}) {$sum5+=$background{"int_pair"}{$i+1}; $count++;}
#        if (exists $background{"int_pair"}{$i+2}) {$sum5+=$background{"int_pair"}{$i+2}; $count++;}
#        $background{"int_pair"}{$i}=$sum5/5;
#    }
#}
##### mfe
$tmp=0;
for (my $i=0; $i <= 70; $i++) {
    if (exists $background{"mfe"}{$i}) { $background{"mfe"}{$i} = $background{"mfe"}{$i} / $count;}
    else { $background{"mfe"}{$i}=0.5/$count; }
    if($background{"mfe"}{$i} > $tmp) {$tmp=$background{"mfe"}{$i};}
}
for (my $i=1; $i < 70; $i++) {
    if ($background{"mfe"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$background{"mfe"}{$i};
        if (exists $background{"mfe"}{$i-2}) {$sum5+=$background{"mfe"}{$i-2}; $count++;}
        if (exists $background{"mfe"}{$i-1}) {$sum5+=$background{"mfe"}{$i-1}; $count++;}
        if (exists $background{"mfe"}{$i+1}) {$sum5+=$background{"mfe"}{$i+1}; $count++;}
        if (exists $background{"mfe"}{$i+2}) {$sum5+=$background{"mfe"}{$i+2}; $count++;}
        $background{"mfe"}{$i}=$sum5/5;
    }
}
##### overhang
$tmp=0;
for (my $i=0; $i <= 50; $i++) {
    if (exists $background{"overhang"}{$i}) { $background{"overhang"}{$i} = $background{"overhang"}{$i} / $count;}
    else { $background{"overhang"}{$i}=0.5/$count; }
    if($background{"overhang"}{$i} > $tmp) {$tmp=$background{"overhang"}{$i};}
}
for (my $i=1; $i < 50; $i++) {
    if ($background{"overhang"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$background{"overhang"}{$i};
        if (exists $background{"overhang"}{$i-2}) {$sum5+=$background{"overhang"}{$i-2}; $count++;}
        if (exists $background{"overhang"}{$i-1}) {$sum5+=$background{"overhang"}{$i-1}; $count++;}
        if (exists $background{"overhang"}{$i+1}) {$sum5+=$background{"overhang"}{$i+1}; $count++;}
        if (exists $background{"overhang"}{$i+2}) {$sum5+=$background{"overhang"}{$i+2}; $count++;}
        $background{"overhang"}{$i}=$sum5/5;
    }
}
##### p5
$tmp=0;
for (my $i=0; $i <= 50; $i++) {
    if (exists $background{"p5"}{$i}) { $background{"p5"}{$i} = $background{"p5"}{$i} / $count;}
    else { $background{"p5"}{$i}=0.5/$count; }
    if($background{"p5"}{$i} > $tmp) {$tmp=$background{"p5"}{$i};}
}
for (my $i=1; $i < 50; $i++) {
    if ($background{"p5"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$background{"p5"}{$i};
        if (exists $background{"p5"}{$i-2}) {$sum5+=$background{"p5"}{$i-2}; $count++;}
        if (exists $background{"p5"}{$i-1}) {$sum5+=$background{"p5"}{$i-1}; $count++;}
        if (exists $background{"p5"}{$i+1}) {$sum5+=$background{"p5"}{$i+1}; $count++;}
        if (exists $background{"p5"}{$i+2}) {$sum5+=$background{"p5"}{$i+2}; $count++;}
        $background{"p5"}{$i}=$sum5/5;
    }
}
##### stab
$tmp=0;
for (my $i=0; $i <= 100; $i++) {
    if (exists $background{"stab"}{$i}) { $background{"stab"}{$i} = $background{"stab"}{$i} / $count;}
    else { $background{"stab"}{$i}=0.005/$count; }
    if($background{"stab"}{$i} > $tmp) {$tmp=$background{"stab"}{$i};}
}
for (my $i=1; $i < 100; $i++) {
    if ($background{"stab"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$background{"stab"}{$i};
        if (exists $background{"stab"}{$i-2}) {$sum5+=$background{"stab"}{$i-2}; $count++;}
        if (exists $background{"stab"}{$i-1}) {$sum5+=$background{"stab"}{$i-1}; $count++;}
        if (exists $background{"stab"}{$i+1}) {$sum5+=$background{"stab"}{$i+1}; $count++;}
        if (exists $background{"stab"}{$i+2}) {$sum5+=$background{"stab"}{$i+2}; $count++;}
        $background{"stab"}{$i}=$sum5/5;
    }
}
##### mat_length
$tmp=0;
for (my $i=15; $i <= 36; $i++) {
    if (exists $background{"mat_len"}{$i}) { $background{"mat_len"}{$i} = $background{"mat_len"}{$i} / $count;}
    else { $background{"mat_len"}{$i}=0.5/$count; }
    if($background{"mat_len"}{$i} > $tmp) {$tmp=$background{"mat_len"}{$i};}
}
for (my $i=16; $i < 36; $i++) {
    if ($background{"mat_len"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$background{"mat_len"}{$i};
        if (exists $background{"mat_len"}{$i-2}) {$sum5+=$background{"mat_len"}{$i-2}; $count++;}
        if (exists $background{"mat_len"}{$i-1}) {$sum5+=$background{"mat_len"}{$i-1}; $count++;}
        if (exists $background{"mat_len"}{$i+1}) {$sum5+=$background{"mat_len"}{$i+1}; $count++;}
        if (exists $background{"mat_len"}{$i+2}) {$sum5+=$background{"mat_len"}{$i+2}; $count++;}
        $background{"mat_len"}{$i}=$sum5/5;
    }
}
##### pre_len
$tmp=0;
for (my $i=40; $i <= 100; $i++) {
    if (exists $background{"pre_len"}{$i}) { $background{"pre_len"}{$i} = $background{"pre_len"}{$i} / $count;}
    else { $background{"pre_len"}{$i}=0.5/$count; }
    if($background{"pre_len"}{$i} > $tmp) {$tmp=$background{"pre_len"}{$i};}
}
for (my $i=41; $i < 100; $i++) {
    if ($background{"pre_len"}{$i} ne $tmp) {
        my $count=1;
        my $sum5=$background{"pre_len"}{$i};
        if (exists $background{"pre_len"}{$i-2}) {$sum5+=$background{"pre_len"}{$i-2}; $count++;}
        if (exists $background{"pre_len"}{$i-1}) {$sum5+=$background{"pre_len"}{$i-1}; $count++;}
        if (exists $background{"pre_len"}{$i+1}) {$sum5+=$background{"pre_len"}{$i+1}; $count++;}
        if (exists $background{"pre_len"}{$i+2}) {$sum5+=$background{"pre_len"}{$i+2}; $count++;}
        $background{"pre_len"}{$i}=$sum5/5;
    }
}


## turn pdf to cdf
#foreach my $type (keys %background) {
#    my $cdf=0;
#    foreach my $id (sort{$a <=> $b}  keys %{$background{$type}}) {
#        $background{$type}{$id}=$background{$type}{$id}/$count+$cdf;
#        $cdf=$background{$type}{$id};
#    }
#    foreach my $id (sort{$a <=> $b}  keys %{$background{$type}}) {
#        print join("\t","bgr",$type,$id,$background{$type}{$id}),"\n";
#    }
#}


# this prior can be changed
my $pre=0.5;
my $bgr=1-$pre;

$count=0;
while(<IN3>) {
    chomp;
    $count++;
    my @a=split("\t",$_);
    my $mfe=round(1-60*$a[30]/$a[1]);
    my $stability=int($a[31]*100);
    if ($stability eq 0) {$stability = 1;}
    my $p5=$a[32];
    my $overhang=abs($a[34] - 2);
    #my $internal_pair=$a[36];
    #my $basepair=$a[39];
    my $mature_length=$a[13];
    if ($a[29] > $a[16]) {$mature_length=$a[26]}
    my $pre_length=$a[1];
    my $pr1=$positive{"mfe"}{$mfe}*$pre/($positive{"mfe"}{$mfe}*$pre + $background{"mfe"}{$mfe}*$bgr);
    #my $pr2=$positive{"stab"}{$stability}*$pre;
    my $pr2=$positive{"stab"}{$stability}*$pre/($positive{"stab"}{$stability}*$pre + $background{"stab"}{$stability}*$bgr);
    my $pr3=$positive{"p5"}{$p5}*$pre/($positive{"p5"}{$p5}*$pre + $background{"p5"}{$p5}*$bgr);
    my $pr4=$positive{"overhang"}{$overhang}*$pre/($positive{"overhang"}{$overhang}*$pre + $background{"overhang"}{$overhang}*$bgr);
    #my $pr5=$positive{"int_pair"}{$internal_pair}*$pre/($positive{"int_pair"}{$internal_pair}*$pre + $background{"int_pair"}{$internal_pair}*$bgr);
    #my $pr6=$positive{"basepair"}{$basepair}*$pre/($positive{"basepair"}{$basepair}*$pre + $background{"basepair"}{$basepair}*$bgr);
    my $pr7=$positive{"mat_len"}{$mature_length}*$pre/($positive{"mat_len"}{$mature_length}*$pre + $background{"mat_len"}{$mature_length}*$bgr);
    my $pr8=$positive{"pre_len"}{$pre_length}*$pre/($positive{"pre_len"}{$pre_length}*$pre + $background{"pre_len"}{$pre_length}*$bgr);
    my $Pf=$pr1*$pr2*$pr3*$pr4*$pr7*$pr8/($pr1*$pr2*$pr3*$pr4*$pr7*$pr8 + (1-$pr1)*(1-$pr2)*(1-$pr3)*(1-$pr4)*(1-$pr7)*(1-$pr8));
    my $score1=$positive{"mfe"}{$mfe}*$positive{"stab"}{$stability}*$positive{"p5"}{$p5}*$positive{"overhang"}{$overhang}*$positive{"mat_len"}{$mature_length}*$positive{"pre_len"}{$pre_length}*$pre+0.0000001;
    my $score2=$background{"mfe"}{$mfe}*$background{"stab"}{$stability}*$background{"p5"}{$p5}*$background{"overhang"}{$overhang}*$background{"mat_len"}{$mature_length}*$background{"pre_len"}{$pre_length}*$bgr+0.0000001;
    my $score=log($score1/$score2);
    print OUT join("\t",$pr1,$pr2,$pr3,$pr4,$pr7,$pr8,$Pf,$score,@a),"\n";
}

sub max2 {
    my($a, $b) = @_;
    return ($a>$b ? $a : $b);
}

sub min2  {
    my($a, $b) = @_;
    return ($a<$b ? $a : $b);
}

sub round {
    my($number) = shift;
    return int($number + .5);
    
}
