#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"*_reversed\"  \"output basename\" " if (@ARGV ne 2);
my $filein=$ARGV[0];
my $fileout=$ARGV[1];
open IN,$filein;    # _reversed
open OUT,">".$fileout.".1";
open OUT2,">".$fileout.".2";
open OUT3,">".$fileout."_reject";
open OUT33,">".$fileout."_accept";
open OUT4,">".$fileout.".2.cf";
open OUT5,">".$fileout.".2.maturefam5p";
open OUT6,">".$fileout.".2.maturefam3p";
open OUT7,">".$fileout.".goodmature";
my $control=0;
my %maturefamily5p;
my %maturefamily3p;

sub round { int(0.5+shift)};

while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    my $prec_length=length($a[1]);
    my $qua='';
    if (scalar(@a) eq 5) {$qua=$a[2]}
    my @b;
    my $TCOUNT=0;
    if ($qua ne '') {
        $TCOUNT=$a[3];
        @b=split(/\#/,$a[4]);
    }
    else {
        $TCOUNT=$a[2];
        @b=split(/\#/,$a[3]);
    }
    #0# record the matures(in good length) and initialize blocks
    my @array;  # 5' start
    my @good_array;
    my $total_unique=0;
    my $good_unique=0;  # sequences of size [17,25]
    my $total_weight=0;
    my $good_weight=0;  # sequences of size [17,25]
    my $peak_length_min=17;
    my $peak_length_max=25;
    my %matures;    # position->length
    my %mature_name;
    for (my $i=0; $i < (36+$prec_length); $i++) {$array[$i]=0; }
    for (my $i=0; $i < $TCOUNT; $i++) {
        my @tmp=split(/\|/,$b[$i]);
        my @info=split(/\_\_/,$tmp[0]);
        my $infocnt=1;
        if (scalar(@info) > 0) { $infocnt=$info[-1]; }
        $total_unique++;
        $total_weight+=$infocnt;
        $array[$tmp[2]]+=$infocnt;
        $matures{$tmp[2]}{$tmp[1]}=$infocnt;
        $mature_name{$tmp[2]}{$tmp[1]}=$tmp[0];
        if (($tmp[1] >= $peak_length_min) and ($tmp[1] <= $peak_length_max)) {
            $good_unique++;
            $good_weight+=$infocnt;
            $good_array[$tmp[2]]+=$infocnt;
        }
    }
    print OUT7 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$good_unique,"\t",$total_weight,"\t",$good_weight,"\t",sprintf("%.3f",$good_unique/$total_unique),"\t",sprintf("%.3f",$good_weight/$total_weight),"\n";
    if (($good_weight/$total_weight < 0.9) or ($good_unique/$total_unique < 0.66)) {
        print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$good_unique,"\t",$total_weight,"\t",$good_weight,"\t",sprintf("%.3f",$good_unique/$total_unique),"\t",sprintf("%.3f",$good_weight/$total_weight),"\t","Rejection-0: mature length.\n";
        next;
    }
    #1# first call blocks using unique reads, then modify block borders resulting less-noised block
    #1# single read at both end of block are omitted, if expression is high enough
    #1# too long block (>= 20) means smear, Rejection 1: smear
    #1# the cut-off (now == 15) may change
    my %block;
    my $block_count=-1;
    my $block_start=0;
    my $block_end=0;
    my $block_readend=0;
    my $block_gap=-1;
    my $block_report='';
    my $BLOCK_conti=5;
    my $BLOCK_smear=15;
    my $smear=0;
    for (my $i=0; $i < (36+$prec_length); $i++) {
        if ($array[$i] ne 0) {
            if ($block_gap eq -1) {
                $block_count++;
                $block_start=$i;
                $block_end=$i;
                $block_gap=$i;
            }
            # end of the previous block, start a new one
            elsif (($i - $block_gap) > $BLOCK_conti) {
                $block_end=$block_gap;
                foreach my $s (sort{$b<=>$a} keys %{$matures{$block_end}}) {
                    $block_readend=$block_end + $s -1;
                    last;
                }
                if ($block_count >= 1) {
                    # if this block and previous block end very close < 10, merge this and previouse one
                    my @tmp=split("\t",$block{$block_count - 1});
                    if (($block_readend - $tmp[2]) < 10) {
                        # 5p-mature, merge block-1 into block-0
                        if (($block_count eq 1) and (($block_start/$prec_length) < 0.6)){
                            $tmp[2]=$block_readend > $tmp[2] ? $block_readend : $tmp[2];
                            $block{$block_count - 1}=join("\t",@tmp);
                            #if (($block_end - $block_start + 1) >= $BLOCK_smear) {$smear=1;}
                            #$block_count++;
                            $block_gap=$i;
                            $block_start=$i;
                            $block_end=$i;
                        }
                        # 3p-mature, merge block-k into block-(k+1)
                        else {
                            $tmp[0]=$tmp[0] > $block_start ? $block_start : $tmp[0];
                            $block{$block_count - 1}=join("\t",$tmp[0],$block_end,$block_readend,($block_end - $block_start + 1));
                            #if (($block_end - $block_start + 1) >= $BLOCK_smear) {$smear=1;}
                            #$block_count++;
                            $block_gap=$i;
                            $block_start=$i;
                            $block_end=$i;
                        }
                    }
                    else {
                        $block{$block_count}=join("\t",$block_start,$block_end,$block_readend,($block_end - $block_start + 1));
                        #if (($block_end - $block_start + 1) >= $BLOCK_smear) {$smear=1;}
                        $block_count++;
                        $block_gap=$i;
                        $block_start=$i;
                        $block_end=$i;
                    }
                }
                # first block, just record
                else {
                    $block{$block_count}=join("\t",$block_start,$block_end,$block_readend,($block_end - $block_start + 1));
                    #if (($block_end - $block_start + 1) >= $BLOCK_smear) {$smear=1;}
                    $block_count++;
                    $block_gap=$i;
                    $block_start=$i;
                    $block_end=$i;
                }
            }
            else {$block_gap=$i;$block_end=$i;}
        }
    }
    # save the last block info.
    foreach my $s (sort{$b<=>$a} keys %{$matures{$block_end}}) {
        $block_readend=$block_end + $s -1;
        last;
    }
    if ($block_count >= 1) {
        my @tmp=split("\t",$block{$block_count - 1});
        if (($block_readend - $tmp[2]) < 10) {
            if (($block_count eq 1) and (($block_start/$prec_length) < 0.6)) {
                # 5p-mature, merge block-1 into block-0
                $tmp[2]=$block_readend > $tmp[2] ? $block_readend : $tmp[2];
                $block{$block_count - 1}=join("\t",@tmp);
                $block_count--;
            }
            else {
                # 3p-mature, merge block-k into block-(k+1)
                $tmp[0]=$tmp[0] > $block_start ? $block_start : $tmp[0];
                $block{$block_count - 1}=join("\t",$tmp[0],$block_end,$block_readend,($block_end - $block_start + 1));
                $block_count--;
            }
        }
        else {
            $block{$block_count}=join("\t",$block_start,$block_end,$block_readend,($block_end - $block_start + 1));
            #else {$block_report=$block_report."block ".$block_count." : ".$block{$block_count}."\n";}
        }
    }
    else {
        $block{$block_count}=join("\t",$block_start,$block_end,$block_readend,($block_end - $block_start + 1));
        #if (($block_end - $block_start + 1) >= $BLOCK_smear) {$smear=1;}
    }
    
    #1.1# remove middle blocks with only one read, 5p-most and 3p-most block untouched
    for(my $i=0; $i <= $block_count; $i++) {
        my $Tsave='';
        my @TMP=split("\t",$block{$i});
        if (($TMP[0] > 5) and (($prec_length - $TMP[2]) > 7 ) and ($TMP[3] == 1)) {
            my $w=0;
            foreach my $mlength (keys %{$matures{$TMP[0]}}) {
                $w+=$matures{$TMP[0]}{$mlength};
            }
            if (($w < 5) or ($w/$total_weight < 0.005)) {$Tsave=1;}
        }
        if ($Tsave ne '') {
            for(my $j=$i; $j<$block_count; $j++) {
                $block{$j}=$block{$j+1};
            }
            $i=$i-1;
            $block_count=$block_count-1;
        }
    }
    
    #1.2# modify blocks borders if exp > 500
    if ($total_weight >= 500) {
        for(my $i=0; $i <= $block_count; $i++) {
            my @TMP=split("\t",$block{$i});
            my $f_left=0;
            for (my $i=$TMP[0]; $i <$TMP[1]; $i++) {
                if (exists $good_array[$i]) {
                    if ($good_array[$i] <= 2) {$f_left=$i;}
                    else {$f_left=$i; last;}
                }
            }
            #if (($f_left - 1) > $TMP[0]) {$TMP[0] = $f_left;}
            if ($f_left > $TMP[0]) {$TMP[0] = $f_left;}
            my $f_right=999;
            for (my $i=$TMP[1]; $i>$TMP[0]; $i--) {
                if (exists $good_array[$i]) {
                    if ($good_array[$i] <= 2) {$f_right=$i;}
                    else {$f_right=$i;last;}
                }
            }
            if ($f_right < $TMP[1]) {$TMP[1] = $f_right;}
            my $b_unique=0;
            my $b_weight=0;
            my $bg_unique=0;
            my $bg_weight=0;
            for (my $j=$TMP[0]; $j <= $TMP[1]; $j++) {
                $b_weight+=$array[$j];
                if (exists $good_array[$j]) {$bg_weight+=$good_array[$j];}
                foreach my $mlength (sort{$a <=> $b} keys %{$matures{$j}}) {
                    $b_unique++;
                    if (($mlength >= $peak_length_min) and ($mlength <= $peak_length_max)) {$bg_unique++;}
                }
            }
            $TMP[3]=$TMP[1] - $TMP[0] + 1;
            $block{$i}=join("\t",@TMP)."\t".$b_unique."\t".$bg_unique."\t".$b_weight."\t".$bg_weight;
        }
    }
    # modify block content as well
    else {
        for(my $i=0; $i <= $block_count; $i++) {
            my @TMP=split("\t",$block{$i});
            my $b_unique=0;
            my $b_weight=0;
            my $bg_unique=0;
            my $bg_weight=0;
            for (my $j=$TMP[0]; $j <= $TMP[1]; $j++) {
                $b_weight+=$array[$j];
                if (exists $good_array[$j]) {$bg_weight+=$good_array[$j];}
                foreach my $mlength (sort{$a <=> $b} keys %{$matures{$j}}) {
                    $b_unique++;
                    if (($mlength >= $peak_length_min) and ($mlength <= $peak_length_max)) {$bg_unique++;}
                }
            }
            $TMP[3]=$TMP[1] - $TMP[0] + 1;
            $block{$i}=join("\t",@TMP)."\t".$b_unique."\t".$bg_unique."\t".$b_weight."\t".$bg_weight;
        }
    }
    #### $block{$i} => start    end    readend    size    unique    good_unique    weight    good_weight
    
    #1.3# report block basic info and ID smear
    $block_report='';
    $block_report="block 0"." : ".$block{0}."\n";
    my @TMP=split("\t",$block{0});
    if($TMP[3] >= $BLOCK_smear) {$smear++;}
    for (my $i=1; $i<=$block_count; $i++) {
        $block_report=$block_report."block ".$i." : ".$block{$i}."\n";
        my @TMP=split("\t",$block{$i});
        if($TMP[3] >= $BLOCK_smear) {$smear++;}
    }
    
    #1.4# remove smear
    if ($smear > 0) {
        print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-1: smear.\n",$block_report,"\n";
        next;
    }
    else {
        if ($block_count > 2) {
            print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-10: too many blocks(chimeric read).\n",$block_report,"\n";
            next;
        }
        else {
            my $smear=0;
            if ($block_count eq 1) {
                my @a0=split("\t",$block{0});
                my @a1=split("\t",$block{1});
                if ((abs($a1[0] - $a0[2]) < 3) and ($a1[3] > 5) and ($a0[3] > 5)){
                    my $uniqread=1;
                    my $endpos=0;
                    for(my $i=$a0[0]; $i<=$a0[1]; $i++) {
                        if (exists $matures{$i}) {
                            foreach my $len (keys %{$matures{$i}}) {
                                $uniqread++;
                                $endpos+=$len;
                            }
                        }
                    }
                    $endpos=round($endpos/$uniqread);
                    if (($a1[0] - $endpos) eq 1) {
                        print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-11: possible chimeric read.\n",$block_report,"\n";
                        next;
                    }
                }
            }
            elsif ($block_count eq 2) {
                my @a0=split("\t",$block{0});
                my @a1=split("\t",$block{1});
                my @a2=split("\t",$block{2});
                if ((abs($a1[0] - $a0[2]) < 3) and ($a1[3] > 5) and ($a0[3] > 5)){
                    my $uniqread=1;
                    my $endpos=0;
                    for(my $i=$a0[0]; $i<=$a0[1]; $i++) {
                        if (exists $matures{$i}) {
                            foreach my $len (keys %{$matures{$i}}) {
                                $uniqread++;
                                $endpos+=$len;
                            }
                        }
                    }
                    $endpos=round($endpos/$uniqread);
                    if (($a1[0] - $endpos) eq 1) {
                        print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-12: possible chimeric read.\n",$block_report,"\n";
                        next;
                    }
                }
                if ((abs($a2[0] - $a1[2]) < 3) and ($a2[3] > 5) and ($a1[3] > 5)){
                    my $uniqread=1;
                    my $endpos=0;
                    for(my $i=$a1[0]; $i<=$a1[1]; $i++) {
                        if (exists $matures{$i}) {
                            foreach my $len (keys %{$matures{$i}}) {
                                $uniqread++;
                                $endpos+=$len;
                            }
                        }
                    }
                    $endpos=round($endpos/$uniqread);
                    if (($a2[0] - $endpos) eq 1) {
                        print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-13: possible chimeric read.\n",$block_report,"\n";
                        next;
                    }
                }
            }
            print OUT33 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\n",$block_report,"\n";
        }
    }
    #else {print $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\n",$block_report,"\n";}
    
    #2# call peak inside blocks
    my $peak_start=-1;
    my $peak_length=-1;
    my $peak_height=-1;
    my $block_unique=0; 
    my $block_weighted=0;
    my $peak_window=0;
    my $peak_window_size=1; # meaning window of [-1,0,+1], size 3
    my $tobe_mature=-1;
    for(my $i=0; $i<=$block_count; $i++) {
        my @c=split("\t",$block{$i});
        for(my $j=$c[0]; $j<=$c[1]; $j++) {
            if (exists $matures{$j}) {
                foreach my $mlength (sort{$a <=> $b} keys %{$matures{$j}}) {
                    #$block_unique++;
                    #$block_weighted+=$matures{$j}{$mlength};
                    if (($mlength <= $peak_length_max) and ($matures{$j}{$mlength} > $peak_height)) {
                    #if ($matures{$j}{$mlength} > $peak_height) {
                        $peak_height=$matures{$j}{$mlength};
                        $peak_start=$j;
                        $peak_length=$mlength;
                    }
                }
            }
        }
        if ($peak_height eq -1) {
            $block{$i}='';
            $peak_height=-1;
            $peak_window=0;
            $block_unique=0;
            $block_weighted=0;
            next;
        }
        for(my $j=($peak_start-$peak_window_size); $j<=($peak_start+$peak_window_size); $j++) {
            if (exists $matures{$j}) {
                foreach my $mlength (sort{$a <=> $b} keys %{$matures{$j}}) {
                    if (($mlength >= $peak_length_min) and ($mlength <= $peak_length_max)) {
                        $peak_window+=$matures{$j}{$mlength};
                    }
                }
            }
        }
        $block{$i}=join("\t",$block{$i},$peak_start,$peak_length,$peak_height,$peak_window);
        #$block{$i}=join("\t",$block{$i},$peak_start,$peak_length,$peak_height,$peak_window,$block_unique,$block_weighted);
        $peak_height=-1;
        $peak_window=0;
        $block_unique=0;
        $block_weighted=0;
    }
    my $tmp_count=0;
    for(my $i=0; $i<=$block_count; $i++) {
        if ($block{$i} ne ''){
            $block{$tmp_count}=$block{$i};
            $tmp_count++;
        }
    }
    $block_count=$tmp_count - 1;
    if ($block_count < 0) {
        print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-20: too long matures.\n",$block_report,"\n";
        next;
    }
    #### $block{$i} => start    end    readend    size    unique    good_unique    weight    good_weight    peak_start    peak_length    peak_height    peak_window
    ####                0       1       2           3       4           5           6           7               8           9               10              11
    
    #3# 5p-most peak start_pos in [0,5] AND 3p-most peak end_pos in [-8,-1]
    #3# otherwise Rejection 2: inconsistent with Dicer processing
    my $leftmost=-1;
    my $leftweight=0;   # unique short reads species
    my @c=split("\t",$block{0});
    if ($c[0] <= 5) {
        if ($c[8] <= 3) {$leftmost=0; $leftweight=$c[4];}   # 5p mature found
        else {$leftmost=-100;}  # wrong 5p mature, should discard
    }
    elsif ($c[0] < 15) {$leftmost=-100;}    # wrong 5p mature, should discard
    else {$leftmost=-10;}   # no 5p mature
    
    my $rightmost=-1;
    my $rightweight=0;
    @c=split("\t",$block{$block_count});
    #print $a[0],"\t",$block_count,"\t",$block{$block_count},"\n";
    if (($prec_length - $c[2]) < 8 ) {
        if (($prec_length - $c[8] - $c[9]) < 8) {$rightmost=$block_count; $rightweight=$c[4];}  # 3p mature found
        else {$rightmost=-100;} # wrong 3p mature, should discard
    }
    else {$rightmost=-10;}  # no 3p mature
    if (($leftmost eq -100) or ($rightmost eq -100)) {
        print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-2: inconsistent with Dicer processing.\n";
        for(my $i=0; $i<=$block_count; $i++) { print OUT3 $block{$i},"\n"; }
        print OUT3 "\n";
        next;
    }
    #else {
    #    print $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\n5pmost :",$leftmost,"\n3pmost : ",$rightmost,"\n";;
    #    for(my $i=0; $i<=$block_count; $i++) { print $block{$i},"\n"; }
    #}
    
    #### $block{$i} => start    end    readend    size    unique    good_unique    weight    good_weight    peak_start    peak_length    peak_height    peak_window
    ####                0       1       2           3       4           5           6           7               8           9               10              11
    #3# middle block should be too close to terminal blocks:
    #3# for 5p, block-1_start > (block-5p_start + peak1_length); for 3p, (block-n_start + peakn_length) < block-3p
    #3# otherwise Rejection 3: inconsistent with Dicer processing and block-smear
    if ($block_count > 0) {
        # see if 5p mature exists
        @c=split("\t",$block{0});
        if ($c[0] <= 5) {
            my @d=split("\t",$block{1});
            if (($d[0] < ($c[0] + $c[9] - 6)) and ($d[7] > 2)) {
                print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-3: inconsistent with Dicer processing and block_smear.\n";
                for(my $i=0; $i<=$block_count; $i++) { print OUT3 $block{$i},"\n"; }
                print OUT3 "\n";
                next;
            }
        }
        else {
            @c=split("\t",$block{$block_count});
            if (($prec_length - $c[2]) < 5) {
                my @d=split("\t",$block{$block_count - 1});
                if (($c[0] < ($d[8] + $d[9] - 6)) and ($d[7] > 2)) {
                    print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-3: inconsistent with Dicer processing and block_smear.\n";
                    for(my $i=0; $i<=$block_count; $i++) { print OUT3 $block{$i},"\n"; }
                    print OUT3 "\n";
                    next;
                }
            }
        }
    }
    
    #### $block{$i} => start    end    readend    size    unique    good_unique    weight    good_weight    peak_start    peak_length    peak_height    peak_window
    ####                0       1       2           3       4           5           6           7               8           9               10              11
    #4# remove middle blocks with only one read, 5p-most and 3p-most block untouched
    #4# If middle peak higher than the lower of 5p-most or 3p-most, Rejection 3: inconsistent with miRNA biogenesis
    my %save;
    if ($block_count > 0) {
        if (($leftweight > 0) and ($rightweight > 0)) {
            my $less_uniq=$leftweight < $rightweight ? $leftweight : $rightweight;
            for(my $i=0; $i<=$block_count; $i++) {
                my @c=split("\t",$block{$i});
                if (($leftmost ne $i) and ($rightmost ne $i)) {
                    if (($c[4] > 2*$less_uniq) or ($c[4] >= 12)) {
                        print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-4: more loop read species than that of matures.\n";
                        for(my $j=0; $j<=$block_count; $j++) { print OUT3 $block{$j},"\n"; }
                        print OUT3 "\n";
                        $less_uniq=-416;
                        next;
                    }
                    if ($c[6] > 1) {$save{$i}=$c[6];}   # $c[8] can be substituted by $c[6] (peak_window_weight instead of block_weight)
                }
            }
            if ($less_uniq eq -416) {next;}
        }
        elsif (($leftweight > 0) and ($rightweight eq 0)) {
            my $less_uniq=$leftweight;
            for(my $i=1; $i<=$block_count; $i++) {
                my @c=split("\t",$block{$i});
                if (($leftmost ne $i) and ($rightmost ne $i)) {
                    if (($c[4] > 2*$less_uniq) or ($c[4] >= 5)) {
                        print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-4: more loop read species than that of 5p-matures.\n";
                        for(my $j=0; $j<=$block_count; $j++) { print OUT3 $block{$j},"\n"; }
                        print OUT3 "\n";
                        $less_uniq=-416;
                        next;
                    }
                    if ($c[6] > 1) {$save{$i}=$c[6];}   # $c[8] can be substituted by $c[6] (peak_window_weight instead of block_weight)
                }
            }
            if ($less_uniq eq -416) {next;}
        }
        elsif (($leftweight eq 0) and ($rightweight > 0)) {
            my $less_uniq=$rightweight;
            for(my $i=0; $i<$block_count; $i++) {
                my @c=split("\t",$block{$i});
                if (($leftmost ne $i) and ($rightmost ne $i)) {
                    if (($c[4] > 2*$less_uniq) or ($c[4] >= 5)) {
                        print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-4: more loop read species than that of 3p-matures.\n";
                        for(my $j=0; $j<=$block_count; $j++) { print OUT3 $block{$j},"\n"; }
                        print OUT3 "\n";
                        $less_uniq=-416;
                        next;
                    }
                    if ($c[6] > 1) {$save{$i}=$c[6];}   # $c[8] can be substituted by $c[6] (peak_window_weight instead of block_weight)
                }
            }
            if ($less_uniq eq -416) {next;}
        }
        else {
            print OUT3 $a[0],"\t",$prec_length,"\t",$total_unique,"\t",$total_weight,"\t","Rejection-5: no mature reads.\n";
            for(my $i=0; $i<=$block_count; $i++) { print OUT3 $block{$i},"\n"; }
            print OUT3 "\n";
            next;
        }
    }
    
    #### $block{$i} => start    end    readend    size    unique    good_unique    weight    good_weight    peak_start    peak_length    peak_height    peak_window
    ####                0       1       2           3       4           5           6           7               8           9               10              11
    #5# output:         id    seq   qua   length       #short_uniq   #short_weight
    my $leftb=$leftmost;
    my $rightb=$rightmost;
    if ($leftb eq $rightb) {$rightb=-200;}
    
    my $leftinfo;
    if ($leftb >= 0) {
        my @c=split("\t",$block{$leftb});
        $leftinfo=join("\t",$block{$leftb},sprintf("%.3f",$c[11]/$total_weight));
    }
    else {
        $leftinfo=join("\t",0,0,0,0,0,0,0,0,0,0,0,0,0);
    }
    my $rightinfo;
    if ($rightb >= 0) {
        my @c=split("\t",$block{$rightb});
        $rightinfo=join("\t",$block{$rightb},sprintf("%.3f",$c[11]/$total_weight));
    }
    else {
        $rightinfo=join("\t",0,0,0,0,0,0,0,0,0,0,0,0,0);
    }
    #5# output:         id    seq    length       #short_uniq   #short_weight
    my $outfmt1;
    if ($qua ne '') {$outfmt1=join("\t",$a[0],$a[1],$qua,$prec_length,$total_unique,$total_weight,$leftinfo,$rightinfo);}
    else {$outfmt1=join("\t",$a[0],$a[1],$prec_length,$total_unique,$total_weight,$leftinfo,$rightinfo);}
    my $outfmt2=join("\t",$a[0],$prec_length,$total_unique,$total_weight,$leftinfo,$rightinfo);
    print OUT $outfmt1,"\n";
    print OUT2 $outfmt2,"\n";
    
    #### id   prec_length    uique    weight   12_left    12 right
    #### $block{$i} => start    end    readend    size    unique    good_unique    weight    good_weight    peak_start    peak_length    peak_height    peak_window  precentage
    ####                4        5        6         7       8          9           10          11              12           13              14              15        16
    ####                17       18       19       20       21         22          23          23              25           26              27              28        29
    #6# check sharpness, prepare for mature_calling
    my $prec_length_min=50;
    my $prec_length_max=94;
    my $block_size=15;
    $peak_length=26;
    my $f=0;
    @a=split("\t",$outfmt2);
    #                                                                 peak_height   block_weight           
    if (($a[1] > $prec_length_min) and ($a[1] < $prec_length_max) and (($a[14]+1)/($a[10]+1) >= 0.2) and (($a[27]+1)/($a[23]+1) >= 0.2)) {
        # 3p-most peak only
        if ($a[14] eq 0) {
            #if (($a[18] < $block_size) and  ($a[20] < $peak_length )) {
            if (($a[20] <= $block_size) and  ($a[26] >= $peak_length_min ) and  ($a[26] <= $peak_length_max )) {
                if ($a[3] >= 5000) {if ($a[29] >= 0.7) {$f=1;}}
                elsif ($a[3] <= 1000) {if ($a[29] >= 0.5) {$f=1;}}
                else {if ($a[29] >= 0.6) {$f=1;}}
            }
        }
        # 5p-most peak only
        elsif ($a[27] eq 0) {
            #if (($a[7] < $block_size)  and ($a[9] < $peak_length )) {
            if (($a[7] <= $block_size)  and ($a[13] >= $peak_length_min ) and ($a[13] <= $peak_length_max )) {
                if ($a[3] >= 5000) {if ($a[16] >= 0.7) {$f=1;}}
                elsif ($a[3] <= 1000) {if ($a[16] >= 0.5) {$f=1;}}
                else {if ($a[16] >= 0.6) {$f=1;}}
            }
        }
        else {
            # both peaks
            #           block_size                                           peak_length      
            #if (($a[7] < $block_size) and ($a[18] < $block_size) and ($a[20] < $peak_length ) and ($a[9] < $peak_length )) {
            if (($a[7] <= $block_size) and ($a[20] <= $block_size) and ((($a[13] >= $peak_length_min ) and ($a[13] <= $peak_length_max )) or (($a[26] >= $peak_length_min ) and  ($a[26] <= $peak_length_max )))) {
                if ($a[3] >= 5000) {if (($a[16] >= 0.7) or ($a[29] >= 0.7) or (($a[16]+$a[29]) >= 0.8)) {$f=1;}}
                elsif ($a[3] <= 1000) {if (($a[16] >= 0.5) or ($a[29] >= 0.5) or (($a[16]+$a[29]) >= 0.6)) {$f=1;}}
                else {if (($a[16] >= 0.5) or ($a[29] >= 0.5) or (($a[16]+$a[29]) >= 0.7)) {$f=1;}}
            }
        }
    }
    if($f eq 1) {
        #print OUT4 $outfmt2,"\n";
        #7# call mature-family
        my $mature_5p="NA";
        my $mature_3p="NA";
        if ($leftweight > 0) {
            @c=split("\t",$block{0});
            $mature_5p=$mature_name{$c[8]}{$c[9]};
            my $tmp="";
            for(my $i=$c[0]; $i <= $c[1]; $i++) {
                if (exists $matures{$i}) {
                    foreach my $j (keys %{$matures{$i}}) {
                        if ($tmp eq '') {$tmp=$mature_name{$i}{$j}}
                        else {$tmp=$tmp."\t".$mature_name{$i}{$j};}
                    }
                }
            }
            print OUT5 $mature_5p,"\t",$tmp,"\n";
        }
        if ($rightweight > 0) {
            @c=split("\t",$block{$block_count});
            $mature_3p=$mature_name{$c[8]}{$c[9]};
            my $tmp="";
            for(my $i=$c[0]; $i <= $c[1]; $i++) {
                if (exists $matures{$i}) {
                    foreach my $j (keys %{$matures{$i}}) {
                        if ($tmp eq '') {$tmp=$mature_name{$i}{$j}}
                        else {$tmp=$tmp."\t".$mature_name{$i}{$j};}
                    }
                }
            }
            print OUT6 $mature_3p,"\t",$tmp,"\n";
        }
        print OUT4 $outfmt2,"\t",$mature_5p,"\t",$mature_3p,"\n";
    }
    else {
        print OUT3 $outfmt2,"\t","Rejection-last: bad support signature.\n";
    }
    #$control++;
    #if ($control > 10) {last;}
}

