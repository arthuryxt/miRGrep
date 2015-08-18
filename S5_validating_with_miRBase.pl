#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"_reversed file\"   \"scoring\"  \"selected score\"  \"short_map_to_miRBase\"  \"output\"  \(optional main_ref_species default=\"mmu\"\)" if (@ARGV <= 5);
my $reversed=$ARGV[0];
my $scoring=$ARGV[1];
my $select=$ARGV[2];
my $maturemap=$ARGV[3];
my $output=$ARGV[4];
my $species="mmu";
if (scalar(@ARGV) > 5) { $species=$ARGV[5]; }
open IN1,$reversed;
open IN2,$scoring;
open IN3,$select;
open IN4,$maturemap;
open OUT,">".$output;
open OUT2,">".$output."_error";

sub fill($$$){   
    my $seq='';
    for(my $i=1; $i<$_[1]; $i++) {$seq=$seq." ";}
    my $tmp=substr($_[2],($_[1]-1),$_[0]);
    $seq=$seq.$tmp;
    #for(my $i=0; $i<$_[0]; $i++) {$seq=$seq."*";}
    for(my $i=($_[0]+$_[1]); $i<=100; $i++) {$seq=$seq." ";}
    my @a=split("",$tmp);
    my $left=0;
    my $right=0;
    for(@a){
        if($_ eq "("){$left++;}
        elsif($_ eq ")"){$right++;}
    }
    my $max=$left > $right ? $left : $right;
    return join("\t",$left,$right,$max,$seq);
}

my %need;
while(<IN3>) {
    chomp;
    my @a=split("\t",$_);
    $need{$a[0]}=$_;
}
close IN3;

print "read need done.\n";

my %Mature;
while(<IN4>) {
    chomp;
    my @a=split("\t",$_);
    $Mature{$a[0]}=$a[7];
}
close IN4;

print "read mapped to miRBase mature done.\n";

my %uniq;
while(<IN1>) {
    chomp;
    my @a=split("\t",$_);
    if (exists $need{$a[0]}) {
        $uniq{$a[0]}=$a[4];
    }
}
close IN1;

print "read reversed done.\n";

my %score;
while(<IN2>) {
    chomp;
    my @a=split("\t",$_);
    my $seq=<IN2>;
    my $str=<IN2>;
    if (exists $need{$a[0]}) {
        $score{$a[0]}=$seq.$str;
    }
}
close IN2;

print "read scoring done.\n";

sub count_left_bracket($) {
    my @a=split('',$_[0]);
    my $count=0;
    for(my $i=0; $i<scalar(@a); $i++) {
        if ($a[$i] eq "(") {$count++;}
        elsif ($a[$i] eq ")") {$count--;} 
    }
    return $count;
}

sub count_right_bracket($) {
    my @a=split('',$_[0]);
    my $count=0;
    for(my $i=0; $i<scalar(@a); $i++) {
        if ($a[$i] eq ")") {$count++;}
        elsif ($a[$i] eq "(") {$count--;} 
    }
    return $count;
}

sub get_left_bracket($$$$) {
    my $start=$_[0];        # 3' mature start
    my $length=$_[1];       # 3' mature length
    my $count=$_[2] - count_right_bracket(substr($_[3],($start+$length-2),2));        # 3' base-pair to be counted
    my @a=split('',$_[3]);
    my $leng=scalar(@a) - 1;
    my $pos=0;
    my $end=0;
    my $tmp_count=0;
    for(my $i=$leng; $i > ($start + $length -1 -2); $i--) {
        if ($a[$i] eq ")") {$tmp_count++;}
        elsif ($a[$i] eq "(") {$tmp_count--;}
    }
    for(my $i=0; ($i < $leng) and ($tmp_count > 0); $i++) {
        if ($a[$i] eq "(") {$tmp_count--;}
        elsif ($a[$i] eq ")") {$tmp_count++;}
        $pos=$i;
    }
    for(my $i=$pos; ($i<$leng) and ($count > 0); $i++) {
        if ($a[$i] eq "(") {$count--;}
        elsif ($a[$i] eq ")") {$count++;}
        $end=$i;
    }
    #$end=$end+2;
    return $pos."\t".$end;
}

sub get_right_bracket($$$$) {
    my $start=$_[0];    # 5' mature start
    my $length=$_[1];   # 5' mature length
    my $count=$_[2] - count_left_bracket(substr($_[3],($start+$length-2),2));   # base-pairs to be counted
    my @a=split('',$_[3]);
    my $leng=scalar(@a) - 1;
    my $end=$leng;  # should be the first right-bracket that coupled with the first left-bracket in the 5' mature part, and plus up to 2 if not exceed $leng;
    my $tmp_count=0;
    for(my $i=0; $i<$start; $i++) {
        if ($a[$i] eq "(") {$tmp_count++;}
    }
    for(my $i=$leng; ($i>0) and ($tmp_count > 0); $i--) {
        #if ($tmp_count eq 0) {last;}
        if ($a[$i] ne ")") {$end--;}
        elsif ($a[$i] eq ")") {$end--; $tmp_count--;}        
    }
    #$tmp_count=0;
    #for(my $i=$start; $i < $leng; $i++) {
    #    if ($a[$i] ne "(") {$end--; $tmp_count++;}
    #    else {last;}
    #}
    my $pos=0;
    for(my $i=$end; ($i>=0) and ($count >= 0); $i--) {
        if ($a[$i] eq "(") {$count++;}
        elsif ($a[$i] eq ")") {$count--;}
        $pos=$i;
    }
    if ($leng > $end) {$end++;}
    if ($leng > $end) {$end++;}
    return $pos."\t".$end;
}

foreach my $id (keys %need) {
    my $error_info='';
    my ($name)=split(/\//,$id);
    #open OUT, ">".$outputdir."\/".$name.".txt";
    my @a=split(/\#/,$uniq{$id});
    my @b=split("\n",$score{$id});
    my @b2=split("\t",$need{$id});
    my %rank;
    for(@a) {
        my @b=split(/\|/,$_);
        $rank{$b[2]}{$b[1]}=$b[0];
    }
    my $weight=0;
    my $count=0;
    my $buffer=$score{$id}."\n";
    my $mature=0;
    my $mature_pair=0;
    foreach my $pos(sort{$a <=> $b} keys %rank) {
        foreach my $len(sort{$a <=> $b} keys %{$rank{$pos}}) {
            my $result=fill($len,$pos,$b[1]);
            my ($left,$right,$mpair,$pat)=split("\t",$result);
            #print join("\t",$pat,$left,$right),"\n";
            my @c=split(/\_\_/,$rank{$pos}{$len});
            if ($c[1] > $mature) {$mature=$c[1]; $mature_pair=sprintf("%.3f",$mpair/$len);}
            if (exists $Mature{$rank{$pos}{$len}}) {
                $buffer=$buffer.$pat.$pos."\t".$len."\t".$c[1]."\t".$rank{$pos}{$len}."\t".$Mature{$rank{$pos}{$len}}."\n";
            }
            else {
                $buffer=$buffer.$pat.$pos."\t".$len."\t".$c[1]."\t".$rank{$pos}{$len}."\n";
            }
            #print OUT $pat,$pos,"\t",$len,"\t",$c[1],"\t",$rank{$pos}{$len},"\n";
            if ($left*$right > 0) {$count++; $weight+=$c[1];}
        }
    }
    #if ($mature_pair < 0.6) {
    #    $error_info="less than 60% nucleotides of mature part of candidate precursor are base-paired";
    #    if ($need{$id}=~m/mmu/) {
    #        my $save=$b2[-1];
    #        $b2[-1]=$count;
    #        print OUT2 join("\t",@b2),"\t",$weight,"\t",sprintf("%.6f",$count/$b2[2]),"\t",sprintf("%.6f",$weight/$b2[3]),"\t",$mature_pair,"\t",$save,"\n";
    #    }
    #    else {
    #        print OUT2 $need{$id},"\t",$count,"\t",$weight,"\t",sprintf("%.6f",$count/$b2[2]),"\t",sprintf("%.6f",$weight/$b2[3]),"\t",$mature_pair,"\n";
    #    }
    #    print OUT2 $buffer;
    #    next;
    #}
    
    # check if matures are properly paired
    #### id   prec_length    uique    weight   12_left    12 right
    #### $block{$i} => start    end    readend    size    unique    good_unique    weight    good_weight    peak_start    peak_length    peak_height    peak_window  precentage
    ####                4        5        6         7       8          9           10          11              12           13              14              15        16
    ####                17       18       19       20       21         22          23          23              25           26              27              28        29
    my $proper_check=1;
    my $pattern=$b[1];
    chomp $pattern;
    my $lenp=length($pattern);
    my @mat5=qw(-100 -100 0);
    my @mat3=qw(-100 -100 0);
    my $corespond_pos=-1;
    my $corespond_end=-1;
    my $mature_bp=0;
    if ($b2[13] ne 0) {
        # 5p-mature exists
        $mat5[0]=$b2[12]-1;  # start 0-based
        $mat5[1]=$b2[13];    # length
        $mat5[2]=$b2[14];    # height
    }
    if ($b2[26] ne 0) {
        # 5p-mature exists
        $mat3[0]=$b2[25]-1;  # start 0-based
        $mat3[1]=$b2[26];    # length
        $mat3[2]=$b2[27];    # height
    }
    $buffer=$buffer."mat5 = ".join("\t",@mat5)."\n";
    $buffer=$buffer."mat3 = ".join("\t",@mat3)."\n";
    # if both mature exists, length difference <= 4, otherwise discard
    if (($mat5[1] > 0) and ($mat3[1] > 0)) {
        if (abs($mat5[1] - $mat3[1]) > 4) {$proper_check=0; $error_info="both mature exists, length distance > 4";}
    }
    # the longer mature sequence length must >= 18, otherwise the mature is too short and should be discarded.
    my $longer_mature_length=$mat5[1] > $mat3[1] ? $mat5[1] : $mat3[1];
    if ($longer_mature_length < 18) {$proper_check=0; $error_info="longer mature length < 18";}
    $buffer=$buffer."proper_check1 = ".$proper_check."\n";
    $buffer=$buffer."error_info1 = ".$error_info."\n";
    
    if ($mat5[2] >= $mat3[2]) {
        # try to find the corresponding 3p-mature according to base-pair info, corrected for 2nt 3'p-overhang
        $mature_bp=count_left_bracket(substr($pattern,$mat5[0],$mat5[1])); # should > 0, as we have annotated 5p-mature
        ($corespond_pos,$corespond_end)=split("\t",get_right_bracket($mat5[0],$mat5[1],$mature_bp,$pattern));
        if ((($corespond_end - $corespond_pos + 1) <= 27) and (($corespond_end - $corespond_pos + 1) > 16)) {
            if ($mat3[1] > 0) {
                # exists 3p-mature, just check if overlap +/- 3nt
                if (abs($corespond_pos - $mat3[0]) <= 3) {$proper_check=1;}
                else {$proper_check=0; $error_info='existed 3p-mature differ with excised mature at 5-prime position > 3';}
            }
            if($proper_check eq 0) {
                # exists 3p-mature, check lenght
                if ($mat3[1] > 0) {
                    if (abs($corespond_end - $corespond_pos + 1 - $mat3[1]) < 5) {
                    #if ((abs($corespond_end - $corespond_pos + 1 - $mat5[1]) < 5) or (abs($corespond_end - $corespond_pos + 1 - $mat3[1]) <= 5)) {
                        $proper_check=1;
                        $error_info='';
                    }
                    else {$proper_check=0; $error_info='existed 3p-mature differ with excised mature on length > 5';}
                }
                else {
                    if (abs($corespond_end - $corespond_pos + 1 - $mat5[1]) < 5) {
                        $proper_check=1;
                        $error_info='';
                    }
                    else {$proper_check=0; $error_info='existed 3p-mature differ with excised mature on length > 5';}
                }
                
            }
        }
        else {$proper_check=0; $error_info='excised 3p-mature length wrong';}
    }
    else {
        # find the corresponding 5p-mature according to base-pair info
        $mature_bp=count_right_bracket(substr($pattern,$mat3[0],$mat3[1])); # should > 0, as we have annotated 3p-mature
        ($corespond_pos,$corespond_end)=split("\t",get_left_bracket($mat3[0],$mat3[1],$mature_bp,$pattern));
        if ((($corespond_end -$corespond_pos + 1) <= 27) and (($corespond_end - $corespond_pos + 1) > 16)) {
            if ($mat5[1] > 0) {
                # exists 5p-mature, just check if overlap +/- 3nt
                if (abs($corespond_pos - $mat5[0]) <= 3) {$proper_check=1;}
                else {$proper_check=0; $error_info='existed 5p-mature differ with excised mature at 5-prime position > 3';}
            }
            if($proper_check eq 0) {
                if ($mat5[1] > 0){
                    if(abs($mat5[1] - ($corespond_end - $corespond_pos + 1)) < 5) {
                    #if((abs($mat3[1] - ($corespond_end - $corespond_pos + 1)) < 5) or (abs($mat5[1] - ($corespond_end - $corespond_pos + 1)) < 5)) {
                        $proper_check=1;
                        $error_info='';
                    }
                    else {$proper_check=0; $error_info='existed 5p-mature differ with excised mature on length > 5';}
                }
                else {
                    if(abs($mat3[1] - ($corespond_end - $corespond_pos + 1)) < 5 ) {
                        $proper_check=1;
                        $error_info='';
                    }
                    else {$proper_check=0; $error_info='existed 5p-mature differ with excised mature on length > 5';}
                }
                
            }
        }
        else {$proper_check=0; $error_info='excised 5p-mature length wrong';}
    }
    $buffer=$buffer.$corespond_pos."\t".$corespond_end."\n";
    $buffer=$buffer."proper_check2 = ".$proper_check."\n";
    $buffer=$buffer."error_info2 = ".$error_info."\n";
    
    if ($error_info ne '') {$proper_check=0;}
    if ($proper_check eq 1) {
        # further check the position of "loop-sequence"
        if ($mat5[2] eq 0) {
            # check if loop-sequence start before 10
            foreach my $pos (sort{$a <=> $b} keys %rank) {
                if ($pos <= 10) {
                    $proper_check=0;
                    $error_info="no 5p mature and loop-sequence start before 10";
                    last;
                }
            }
        }
        elsif ($mat3[2] eq 0) {
            # check if loop-sequence end after $corespond_pos+5
            my $tmp_w_end=0;
            my $tmp_w_count=0;
            foreach my $pos (sort{$a <=> $b} keys %rank) {
                if ($pos > ($mat5[0] + $mat5[1] + 1)) {
                    foreach my $len (sort{$a <=> $b} keys %{$rank{$pos}}) {
                        my @c=split(/\_\_/,$rank{$pos}{$len});
                        $tmp_w_count+=$c[-1];
                        $tmp_w_end=$tmp_w_end + ($pos+$len-1)*$c[-1];
                    }
                }
            }
            if ($tmp_w_end > 0) {$tmp_w_end=$tmp_w_end/$tmp_w_count;}
            if (($tmp_w_end - $corespond_pos) > 3){$proper_check=0; $error_info="no 3p mature and loop-sequence start after 3nt of corresponding 3p-start";}
        }
    }
    if ($error_info eq '') {$proper_check = 1}
    if ($proper_check eq 1) {
        foreach my $pos(sort{$a <=> $b} keys %rank) {
            if (($pos > (5+$b2[5])) and ($pos < ($b2[17]-5))) {
                foreach my $len(sort{$a <=> $b} keys %{$rank{$pos}}) {
                    if (($len > 18) and ($len < 26) and (($pos + $len - 6) < $b2[17])) {
                        $mature_bp=count_left_bracket(substr($pattern,($pos - 1),$len));
                        if (abs($mature_bp) > 8) {
                            $proper_check=0;
                            $error_info="loop extends into stem > 8nt: pos = $pos; len = $len";
                            last;
                        }
                    }
                }
                if ($error_info ne ''){last;}
            }
        }
    }
    if ($error_info eq '') {$proper_check = 1}
    $buffer=$buffer."proper_check3 = ".$proper_check."\n";
    $buffer=$buffer."error_info3 = ".$error_info."\n";
    
    if ($proper_check eq 1) {
        if ($need{$id}=~m/$species/) {
            my $save=$b2[-1];
            $b2[-1]=$count;
            print OUT join("\t",@b2),"\t",$weight,"\t",sprintf("%.6f",$count/$b2[2]),"\t",sprintf("%.6f",$weight/$b2[3]),"\t",$mature_pair,"\t",$save,"\n";
        }
        else {
            print OUT $need{$id},"\t",$count,"\t",$weight,"\t",sprintf("%.6f",$count/$b2[2]),"\t",sprintf("%.6f",$weight/$b2[3]),"\t",$mature_pair,"\n";
        }
        print OUT $buffer;
        print OUT join("\t",@mat5),"\t",join("\t",@mat3),"\t",$mature_bp,"\t",$corespond_pos,"\t",$corespond_end,"\t",$error_info,"\n";
    }
    else {
        if ($need{$id}=~m/$species/) {
            my $save=$b2[-1];
            $b2[-1]=$count;
            print OUT2 join("\t",@b2),"\t",$weight,"\t",sprintf("%.6f",$count/$b2[2]),"\t",sprintf("%.6f",$weight/$b2[3]),"\t",$mature_pair,"\t",$save,"\n";
        }
        else {
            print OUT2 $need{$id},"\t",$count,"\t",$weight,"\t",sprintf("%.6f",$count/$b2[2]),"\t",sprintf("%.6f",$weight/$b2[3]),"\t",$mature_pair,"\n";
        }
        print OUT2 $buffer;
        print OUT2 join("\t",@mat5),"\t",join("\t",@mat3),"\t",$mature_bp,"\t",$corespond_pos,"\t",$corespond_end,"\t",$error_info,"\n";
    }
    
}
