#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"randfold file\"  \"rnafold file\"  \"long_miRBase\"  \"reversed.2\"  \"output name\" " if (@ARGV ne 5);
my $rand=$ARGV[0];  # .rand
my $rna=$ARGV[1];   # .str
my $kmirid=$ARGV[2];    # long mapped to miRBase, raw file, optional : place a non-existing file name here
my $reversed=$ARGV[3];  # _reversed.2
my $fileout=$ARGV[4];   # scoring
open IN1,$rand;
open IN2,$rna;
if (-e $kmirid) { open IN3,$kmirid; }
open IN4,$reversed;
open OUT,">".$fileout;
open OUT2,">".$fileout.".fail";
open OUT3,">".$fileout.".supportless";

sub EXP($$){
    my $result=1;
    if ($_[1] <= 0) {return $result;}
    else {
        for(my $i=1; $i<=$_[1]; $i++) {
            $result=$result*$_[0];
        }
        return $result;
    }
}

my %ID;
if (-e $kmirid) {
    while(<IN3>) {  # long reads annotation, optional
        chomp;
        my @a=split("\t",$_);
        if ($a[6] eq "+") {$ID{$a[0]}=$a[7]."\_".$a[8];}
    }
}
close IN3;
my %uniq;
while(<IN1>) {  # randfold 
    chomp;
    my @a=split("\t",$_);
    $uniq{$a[0]}=join("\t",$a[1],$a[2]);
}
close IN1;
my %rev;
while(<IN4>) {  # reversed.2 like line-by-line info
    chomp;
    my @a=split("\t",$_);
    $rev{$a[0]}=join("\t",@a);
}
close IN4;
while(<IN2>) {  # rnafold
    chomp;
    s/^>//;
    my $id=$_;
    my $seq=<IN2>;
    chomp $seq;
    my $line=<IN2>;
    # sanity check
    if ((exists $rev{$id}) and (exists $uniq{$id})) {
        # good, proceed
    }
    else {
        next;
    }
    chomp $line;
    my ($bd)=split(" ",$line);
    my @a=split("",$bd);
    my $left=0;
    for(my $i=0; ($a[$i] eq ".")and($i < scalar(@a)-1); $i++){$left++;}
    my $right=0;
    for(my $i=scalar(@a)-1; ($a[$i] eq ".")and($i > 0); $i--){$right++;}
    my $pairs=0;
    my $flag=0;
    my $count=0;
    my @PAIR;   # for (
    $PAIR[0]=0;
    my $flag2=0;
    my @PBIR;   # for )
    $PBIR[0]=0;
    my @PCIR;   # smaller ones between ( and )
    for(my $i=0; $i < scalar(@a); $i++) {
        if ($a[$i] eq ")") {
            if ($flag eq 0) {
                $flag=1;
            }
            $PBIR[$pairs]++;
        }
        elsif ($a[$i] eq "(") {
            #$count++;
            if ($flag eq 1) {
                $flag = 0;
                $pairs++;
            }
            if ($pairs > 0) { $PAIR[$pairs-1]++; }
        }
    }
    for(my $i=0; $i<$pairs; $i++) {
        if ($PAIR[$i] <= $PBIR[$i]) {$PCIR[$i] = $PAIR[$i];}
        else {$PCIR[$i] = $PBIR[$i];}
    }
    if ($pairs eq 0){$PCIR[0]=0;}
    # find the loop borders
    my $loop_left=-10;
    my $loop_right=-10;
    for(my $i=0; $i < scalar(@a); $i++) {
        if ($a[$i] eq ")") {last;}
        elsif ($a[$i] eq "(") {$loop_left=$i;}
    }
    $loop_left++;
    for(my $i=scalar(@a)-1; $i > 0; $i--) {
        if ($a[$i] eq "(") {last;}
        elsif ($a[$i] eq ")") {$loop_right=$i;}
    }
    $loop_right++;
    if (exists $ID{$id}) {
        if ($pairs > 1) {
            print OUT2 join("\t",$id,$uniq{$id},length($seq),$left,$right,($right-$left),$pairs,join("\,",@PCIR),$ID{$id}),"\n",$seq,"\n",$bd,"\n";
            next;
        }
        # only one loop, the mature MUST NOT cover loop
        if (($loop_left > 0) and ($loop_right > 0) and ($pairs eq 0)) {
            my $fail=0;
            my @info=split("\t",$rev{$id});
            # 5p mature 
            if (($info[13] > 0) and (($info[12] + $info[13] - 1 - $loop_left) > 3)){$fail++;}
            # 3p mature
            if (($info[26] > 0) and (($loop_right - $info[25]) > 3)) {$fail++;}
            if ($fail > 0) {
                print OUT2 join("\t",$id,$uniq{$id},length($seq),$left,$right,($right-$left),$pairs,join("\,",@PCIR),$ID{$id}),"\n",$seq,"\n",$bd,"\n";
                print OUT2 "loop_left = ",$loop_left,"\t","loop_right = ",$loop_right,"\n";
                next;
            }
        }
        # two loops, the mature MUST NOT cover shrinked merged loop
        elsif (($loop_left > 0) and ($loop_right > 0) and ($pairs eq 1)) {
            my $fail=0;
            my @info=split("\t",$rev{$id});
            if ($loop_left < 25) {$loop_left = 25;}
            # 5p mature 
            if (($info[13] > 0) and (($info[12] + $info[13] - 1 - $loop_left - 5) > 0)){
                if ($info[12] <= 2) {}
                else {$fail++;}
            }
            # 3p mature
            if (($info[26] > 0) and (($loop_right -5 - $info[25]) > 0)) {
                if (($info[1] - $info[25] - $info[26]) <= 3) {}
                else {$fail++;}
            }
            if ($fail > 0) {
                print OUT2 join("\t",$id,$uniq{$id},length($seq),$left,$right,($right-$left),$pairs,join("\,",@PCIR),$ID{$id}),"\n",$seq,"\n",$bd,"\n";
                print OUT2 "Shrinked_loop_left = ",$loop_left,"\t","Shrinked_loop_right = ",$loop_right,"\n";
                next;
            }
        }
        my $tmp=1;
        if (($left < 4) and ($right < 5)) {$tmp=1;}
        else {
            my @info=split("\t",$rev{$id});
            my $revinfo=$info[0];
            for(my $i=1; $i<=29; $i++) {$revinfo=$revinfo."\t".$info[$i];}
            if ($left >= 4) {
                if ($info[12] <= 5) {    # 5p mature
                    if (($left - $info[12]) < 4) {$tmp=1;}
                    else {$tmp=$left - $info[12] - 4;}
                }
                else {$tmp+=$left - 4;}
            }
            if ($right >= 5) {
                if (($info[1] - $info[25] ) <8) {   # 3p mature
                    if (($right - $info[1] + $info[25]) < 5) {}
                    else {$tmp+=$right - $info[1] + $info[25];}
                }
                else {$tmp+=$right - 4;}
            }
            $tmp=EXP(1.2,$tmp);
        }
        my ($mfe,$pv)=split("\t",$uniq{$id});
        my $score=$mfe/$pv/($pairs+1)/($PCIR[0]*3+1)/$tmp;
        my $sharp=1;
        if (exists $rev{$id}) {
            my @info=split("\t",$rev{$id});
            my $basepair=0;
            my $leftbp=0;
            my $rightbp=0;
            if ($info[13] ne 0) {    # 5p mature length
                my $mat_seq=substr($bd,($info[12]-1),$info[13]);
                my @mat_bp=split("",$mat_seq);
                for(@mat_bp){ if ($_ eq "\(") {$leftbp++;} }
            }
            if ($info[26] ne 0) {   # 3p mature length
                my $mat_seq=substr($bd,($info[25]-1),$info[26]);
                my @mat_bp=split("",$mat_seq);
                for(@mat_bp){ if ($_ eq "\)") {$rightbp++;} }
            }
            $basepair=$leftbp > $rightbp ? $leftbp : $rightbp;
            my $revinfo=$info[0];
            for(my $i=1; $i<=29; $i++) {$revinfo=$revinfo."\t".$info[$i];}
            if ((scalar(@info) eq 32) and ($info[30] ne 'NA') and ($info[31] ne 'NA')) {$sharp=2;}
            my $bonus=1;
            if ($mfe/$info[1] < -20/60 ) {$bonus=2;}
            my $score2=$score*(log($info[3])+1)*($info[14]+$info[27])/$info[3]*$sharp*$bonus;   # 5p and 3p peak height
            #my $score2=$score*(log($info[3])+1)*($info[4]+$info[6])/$info[3]/(1+abs($right-$left-2)/2)*$sharp;
            print OUT join("\t",$revinfo,$uniq{$id},$left,$right,($right-$left),$pairs,join("\,",@PCIR),$score,$score2,$basepair,$ID{$id}),"\n",$seq,"\n",$bd,"\n";
        }
        else {
            print OUT3 join("\t",$id,$uniq{$id},$left,$right,($right-$left),$pairs,join("\,",@PCIR),$score,$ID{$id}),"\n",$seq,"\n",$bd,"\n";
        }
    }
    else {
        if ($pairs > 1) {
            print OUT2 join("\t",$id,$uniq{$id},length($seq),$left,$right,($right-$left),$pairs,join("\,",@PCIR)),"\n",$seq,"\n",$bd,"\n";
            next;
        }
        # only one loop, the mature MUST NOT cover loop
        if (($loop_left > 0) and ($loop_right > 0) and ($pairs eq 0)) {
            my $fail=0;
            my @info=split("\t",$rev{$id});
            # 5p mature 
            if (($info[13] > 0) and (($info[12] + $info[13] - 1 - $loop_left) > 3)){$fail++;}
            # 3p mature
            if (($info[26] > 0) and (($loop_right - $info[25]) > 3)) {$fail++;}
            if ($fail > 0) {
                print OUT2 join("\t",$id,$uniq{$id},length($seq),$left,$right,($right-$left),$pairs,join("\,",@PCIR)),"\n",$seq,"\n",$bd,"\n";
                print OUT2 "loop_left = ",$loop_left,"\t","loop_right = ",$loop_right,"\n";
                next;
            }
        }
        # two loops, the mature MUST NOT cover shrinked merged loop
        elsif (($loop_left > 0) and ($loop_right > 0) and ($pairs eq 1)) {
            my $fail=0;
            my @info=split("\t",$rev{$id});
            if ($loop_left < 25) {$loop_left = 25;}
            # 5p mature 
            if (($info[13] > 0) and (($info[12] + $info[13] - 1 - $loop_left - 5) > 0)){
                if ($info[12] <= 2) {}
                else {$fail++;}
            }
            # 3p mature
            if (($info[26] > 0) and (($loop_right -5 - $info[25]) > 0)) {
                if (($info[1] - $info[25] - $info[26]) <= 3) {}
                else {$fail++;}
            }
            if ($fail > 0) {
                print OUT2 join("\t",$id,$uniq{$id},length($seq),$left,$right,($right-$left),$pairs,join("\,",@PCIR)),"\n",$seq,"\n",$bd,"\n";
                print OUT2 "Shrinked_loop_left = ",$loop_left,"\t","Shrinked_loop_right = ",$loop_right,"\n";
                next;
            }
        }
        my $tmp=1;
        if (($left < 4) and ($right < 5)) {$tmp=1;}
        else {
            my @info=split("\t",$rev{$id});
            my $revinfo=$info[0];
            for(my $i=1; $i<=29; $i++) {$revinfo=$revinfo."\t".$info[$i];}
            if ($left >= 4) {
                if ($info[12] <= 5) {
                    if (($left - $info[12]) < 4) {$tmp=1;}
                    else {$tmp=$left - $info[12] - 4;}
                }
                else {$tmp+=$left - 4;}
            }
            if ($right >= 5) {
                if (($info[1] - $info[25] ) <8) {
                    if (($right - $info[1] + $info[25]) < 5) {}
                    else {$tmp+=$right - $info[1] + $info[25];}
                }
                else {$tmp+=$right - 4;}
            }
            $tmp=EXP(1.2,$tmp);
        }
        my ($mfe,$pv)=split("\t",$uniq{$id});
        my $score=$mfe/$pv/($pairs+1)/($PCIR[0]*3+1)/$tmp;
        my $sharp=1;
        if (exists $rev{$id}) {
            my @info=split("\t",$rev{$id});
            my $basepair=0;
            my $leftbp=0;
            my $rightbp=0;
            if ($info[13] ne 0) {
                my $mat_seq=substr($bd,($info[12]-1),$info[13]);
                my @mat_bp=split("",$mat_seq);
                for(@mat_bp){ if ($_ eq "\(") {$leftbp++;} }
            }
            if ($info[26] ne 0) {
                my $mat_seq=substr($bd,($info[25]-1),$info[26]);
                my @mat_bp=split("",$mat_seq);
                for(@mat_bp){ if ($_ eq "\)") {$rightbp++;} }
            }
            $basepair=$leftbp > $rightbp ? $leftbp : $rightbp;
            my $revinfo=$info[0];
            for(my $i=1; $i<=29; $i++) {$revinfo=$revinfo."\t".$info[$i];}
            if ((scalar(@info) eq 32) and ($info[30] ne 'NA') and ($info[31] ne 'NA')) {$sharp=2;}
            my $bonus=1;
            if ($mfe/$info[1] < -20/60 ) {$bonus=2;}
            my $score2=$score*(log($info[3])+1)*($info[14]+$info[27])/$info[3]*$sharp*$bonus;
            #my $score2=$score*(log($info[3])+1)*($info[4]+$info[6])/$info[3]/(1+abs($right-$left-2)/2)*$sharp;
            print OUT join("\t",$revinfo,$uniq{$id},$left,$right,($right-$left),$pairs,join("\,",@PCIR),$score,$score2,$basepair),"\n",$seq,"\n",$bd,"\n";
        }
        else {
            print OUT3 join("\t",$id,$uniq{$id},$left,$right,($right-$left),$pairs,join("\,",@PCIR),$score,$ID{$id}),"\n",$seq,"\n",$bd,"\n";
        }
    }
}
