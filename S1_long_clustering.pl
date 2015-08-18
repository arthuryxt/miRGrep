#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"*_reversed.2.cfs\"  \"long_reads.fa\"  \"output basename\" " if (@ARGV ne 3);
my $filein=$ARGV[0];
my $fasta=$ARGV[1];
my $fileout=$ARGV[2];
open IN,$filein;    	# S1_reversed.2.cfs
open OUT,">".$fileout;  # S1_reversed.2.cfsc
open OUT2,">".$fileout."_reject";
open OUT3,">".$fileout."_mid";
my %matboth;
my %matrev;
my %uniq;
my %mat3p;
my %mat5p;

sub levenshtein($$){
    my @A=split //, lc shift; # lower case
    my @B=split //, lc shift;
    my @W=(0..@B);
    my ($i, $j, $cur, $next);
    for $i (0..$#A){
	$cur=$i+1;
	for $j (0..$#B){
	    $next=min( $W[$j+1]+1, $cur+1, ($A[$i] ne $B[$j])+$W[$j] );
            $W[$j]=$cur;
            $cur=$next;
	}
	$W[@B]=$next;
    }
    return $next;
}

sub min($$$){
    if ($_[0] < $_[2]){ pop @_; } else { shift @_; }
    return $_[0] < $_[1]? $_[0]:$_[1];
}

while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    $uniq{$a[0]}=1;
    if (($a[30] ne "NA") and ($a[31] ne "NA")) {
        # both matures exist
        if (exists $matboth{$a[30]}{$a[31]}) {
            my @b=split("\t",$matboth{$a[30]}{$a[31]});
	    
            if ($a[3] > 1.5*$b[3]) {
		print OUT2 "1#discard the second one, which has much less mature reads.\n";
		print OUT2 join("\t",@a),"\n";
                print OUT2 $matboth{$a[30]}{$a[31]},"\n";
		$matboth{$a[30]}{$a[31]}=join("\t",@a);
		$matrev{$a[31]}{$a[30]}=join("\t",@a);
	    }
            elsif ($a[3]*1.5 < $b[3]) {
		print OUT2 "2#discard the second one, which has much less mature reads.\n";
                print OUT2 $matboth{$a[30]}{$a[31]},"\n";
                print OUT2 join("\t",@a),"\n";
            }
            else {
                my $pointa=0;
		my $pointb=0;
                if ($a[12] < $b[12]) {$pointa++;}	# 5p mature start - 0
		elsif ($a[12] > $b[12]) {$pointb++;}
                if (($a[1] - $a[19]) < ($b[1] - $b[19])) {$pointa++;}	# prec_length - 3p mature end
		elsif (($a[1] - $a[19]) > ($b[1] - $b[19])) {$pointb++;}
		#print join("\t",@a),"\n",$pointa,"\n",$matboth{$a[30]}{$a[31]},"\n",$pointb,"\n";
                if ($pointa > $pointb) {
		    print OUT2 "3#discard the second one, which has worse support mature reads.\n";
		    print OUT2 join("\t",@a),"\n";
		    print OUT2 $matboth{$a[30]}{$a[31]},"\n";
		    $matboth{$a[30]}{$a[31]}=join("\t",@a);
		    $matrev{$a[31]}{$a[30]}=join("\t",@a);
		}
                elsif ($pointa < $pointb) {
		    print OUT2 "4#discard the second one, which has worse support mature reads.\n";
		    print OUT2 $matboth{$a[30]}{$a[31]},"\n";
		    print OUT2 join("\t",@a),"\n";
                }
		else {
		    $pointa=0;
		    $pointb=0;
		    if ($a[12] < $b[12]) {$pointa=$b[12] - $a[12];}
		    elsif ($a[12] > $b[12]) {$pointb=$a[12] - $b[12];}
		    if (($a[1] - $a[19]) < ($b[1] - $b[19])) {$pointa+=($b[1] - $b[19] - $a[1] + $a[19])/2;}
		    elsif (($a[1] - $a[19]) > ($b[1] - $b[19])) {$pointb+=($a[1] - $a[19] - $b[1] + $b[19])/2;}
		    if ($pointa > $pointb) {
			print OUT2 "5#discard the second one, which has a little worse support mature reads.\n";
			print OUT2 join("\t",@a),"\n";
			print OUT2 $matboth{$a[30]}{$a[31]},"\n";
			$matboth{$a[30]}{$a[31]}=join("\t",@a);
			$matrev{$a[31]}{$a[30]}=join("\t",@a);
		    }
		    elsif ($pointa < $pointb) {
			print OUT2 "6#discard the second one, which has a little worse support mature reads.\n";
			print OUT2 $matboth{$a[30]}{$a[31]},"\n";
			print OUT2 join("\t",@a),"\n";
		    }
		    elsif ($a[3] > $b[3]) {
			print OUT2 "7#discard the second one, slightly worse.\n";
			print OUT2 join("\t",@a),"\n";
			print OUT2 $matboth{$a[30]}{$a[31]},"\n";
			$matboth{$a[30]}{$a[31]}=join("\t",@a);
			$matrev{$a[31]}{$a[30]}=join("\t",@a);
		    }
		    elsif($a[3] < $b[3]) {
			#print $a[3],"\t",$b[3],"\n";;
			print OUT2 "8#discard the second one, slightly worse.\n";
			print OUT2 $matboth{$a[30]}{$a[31]},"\n";
			print OUT2 join("\t",@a),"\n";
			#$matboth{$a[30]}{$a[31]}=join("\t",@a);
			#$matrev{$a[31]}{$a[30]}=join("\t",@a);
		    }
		    elsif($a[1] > $b[1]) {
			# given solexa produce deletions more likely than incertions, choose the longer
			print OUT2 "8a#discard the second one, slightly worse.\n";
			print OUT2 join("\t",@a),"\n";
			print OUT2 $matboth{$a[30]}{$a[31]},"\n";
			$matboth{$a[30]}{$a[31]}=join("\t",@a);
			$matrev{$a[31]}{$a[30]}=join("\t",@a);
		    }
		    else {
			print OUT2 "8b#discard the second one, slightly worse.\n";
			print OUT2 $matboth{$a[30]}{$a[31]},"\n";
			print OUT2 join("\t",@a),"\n";
		    }
		}
            }
        }
        else {
            $matboth{$a[30]}{$a[31]}=join("\t",@a);
            $matrev{$a[31]}{$a[30]}=join("\t",@a);
        }
    }
    elsif ($a[30] ne "NA") {
        # only 5p-mature exists
	$mat5p{$a[30]}{$a[0]}=join("\t",@a);
    }
    elsif ($a[31] ne "NA") {
        # only 3p-mature exists
        $mat3p{$a[31]}{$a[0]}=join("\t",@a);
    }
}
close IN;

#foreach my $id (keys %mat5p) {
#    print OUT $mat5p{$id},"\n";
#}
#foreach my $id (keys %mat3p) {
#    print OUT $mat3p{$id},"\n";
#}

open IN1,$fasta;
while(<IN1>) {
    chomp;
    s/^>//;
    s/^@//;
    my $id=$_;
    my $seq=<IN1>;
    if (exists $uniq{$id}) {
        chomp $seq;
        $uniq{$id}=$seq;
    }   
}
close IN1;
foreach my $left (keys %matboth) {
    foreach my $right (keys %{$matboth{$left}}) {
        print OUT3 $matboth{$left}{$right},"\n";
    }
}
# merge all long reads with the same 5p-mature
foreach my $left (keys %matboth) {
    if ($left eq "NA") {next;}
    my @a;
    my $count=0;
    foreach my $right (sort{$a cmp $b} keys %{$matboth{$left}}) {
        $a[$count]=$left."\t".$right;
        $count++;
    }
    $count--;
    if ($count eq 0){
        my @b=split("\t",$a[0]);
        print OUT $matboth{$b[0]}{$b[1]},"\n";
        next;
    }
    for(my $i=0; $i<$count; $i++) {
        for(my $j=$i+1; $j<=$count; $j++) {
            if (($a[$i] eq '') or ($a[$j] eq '')) {next;}
            my @first=split("\t",$a[$i]);
            my @second=split("\t",$a[$j]);
            my @F=split("\t",$matboth{$first[0]}{$first[1]});
            my @S=split("\t",$matboth{$second[0]}{$second[1]});
            my $dist=levenshtein($uniq{$F[0]},$uniq{$S[0]}) - abs(length($uniq{$F[0]}) - length($uniq{$S[0]}));
            my $len=min(length($uniq{$F[0]}), length($uniq{$S[0]}), 1000);
            if (($dist/$len) < 0.1) {
                # similar sequence found, find the one with more matures and delete the one with less
                if ($F[24] > $S[24]) {
                    $a[$j]='';
                    print OUT2 "9#similar sequence found, delete the second(less mature reads)\n";
                    print OUT2 $matboth{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                    print OUT2 $matboth{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                }
                elsif ($F[24] < $S[24]) {
                    $a[$i]='';
                    print OUT2 "10#similar sequence found, delete the second(less mature reads)\n";
                    print OUT2 $matboth{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                    print OUT2 $matboth{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                }
                else {
		    #my $pointa=0;
		    #my $pointb=0;
		    #if ($F[8] < $S[8]) {$pointa++;}
		    #elsif ($F[8] > $S[8]) {$pointb++;}
		    #if (($F[1] - $F[19]) < ($S[1] - $S[19])) {$pointa++;}
		    #elsif (($F[1] - $F[19]) > ($S[1] - $S[19])) {$pointb++;}
                    my $points=0;
                    if ($F[8] < $S[8]) {$points++;}
                    if (($F[1] - $F[19]) < ($S[1] - $S[19])) {$points++;}
                    if ($points > 0) {
                        $a[$j]='';
                        print OUT2 "11#similar sequence found, delete the second(less mature reads)\n";
                        print OUT2 $matboth{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                        print OUT2 $matboth{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                    }
                    else {
                        $a[$i]='';
                        print OUT2 "12#similar sequence found, delete the second(less mature reads)\n";
                        print OUT2 $matboth{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                        print OUT2 $matboth{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                    }
                }
            }
        }
    }
    for(my $i=0; $i<=$count; $i++) {
        if ($a[$i] ne '') {
            my @b=split("\t",$a[$i]);
            print OUT $matboth{$b[0]}{$b[1]},"\n";
        }
    }
}

foreach my $left (keys %matrev) {
    if ($left eq "NA") {next;}
    my @a;
    my $count=0;
    foreach my $right (sort{$a cmp $b} keys %{$matrev{$left}}) {
        $a[$count]=$left."\t".$right;
        $count++;
    }
    $count--;
    if ($count eq 0){
        my @b=split("\t",$a[0]);
        print OUT $matrev{$b[0]}{$b[1]},"\n";
        next;
    }
    for(my $i=0; $i<$count; $i++) {
        for(my $j=$i+1; $j<=$count; $j++) {
            if (($a[$i] eq '') or ($a[$j] eq '')) {next;}
            my @first=split("\t",$a[$i]);
            my @second=split("\t",$a[$j]);
            my @F=split("\t",$matrev{$first[0]}{$first[1]});
            my @S=split("\t",$matrev{$second[0]}{$second[1]});
            my $dist=levenshtein($uniq{$F[0]},$uniq{$S[0]}) - abs(length($uniq{$F[0]}) - length($uniq{$S[0]}));
            my $len=min(length($uniq{$F[0]}), length($uniq{$S[0]}), 1000);
            if (($dist/$len) < 0.1) {
                # similar sequence found, find the one with more matures and delete the one with less
                if ($F[13] > $S[13]) {
                    $a[$j]='';
                    print OUT2 "13#similar sequence found, delete the second(less mature reads)\n";
                    print OUT2 $matrev{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                    print OUT2 $matrev{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                }
                elsif ($F[13] < $S[13]) {
                    $a[$i]='';
                    print OUT2 "14#similar sequence found, delete the second(less mature reads)\n";
                    print OUT2 $matrev{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                    print OUT2 $matrev{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                }
                else {
                    my $points=0;
                    if ($F[8] < $S[8]) {$points++;}
                    if (($F[1] - $F[19]) < ($S[1] - $S[19])) {$points++;}
                    if ($points > 0) {
                        $a[$j]='';
                        print OUT2 "15#similar sequence found, delete the second(less mature reads)\n";
                        print OUT2 $matrev{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                        print OUT2 $matrev{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                    }
                    else {
                        $a[$i]='';
                        print OUT2 "16#similar sequence found, delete the second(less mature reads)\n";
                        print OUT2 $matrev{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                        print OUT2 $matrev{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                    }
                }
            }
        }
    }
    for(my $i=0; $i<=$count; $i++) {
        if ($a[$i] ne '') {
            my @b=split("\t",$a[$i]);
            print OUT $matrev{$b[0]}{$b[1]},"\n";
        }
    }
}

foreach my $left(keys %mat5p) {
    my @a;
    my $count=0;
    foreach my $right (sort{$a cmp $b} keys %{$mat5p{$left}}) {
        $a[$count]=$left."\t".$right;
        $count++;
    }
    $count--;
    if ($count eq 0){
        my @b=split("\t",$a[0]);
        print OUT $mat5p{$b[0]}{$b[1]},"\n";
        next;
    }
    for(my $i=0; $i<$count; $i++) {
        for(my $j=$i+1; $j<=$count; $j++) {
            if (($a[$i] eq '') or ($a[$j] eq '')) {next;}
            my @first=split("\t",$a[$i]);
            my @second=split("\t",$a[$j]);
            my @F=split("\t",$mat5p{$first[0]}{$first[1]});
            my @S=split("\t",$mat5p{$second[0]}{$second[1]});
            my $dist=levenshtein($uniq{$F[0]},$uniq{$S[0]}) - abs(length($uniq{$F[0]}) - length($uniq{$S[0]}));
            my $len=min(length($uniq{$F[0]}), length($uniq{$S[0]}), 1000);
            if (($dist/$len) < 0.1) {
                # similar sequence found, save the longer
                if ($F[1] > $S[1]) {
                    $a[$j]='';
                    print OUT2 "17#similar sequence found, delete the second(less mature reads)\n";
                    print OUT2 $mat5p{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                    print OUT2 $mat5p{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                }
                else {
                    $a[$j]='';
                    print OUT2 "18#similar sequence found, delete the second(less mature reads)\n";
                    print OUT2 $mat5p{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                    print OUT2 $mat5p{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                }
            }
        }
    }
    for(my $i=0; $i<=$count; $i++) {
        if ($a[$i] ne '') {
            my @b=split("\t",$a[$i]);
            print OUT $mat5p{$b[0]}{$b[1]},"\n";
        }
    }
}

foreach my $left(keys %mat3p) {
    my @a;
    my $count=0;
    foreach my $right (sort{$a cmp $b} keys %{$mat3p{$left}}) {
        $a[$count]=$left."\t".$right;
        $count++;
    }
    $count--;
    if ($count eq 0){
        my @b=split("\t",$a[0]);
        print OUT $mat3p{$b[0]}{$b[1]},"\n";
        next;
    }
    for(my $i=0; $i<$count; $i++) {
        for(my $j=$i+1; $j<=$count; $j++) {
            if (($a[$i] eq '') or ($a[$j] eq '')) {next;}
            my @first=split("\t",$a[$i]);
            my @second=split("\t",$a[$j]);
            my @F=split("\t",$mat3p{$first[0]}{$first[1]});
            my @S=split("\t",$mat3p{$second[0]}{$second[1]});
            my $dist=levenshtein($uniq{$F[0]},$uniq{$S[0]}) - abs(length($uniq{$F[0]}) - length($uniq{$S[0]}));
            my $len=min(length($uniq{$F[0]}), length($uniq{$S[0]}), 1000);
            if (($dist/$len) < 0.1) {
                # similar sequence found, save the longer
                if ($F[1] > $S[1]) {
                    $a[$j]='';
                    print OUT2 "19#similar sequence found, delete the second(less mature reads)\n";
                    print OUT2 $mat3p{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                    print OUT2 $mat3p{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                }
                else {
                    $a[$j]='';
                    print OUT2 "20#similar sequence found, delete the second(less mature reads)\n";
                    print OUT2 $mat3p{$first[0]}{$first[1]},"\n",$uniq{$F[0]},"\n";
                    print OUT2 $mat3p{$second[0]}{$second[1]},"\n",$uniq{$S[0]},"\n";
                }
            }
        }
    }
    for(my $i=0; $i<=$count; $i++) {
        if ($a[$i] ne '') {
            my @b=split("\t",$a[$i]);
            print OUT $mat3p{$b[0]}{$b[1]},"\n";
        }
    }
}