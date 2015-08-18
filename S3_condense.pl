#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"score file\"  \"blat result(short list as Query)\"  \"output\"" if (@ARGV < 3);
my $filein1=$ARGV[0];   # *score
my $filein2=$ARGV[1];   # mutual blat file
my $fileout=$ARGV[2];
open IN1, $filein1;
open IN2, $filein2;
open OUT, ">".$fileout;
open OUT1, ">".$fileout.".uniq";
open OUT2, ">".$fileout.".log";
my %uniq;
my %expr;
while(<IN1>) {
    chomp;
    my @a=split("\t",$_);
    $uniq{$a[0]}=join("\t",$_);
    $expr{$a[0]}=$a[3];
}
close IN1;
my %save;
my %reject;
my %block;
my %TMP;
my %Redundant;
while(<IN2>) {
    chomp;
    if ((m/^#/) or (m/^$/) or (m/^psL/) or (m/^match/)){next;}
    else {
        my @a=split("\t",$_);
        if (($a[8] eq "+")  and ($a[17] eq 1) and ($a[0] > 40) and ($a[0]/$a[10] >= 0.9)) {
            if (!exists $TMP{$a[9]}{$a[13]}) {$TMP{$a[9]}{$a[13]}=join("\t",@a);}
            else {$Redundant{$a[9]}{$a[13]}=1;}
        }
    }
}
foreach my $Query (keys %TMP) {
    foreach my $Target (keys %{$TMP{$Query}}) {
        if (!exists $Redundant{$Query}{$Target}) {
            print OUT1 $TMP{$Query}{$Target},"\n";
        }
    }
}
close IN2;
close OUT1;
open IN3, $fileout.".uniq";
while(<IN3>) {
    chomp;
    if ((m/^#/) or (m/^$/) or (m/^psL/) or (m/^match/)){next;}
    else {
        my @a=split("\t",$_);
        if (($a[8] eq "+")  and ($a[17] eq 1)) {
            if ($a[9] eq $a[13]) {
                $save{$a[9]}=1;
            }
            else {
                # cluster similar sequences, select the one with more mature reads
                #if (abs($a[10] - $a[14]) <= 5) {
                if ((abs($a[10] - $a[14]) <= 5) and (abs($a[0] - $a[10]) <= 0.1*$a[10])) {
                    if (!exists $expr{$a[13]}) {$save{$a[9]}++;}
                    else {
                        if ($expr{$a[9]} > $expr{$a[13]}) {
                            $save{$a[9]}++;
                            $reject{$a[13]}=$a[9];
                            print OUT2 join("\t","save",$a[9],$expr{$a[9]},"discard",$a[13],$expr{$a[13]}),"\n";
                            print OUT2 join("\t",@a),"\n";
                        }
                        elsif ($expr{$a[9]} < $expr{$a[13]}) {
                            #$save{$a[13]}++;
                            $reject{$a[9]}=$a[13];
                            print OUT2 join("\t","save",$a[13],$expr{$a[13]},"discard",$a[9],$expr{$a[9]}),"\n";
                            print OUT2 join("\t",@a),"\n";
                        }
                        else {
                            my @A=split(/\_\_/,$a[9]);
                            my @B=split(/\_\_/,$a[13]);
                            if ($A[-1] >= $B[-1]) {$save{$a[9]}++;$reject{$a[13]}=$a[9];}
                        }
                    }
                    next;
                }
                if ((exists $reject{$a[9]}) or (exists $reject{$a[13]})) {next;}
                
                #if ((abs($a[0] - $a[10]) <= 5) or (abs($a[0] - $a[14]) <= 5)) {
                #    #  A include B
                #    $save{$a[9]}++;
                #    next;
                #}
                
                # only mature-read similar
                if (($a[0] >= 16) and ($a[0] <= 27)) {
                    my $Query='';
                    if ($a[10] > ($a[11]+$a[12])) {$Query="left";}
                    else{$Query="right";}
                    my $Target='';
                    if ($a[14] > ($a[15]+$a[16])) {$Target="left";}
                    else{$Target="right";}
                    # mature at the same 5p or 3p end
                    if ($Target eq $Query) {
                        $save{$a[9]}++;
                        #$save{$a[13]}++;
                        #print OUT2 join("\t","save",$a[9],$expr{$a[9]},"save",$a[13],$expr{$a[13]}),"\n";
                    }
                    # mature at different end, one is chimera
                    else {
                        # save the block alignments
                        if (exists $block{$a[9]}) {
                            my $tmp=$a[11]."\-".$a[12];
                            if ($block{$a[9]} !~ m/$tmp/) {$block{$a[9]}=$block{$a[9]}."\t".$tmp;}
                        }
                        else {
                            $block{$a[9]}=$a[10]."\t".$a[11]."\-".$a[12];
                        }
                    }
                }
                #else {
                #    # save the block alignments
                #    if (exists $block{$a[9]}) {
                #        #$block{$a[9]}=$block{$a[9]}."\#".$a[11]."\t".$a[12];
                #        my $tmp=$a[11]."\-".$a[12];
                #        if ($block{$a[9]} !~ m/$tmp/) {$block{$a[9]}=$block{$a[9]}."\t".$tmp;}
                #    }
                #    else {
                #        $block{$a[9]}=$a[10]."\t".$a[11]."\-".$a[12];
                #    }
                #}
            }
        }
    }
}
close IN3;
foreach my $id (keys %save) {
    if (exists $reject{$id}) {}
    elsif (exists $block{$id}) {
        print OUT2 $id,"\t",$block{$id},"\n";
        my @a=split("\t",$block{$id});
        my $count=scalar(@a);
        my %span;
        my $f=0;
        for (my $i=1; $i<$count; $i++) {
            $span{$a[$i]}=1;
            my @b=split(/\-/,$a[$i]);
            if (($b[0] < 3) and (($a[0] - $b[1]) > 10)) {
                # find 0-50 and looking for 50-100
                if (exists $span{$b[1]."-".$a[0]}) {$f=50; last;}
                if (exists $span{$b[1]."-".($a[0]-1)}) {$f=50; last;}
                if (exists $span{$b[1]."-".($a[0]-2)}) {$f=50; last;}
            }
            elsif (($b[0] > 10) and (($a[0] - $b[1]) < 3)) {
                # find 50-100 and looking for 0-50
                if (exists $span{"0-".$b[0]}) {$f=100; last;}
                if (exists $span{"1-".$b[0]}) {$f=100; last;}
                if (exists $span{"2-".$b[0]}) {$f=100; last;}
            }
        }
        if ($f > 10) {$reject{$id}=$id;}
    }
}
foreach my $id (keys %save) {
    if (exists $reject{$id}) {}
    else {
        print OUT $uniq{$id},"\n";
    }
}
