#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"_reversed file\"  \"long reads fasta\"  \"long_selected\"  \"structure\"  \"short_map_to_miRBase\"  \"long_map_to_miRBase\"   \"output\"  " if (@ARGV ne 7);
my $reversed=$ARGV[0];
my $fasta=$ARGV[1];
my $select=$ARGV[2];    # selected long reads ID
my $structure=$ARGV[3];
my $maturemap=$ARGV[4];
my $longmap=$ARGV[5];
my $fileout=$ARGV[6];
open IN1,$reversed;
#open IN11,$reversed2;
open IN2,$fasta;
open IN3,$select;
open IN4,$maturemap;
open IN5,$longmap;
open IN6,$structure;
open OUT,">".$fileout;
open OUT2,">".$fileout.".MIR";

sub fill($$$$){   
    my $seq='';
    for(my $i=1; $i<$_[1]; $i++) {$seq=$seq." ";}
    my $tmp=substr($_[2],($_[1]-1),$_[0]);
    $seq=$seq.$tmp;
    #for(my $i=0; $i<$_[0]; $i++) {$seq=$seq."*";}
    for(my $i=($_[0]+$_[1]); $i<=$_[3]; $i++) {$seq=$seq." ";}
    #my @a=split("",$tmp);
    #my $left=0;
    #my $right=0;
    #for(@a){
    #    if($_ eq "("){$left++;}
    #    elsif($_ eq ")"){$right++;}
    #}
    #my $max=$left > $right ? $left : $right;
    #return join("\t",$left,$right,$max,$seq);
    return $seq;
}

my %need;
while(<IN3>) {  # selected score file, store ID and long-read-length only
    chomp;
    my @a=split("\t",$_);
    $need{$a[0]}=$a[0];
}
close IN3;
print "read selected_ID and long_read_length done.\n";

my %RNAfold;
while(<IN6>) {
    chomp;
    if (m/^>/){
        s/^>//;
        my $id=$_;
        my $seq=<IN6>;
        $seq=<IN6>;
        chomp $seq;
        if (exists $need{$id}) {$RNAfold{$id}=$seq;}
    }
}
close IN6;

my %Mature;
while(<IN4>) {  # short reads annotation(miRBase, SOAP2)
    chomp;
    my @a=split("\t",$_);
    if ($a[6] eq "+") {$Mature{$a[0]}=$a[7];}
}
close IN4;

print "read mapped to miRBase mature done.\n";

my %MIR;
while(<IN5>) {  # long reads annotation(miRBase, SOAP2)
    chomp;
    my @a=split("\t",$_);
    if ($a[6] eq "+") {$MIR{$a[0]}=$a[7];}
}
close IN5;

print "read mapped to miRBase precursor done.\n";

my %Fasta;
while(<IN2>) {  # read fasta long reads
    chomp;
    s/^>//;
    s/^@//;
    my @a=split("\t",$_);
    my $seq=<IN2>;
    chomp $seq;
    if (exists $need{$a[0]}) { $Fasta{$a[0]}=$seq; }
}
close IN2;

print "read long read done.\n";

my %uniq;
my %Len;
while(<IN1>) {  # reversed file. FULL mapping pattern
    chomp;
    my @a=split("\t",$_);
    if (exists $need{$a[0]}) {
        $uniq{$a[0]}=$a[4];
        $Len{$a[0]}=length($a[1]);
    }
}
close IN1;

#while(<IN11>) { # reversed.2-like file, detailed mapping information
#    chomp;
#    my @a=split("\t",$_);
#    if (exists $need{$a[0]}) {
#        $need{$a[0]}=join("\t",@a);
#    }
#}
#close IN11;
print "read reversed done.\n";

foreach my $id (keys %need) {
    if (exists $uniq{$id}) {
        my @a=split(/\#/,$uniq{$id});   # reversed mapping info
        my $bseq=$Fasta{$id};
        #my @b2=split("\t",$need{$id});  # reversed.2 detail info
        my %rank;
        for(@a) {
            my @tmp=split(/\|/,$_);
            $rank{$tmp[2]}{$tmp[1]}=$tmp[0];
        }
        my $weight=0;
        my $count=0;
        my $buffer=$Fasta{$id}."\n";
        if (exists $RNAfold{$id}) { $buffer=$buffer.$RNAfold{$id}."\n"; }
        #my $mature=0;
        #my $mature_pair=0;
        my $fit_size=0;
        my $fit_size_weight=0;
        foreach my $pos(sort{$a <=> $b} keys %rank) {
            foreach my $len(sort{$a <=> $b} keys %{$rank{$pos}}) {
                my $result=fill($len,$pos,$bseq,$Len{$id});
                #my ($left,$right,$mpair,$pat)=split("\t",$result);
                #print join("\t",$pat,$left,$right),"\n";
                my @c=split(/\_\_/,$rank{$pos}{$len});
                if (($len >= 17) and ($len <= 25)) {$fit_size++; $fit_size_weight+=$c[1];}
                #if ($c[1] > $mature) {$mature=$c[1]; $mature_pair=sprintf("%.3f",$mpair/$len);}
                if (exists $Mature{$rank{$pos}{$len}}) {
                    #$buffer=$buffer.$pat.$pos."\t".$len."\t".$c[1]."\t".$rank{$pos}{$len}."\t".$Mature{$rank{$pos}{$len}}."\n";
                    $buffer=$buffer.$result."\t".$c[1]."\t".$Mature{$rank{$pos}{$len}}."\n";
                }
                else {
                    #$buffer=$buffer.$pat.$pos."\t".$len."\t".$c[1]."\t".$rank{$pos}{$len}."\n";
                    $buffer=$buffer.$result."\t".$c[1]."\n";
                }
                #print OUT $pat,$pos,"\t",$len,"\t",$c[1],"\t",$rank{$pos}{$len},"\n";
                #if ($left*$right > 0) {$count++; $weight+=$c[1];}
            }
        }
        if ((exists $MIR{$id}) and (($MIR{$id}=~m/mir/) or ($MIR{$id}=~m/let/) or ($MIR{$id}=~m/lin/) or ($MIR{$id}=~m/bantam/))) {
            #my $save=$b2[-1];
            #$b2[-1]=$count;
            print OUT2 $id,"\t",$MIR{$id},"\n";
            #print OUT2 $need{$id},"\t",$MIR{$id},"\n";
            print OUT2 $buffer,"###\n";
        }
        else {
            print OUT $id,"\n";
            #print OUT $need{$id},"\n";
            print OUT $buffer,"###\n";
        }
    }
}
