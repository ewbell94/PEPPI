#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min);
use Time::HiRes qw(time);

my $peppidir="/nfs/amino-home/ewbell/PEPPI";
my $outdir="/nfs/amino-home/ewbell/PEPPI/test/PPI";
my $pairname="prot1-prot4";
my $user=`whoami`;
chomp($user);

my $randomTag=int(rand(1000000));
#my $uniprotdb="/nfs/amino-library/local/hhsuite/uniprot20_2016_02/uniprot20_2016_02";
#my $blitsdb="$peppidir/lib/SEQ/hhblitsdb/db";
my $interactfile="$peppidir/lib/SEQ/psicquic.txt";
my $aliasfile="$peppidir/lib/SEQ/100_psicquic.fasta.clstr.aliases";

my $tempdir="/tmp/$user/PEPPI_SEQ_$pairname\_$randomTag";
if (! -e "$tempdir"){
    print `mkdir -p $tempdir`;
} else {
    print `rm -rf $tempdir/*`;
}
#print `cp $interactfile $tempdir/ints.txt`;

my $maxtargets=200;
my $tophits=100;
my @chains=split("-",$pairname);

print `mkdir -p $outdir/$pairname/SEQ`;
open(my $resfile,">","$outdir/$pairname/SEQ/res.txt");
print `date`;

blastSearch();
#psiNW();

print `sync`;
print `rm -rf $tempdir`;
print `date`;

sub blastSearch{
    my $t0=time();
    my $blast1res=`cat $outdir/../mono/$chains[0]/$chains[0].seq`;
    print $blast1res;
    if ($blast1res eq ""){
	print "BLAST for sequence A returned no hits in database.\n";
	print $resfile "0.0\n";
	print `echo "$pairname,0.0" >> $outdir/SEQres.txt`;
	print `rm -rf $tempdir`;
	exit(1);
    }
    print "\n";
    
    my $blast2res=`cat $outdir/../mono/$chains[1]/$chains[1].seq`;
    print $blast2res;
    if ($blast2res eq ""){
	print "BLAST for sequence B returned no hits in STRING.\n";
	print $resfile "0.0\n";
	print `echo "$pairname,0.0" >> $outdir/SEQres.txt`;
	print `rm -rf $tempdir`;
	exit(1);
    }
    print "\n";
    
    my @ids1=extractID($blast1res);
    my @ids2=extractID($blast2res);
    my @results=();

    my $t1=time();
    print "Total id extraction time: ".($t1-$t0)."\n";
    open(my $complexfile,"<",$interactfile);
    my %complexlist=(); #Find all potential partner chains given some core chain
    while (my $line=<$complexfile>){
	chomp($line);
	my @chains=split(" ",$line);
	if (exists($complexlist{$chains[0]})){
	    push(@{$complexlist{$chains[0]}},$chains[1]);
	} else {
	    my @value=($chains[1]);
	    $complexlist{$chains[0]}=\@value;
	}

	if (exists($complexlist{$chains[1]})){
	    push(@{$complexlist{$chains[1]}},$chains[0]);
	} else {
	    my @value=($chains[0]);
	    $complexlist{$chains[1]}=\@value;
	}
    }
    close($complexfile);

    for my $i (0..scalar(@ids1)-1){
	my $id1=$ids1[$i][0];
	my @partners=();
	if (exists($complexlist{$id1})){
	    @partners=@{$complexlist{$id1}};
	} else {
	    next;
	}
	for my $j (0..scalar(@ids2)-1){
	    my $id2=$ids2[$j][0];
	    my $pflag=0;
	    for my $partner (@partners){
		if ($id2 eq $partner){
		    $pflag=1;
		    last;
		}
	    }
	    if ($pflag){
		print "Dimer found!\n";
		my $seqid1=$ids1[$i][1];
		my $seqid2=$ids2[$j][1];
		my @topush=($id1,$id2,2/(1/$seqid1+1/$seqid2),$seqid1,$seqid2);
		push(@results,\@topush);
	    }
	}
    }
    @results=sort{$b->[2]<=>$a->[2]} @results;
    open(my $summary,">","$outdir/$pairname/SEQ/SearchSummary.txt");
    for my $i (0..min(scalar(@results),$tophits)-1){
	print $summary "$results[$i][0]\t$results[$i][1]\t$results[$i][2]\t$results[$i][3]\t$results[$i][4]\n";
    }
    if (scalar(@results)>0){
	print $resfile "$results[0][2]\n";
	print `echo "$pairname,$results[0][2]" >> $outdir/SEQres.txt`;
    } else {
	print $resfile "0.0\n";
	print `echo "$pairname,0.0" >> $outdir/SEQres.txt`;
    }
    my $t2=time();
    print "Total search time: ".($t2-$t1)."\n";
}

=pod
sub psiNW{
    print `$peppidir/bin/hhsuite/bin/hhblits -i $outdir/$pairname/$chains[0].seq -oa3m $tempdir/$chains[0].a3m -d $uniprotdb -n 2 -e 0.001`;
    print `$peppidir/bin/hhsuite/bin/hhblits -i $outdir/$pairname/$chains[1].seq -oa3m $tempdir/$chains[1].a3m -d $uniprotdb -n 2 -e 0.001`;
    print `$peppidir/bin/hhsuite/bin/hhblits -i $tempdir/$chains[0].a3m -d $blitsdb`;
    print `$peppidir/bin/hhsuite/bin/hhblits -i $tempdir/$chains[1].a3m -d $blitsdb`;
}
=cut

sub extractID{
    my $blastres=$_[0];
    my @blines=split("\n",$blastres);
    my @ids=();

    my %aliases=();
    open(my $aliasfp,"<","$aliasfile");
    while (my $line=<$aliasfp>){
	chomp($line);
	my @bigparts=split(";",$line);
	my @littleparts=split(",",$bigparts[1]);
	$aliases{$bigparts[0]}=\@littleparts;
    }

    for my $i (0..scalar(@blines)-1){
	my @bparts=split(" ",$blines[$i]);
	my $skipflag=0;
	for my $j (0..scalar(@ids)-1){
	    if ($bparts[0] eq $ids[$j][0]){
		$skipflag=1;
		last;
	    }
	}
	next if ($skipflag);
	my $seqid=min($bparts[1]/$bparts[2],$bparts[1]/$bparts[3]);
	for my $alias (@{$aliases{$bparts[0]}}){
	    my @pair=($alias,$seqid);
	    push(@ids,\@pair);
	}
    }

    return @ids;
}