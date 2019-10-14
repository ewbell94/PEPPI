#!/usr/bin/env perl

use strict;
use warnings;

my $peppidir="!PEPPIDIR!";
my $outdir="!OUTDIR!";
my $hpc=!HPC!;
my $maxjobs=!MAXJOBS!;
my $nohomo=1;
my $benchmarkflag=0;

my $user=`whoami`;
chomp($user);

print `python $peppidir/bin/splitFUD.py $outdir`;

my @domainsA=();
open(my $protcodeA,"<","$outdir/protcodeA.csv");
while (my $line=<$protcodeA>){
    my @parts=split(",",$line);
    my $prot=$parts[0];
    my $i=1;
    while (-e "$outdir/fasta/$prot/$i.fasta"){
	push(@domainsA,"$prot\_$i");
	$i++;
    }
}
close($protcodeA);

my @domainsB=();
open(my $protcodeB,"<","$outdir/protcodeB.csv");
while (my $line=<$protcodeB>){
    my @parts=split(",",$line);
    my $prot=$parts[0];
    my $i=1;
    while (-e "$outdir/fasta/$prot/$i.fasta"){
	push(@domainsB,"$prot\_$i");
	$i++;
    }
}
close($protcodeB);

print `mkdir $outdir/PPI`;
print `mkdir $outdir/hhr`;
for my $i (0..scalar(@domainsA)-1){
    for my $j (0..scalar(@domainsB)-1){
	next if (-e "$outdir/PPI/$domainsB[$j]-$domainsA[$i]" || ($nohomo && $domainsB[$j] eq $domainsA[$i]));
	my $pairdir="$outdir/PPI/$domainsA[$i]-$domainsB[$j]";
	print `mkdir $pairdir`;
	my @domainparts=split("_",$domainsA[$i]);
	print `cp $outdir/fasta/$domainparts[0]/$domainparts[1].fasta $pairdir/$domainsA[$i].seq`;
	@domainparts=split("_",$domainsB[$j]);
	print `cp $outdir/fasta/$domainparts[0]/$domainparts[1].fasta $pairdir/$domainsB[$j].seq`;
    }
}

my @supported=("SPRING","STRING");

for my $int (glob("$outdir/PPI/*/")){
    my @parts=split("/",$int);
    my $pairname=$parts[-1];
    for my $prog (@supported){
	my $modtext=`cat $peppidir/bin/${prog}mod`;
	$modtext=~s/\!PEPPIDIR\!/$peppidir/;
	$modtext=~s/\!OUTDIR\!/$outdir/;
	$modtext=~s/\!PAIRNAME\!/$pairname/;
	$modtext=~s/\!BENCHMARK\!/$benchmarkflag/;
	open(my $jobscript,">","$int/$pairname-$prog.pl");
	print $jobscript $modtext;
	close($jobscript);
	print `chmod +x $int/$pairname-$prog.pl`;
	if ($hpc){
	    my $jobname="PEPPI_$prog\_$pairname";
	    my $errloc="$int/err_$prog.log";
	    my $outloc="$int/out_$prog.log";
	    print `qsub -e $errloc -o $outloc -N $jobname -l walltime="24:00:00" $int/$pairname-$prog.pl`;
	} else {
	    print `$int/$pairname-$prog.pl`;
	}
    }
}

open(my $peppi3script,">","$outdir/PEPPI3.pl");
my $peppi3=`cat $peppidir/bin/PEPPI3temp.pl`;
$peppi3=~s/\!OUTDIR\!/$outdir/;
print $peppi3script $peppi3;
close($peppi3script);
print `chmod +x $outdir/PEPPI3.pl`;
