#!/usr/bin/env perl

use strict;
use warnings;

my $peppidir="/nfs/amino-home/ewbell/PEPPI";
my $outdir="/nfs/amino-home/ewbell/PEPPI/PEPPI";
my $hpc=0;
my $maxjobs=300;

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
	my $pairdir="$outdir/PPI/$domainsA[$i]-$domainsB[$j]";
	print `mkdir $pairdir`;
	my @domainparts=split("_",$domainsA[$i]);
	print `cp $outdir/fasta/$domainparts[0]/$domainparts[1].fasta $pairdir/$domainsA[$i].seq`;
	@domainparts=split("_",$domainsB[$j]);
	print `cp $outdir/fasta/$domainparts[0]/$domainparts[1].fasta $pairdir/$domainsB[$j].seq`;
    }
}

my @supported=("SPRING","STRING");
for my $prog (@supported){
    my $modtext=`cat $peppidir/bin/${prog}mod`;
    $modtext=~s/\!PEPPIDIR\!/$peppidir/;
}

=pod
my $batchnum=1;
my @querylist=();
for my $pair (glob("$outdir/PPI/*/")){
    next if (-e "$pair/SPRING/TemplateSummary.txt");
    my @pairparts=split("/",$pair);
    my $pairname=$pairparts[-1];
    push(@querylist,$pairname);
    if (scalar(@querylist) >= $batchsize){
	my $args=join(' ',@querylist);
	if ($hpc){
	    my $jobname="SPRINGbatch$batchnum";
	    my $errloc="$outdir/SPRING/SPRINGbatch$batchnum\_err.log";
	    my $outloc="$outdir/SPRING/SPRINGbatch$batchnum\_out.log";
	    if ($batchsize == 1){
		$jobname="SPRING_$args";
		$errloc="$pair/err.log";
		$outloc="$pair/out.log";
	    }
	    while (`qstat -u $user | wc -l`-5 >= $maxjobs){
		sleep(60);
	    }
	    print `qsub -e $errloc -o $outloc -l walltime="24:00:00" -N $jobname $peppidir/bin/springwrapper.pl -F "$springdir $outdir/SPRING $args"`; 
	} else {
	    print `$peppidir/bin/springwrapper.pl $springdir $outdir/SPRING $args`;
	}
	@querylist=();
	$batchnum++;
    }
}
if (scalar(@querylist) > 0){
    my $args=join(' ',@querylist);
    if ($hpc){
	my $jobname="SPRINGbatch$batchnum";
	my $errloc="$outdir/SPRING/SPRINGbatch$batchnum\_err.log";
	my $outloc="$outdir/SPRING/SPRINGbatch$batchnum\_out.log";
	while (`qstat -u $user | wc -l`-5 >= $maxjobs){
	    sleep(60);
	}
	print `qsub -e $errloc -o $outloc -l walltime="24:00:00" -N $jobname $peppidir/bin/springwrapper.pl -F "$springdir $outdir/SPRING $args"`;
    } else {
	print `$peppidir/bin/springwrapper.pl $springdir $outdir/SPRING $args`;
    }
}
=cut
open(my $peppi3script,">","$outdir/PEPPI3.pl");
my $peppi3=`cat $peppidir/bin/PEPPI3temp.pl`;
$peppi3=~s/\!OUTDIR\!/$outdir/;
print $peppi3script $peppi3;
close($peppi3script);
print `chmod +x $outdir/PEPPI3.pl`;
