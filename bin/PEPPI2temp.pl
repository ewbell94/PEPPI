#!/usr/bin/env perl

use strict;
use warnings;

my $peppidir="!PEPPIDIR!";
my $outdir="!OUTDIR!";
my $maxjobs=!MAXJOBS!;
my $nohomo=0;
my $benchmarkflag=!BENCHMARKFLAG!;
my $batchsize=5;

my $user=`whoami`;
chomp($user);

my $singleflag=0;
$singleflag=1 if (scalar(@ARGV) > 0 && $ARGV[0] eq "s");

print `python $peppidir/bin/splitFUD.py $outdir`;

my @protsA=();
open(my $protcodeA,"<","$outdir/protcodeA.csv");
while (my $line=<$protcodeA>){
    my @parts=split(",",$line);
    push(@protsA,"$parts[0]");
}
close($protcodeA);

my @protsB=();
open(my $protcodeB,"<","$outdir/protcodeB.csv");
while (my $line=<$protcodeB>){
    my @parts=split(",",$line);
    push(@protsB,"$parts[0]");
}
close($protcodeB);

print `mkdir $outdir/PPI`;
#print `mkdir $outdir/hhr`;
#print `mkdir $outdir/model`;
for my $i (0..scalar(@protsA)-1){
    for my $j (0..scalar(@protsB)-1){
	next if (-e "$outdir/PPI/$protsB[$j]-$protsA[$i]" || ($nohomo && $protsB[$j] eq $protsA[$i]));
	my $pairdir="$outdir/PPI/$protsA[$i]-$protsB[$j]";
	print `mkdir $pairdir`;
	print `cp $outdir/fasta/$protsA[$i]/seq.fasta $pairdir/$protsA[$i].seq`;
	my $n=1;
	while (-f "$outdir/fasta/$protsA[$i]/$n.fasta"){
	    print `cp $outdir/fasta/$protsA[$i]/$n.fasta $pairdir/$protsA[$i]\_$n.seq`;
	    $n++;
	}
	print `cp $outdir/fasta/$protsB[$j]/seq.fasta $pairdir/$protsB[$j].seq`;
	$n=1;
	while (-f "$outdir/fasta/$protsB[$j]/$n.fasta"){
	    print `cp $outdir/fasta/$protsB[$j]/$n.fasta $pairdir/$protsB[$j]\_$n.seq`;
	    $n++;
	}
    }
}

my @supported=("SPRING");
my @intset=();
#my @supported=("COTHPPI");
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
	print `$int/$pairname-$prog.pl` if ($singleflag);
=pod
	while (`squeue -u $user | wc -l`-1 >= $maxjobs){
	    print "Queue is currently full, waiting for submission...\n";
	    sleep(60);
	}
	my $jobname="PEPPI_$prog\_$pairname";
	my $errloc="$int/err_$prog.log";
	my $outloc="$int/out_$prog.log";
	print `sbatch -J $jobname -o $outloc -e $errloc -t 20:00:00 $int/$pairname-$prog.pl`;
=cut
    }
    push(@intset,$int) if (!$singleflag);
    if (scalar(@intset) == $batchsize){
	submitBatch(\@intset);
	@intset=();
    }
}

if (scalar(@intset) > 0){
    submitBatch(\@intset);
}


open(my $peppi3script,">","$outdir/PEPPI3.pl");
my $peppi3=`cat $peppidir/bin/PEPPI3temp.pl`;
$peppi3=~s/\!OUTDIR\!/$outdir/;
print $peppi3script $peppi3;
close($peppi3script);
print `chmod +x $outdir/PEPPI3.pl`;

sub submitBatch{
    my @intset=@{$_[0]};
    my $args=join(",",@supported);
    $args="$args ".join(",",@intset);
    
    while (`squeue -u $user | wc -l`-1 >= $maxjobs){
	print "Queue is currently full, waiting for submission...\n";
	sleep(60);
    }
    print `sbatch -J PEPPI2batch -o /dev/null -t 24:00:00 $peppidir/bin/multiwrapper.pl $args`;
   
}
