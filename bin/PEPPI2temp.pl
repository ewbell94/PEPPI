#!/usr/bin/env perl

use strict;
use warnings;

my $peppidir="!PEPPIDIR!";
my $outdir="!OUTDIR!";
my $maxjobs=!MAXJOBS!;
my $nohomo=0;
my $keepflag=1;
my $benchmarkflag=!BENCHMARKFLAG!;
my $batchsize=1;

my $user=`whoami`;
chomp($user);

my $singleflag=0;
$singleflag=1 if (scalar(@ARGV) > 0 && $ARGV[0] eq "s");

#load all proteins from the protcodes and number domains in order
my @protsA=();
open(my $protcodeA,"<","$outdir/protcodeA.csv");
while (my $line=<$protcodeA>){
    my @parts=split(",",$line);
    push(@protsA,"$parts[0]");
    my @domains=treeSearch($parts[0],"$outdir/mono/$parts[0]");
    for my $i (0..scalar(@domains)-1){
	print "$domains[$i]\n";
	my $n=$i+1;
	print `cp $outdir/mono/$parts[0]/$domains[$i].fasta $outdir/mono/$parts[0]/$parts[0]\_$n.fasta`;
	print `cp $outdir/mono/$parts[0]/$domains[$i].hhr.gz $outdir/mono/$parts[0]/$parts[0]\_$n.hhr.gz`;
    }
}
close($protcodeA);

my @protsB=();
open(my $protcodeB,"<","$outdir/protcodeB.csv");
while (my $line=<$protcodeB>){
    my @parts=split(",",$line);
    push(@protsB,"$parts[0]");
    my @domains=treeSearch($parts[0],"$outdir/mono/$parts[0]");
    for my $i (0..scalar(@domains)-1){
        my $n=$i+1;
	print "$domains[$i]\n";
        print `cp $outdir/mono/$parts[0]/$domains[$i].fasta $outdir/mono/$parts[0]/$parts[0]\_$n.fasta`;
	print `cp $outdir/mono/$parts[0]/$domains[$i].hhr.gz $outdir/mono/$parts[0]/$parts[0]\_$n.hhr.gz`;
    }
}
close($protcodeB);

my @supported=("AC"); #Change this to change which modules are available
my @intset=();

#Create PPI folders and run the pairs
print `mkdir -p $outdir/PPI`;
for my $i (0..scalar(@protsA)-1){
    for my $j (0..scalar(@protsB)-1){
	next if ((-e "$outdir/PPI/$protsB[$j]-$protsA[$i]" && $protsA[$i] ne $protsB[$j]) || ($nohomo && $protsB[$j] eq $protsA[$i]));
	my $pairdir="$outdir/PPI/$protsA[$i]-$protsB[$j]";
	print `mkdir -p $pairdir`;
	my $pairname="$protsA[$i]-$protsB[$j]";
	my $int="$outdir/PPI/$pairname";
	print `echo "$int" > $int/out.log` if ($singleflag);
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
	    if ($singleflag){
		print `echo "$prog" >> $int/out.log`;
		print `$int/$pairname-$prog.pl >> $int/out.log`;
	    }
            while (`squeue -u $user | wc -l`-1 >= $maxjobs){
                print "Queue is currently full, waiting for submission...\n";
                sleep(60);
            }
        }
        push(@intset,$int) if (!$singleflag);
        if (scalar(@intset) == $batchsize){
            submitBatch(\@intset);
            @intset=();
        }
    }
}

if (scalar(@intset) > 0){
    submitBatch(\@intset);
}

#Create new 
open(my $peppi3script,">","$outdir/PEPPI3.py");
my $peppi3=`cat $peppidir/bin/PEPPI3temp.py`;
$peppi3=~s/\!OUTDIR\!/$outdir/;
$peppi3=~s/\!PEPPIDIR\!/$peppidir/;
print $peppi3script $peppi3;
close($peppi3script);
print `chmod +x $outdir/PEPPI3.py`;

#Given a set of interactions, submit a wrapper script which runs them all
sub submitBatch{
    my @intset=@{$_[0]};
    my $args=join(",",@supported);
    $args="$args ".join(",",@intset);
    $args="$args $keepflag";
    while (`squeue -u $user | wc -l`-1 >= $maxjobs){
	print "Queue is currently full, waiting for submission...\n";
	sleep(60);
    }
    print `sbatch -J PEPPI2batch -o /dev/null -t 24:00:00 $peppidir/bin/multiwrapper.pl $args`;

    #print `$peppidir/bin/multiwrapper.pl $args`;
}

#Tree traversing algorithm which grabs a list of domains
sub treeSearch{
    my $prot=$_[0];
    my $dir=$_[1];
    
    my @domlist=();

    if (`ls $dir/$prot\_A*.fasta 2> /dev/null | wc -l` > 0){
	@domlist=treeSearch("$prot\_A",$dir);
	@domlist=(@domlist,treeSearch("$prot\_B",$dir));
    } else {
	@domlist=($prot);
    }
    
    return @domlist;
}
