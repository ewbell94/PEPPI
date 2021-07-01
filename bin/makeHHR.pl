#!/usr/bin/perl
#SBATCH -t 24:00:00
#SBATCH --mem=5G
#SBATCH -J makeHHR.pl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw(sum0);
use POSIX qw(ceil);

$ENV{'PATH'}="/nfs/amino-home/zhanglabs/bin:$ENV{'PATH'}";

my $user="$ENV{USER}";
my $outdir="";
my $target="";
my $peppidir="/nfs/amino-home/ewbell/PEPPI";
my $domaindiv=0;
my $benchflag=0;
my $maxjobs=300;

GetOptions(
    "benchmark" => \$benchflag,
    "domains" => \$domaindiv,
    "outdir=s" => \$outdir,
    "target=s" => \$target,
    "peppidir=s" => \$peppidir
    ) or die "Invalid arguments were passed into makeHHR\n";

#User-set parameters
my $bindir="$peppidir/bin";
my $uniprotdb="/nfs/amino-library/local/hhsuite/uniprot20_2016_02/uniprot20_2016_02"; #location of Uniprot database for HHblits search
my $springdb="/nfs/amino-home/ewbell/SPRINGDB/";
my $dimerdb="$springdb/70negpos.db";
(my $sourcefasta=$target)=~s/\_[AB]//g;
my $hhdir="$outdir/$sourcefasta";
my $zthresh=8.5;
my $homologthresh=0.5;
my $hhsuitedir="$bindir/../lib/hhsuite";

print "benchmark: $benchflag\n";
#DO NOT CHANGE BENEATH THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING
$ENV{'HHLIB'}="$hhsuitedir"; #necessary for proper function of HHsearch

my $randomTag=int(rand(1000000)); #This is to prevent multiple instances from deleting eachother's directories
print "Target: $target\n";
my $tempdir="/tmp/$user/makeHHR\_$target\_$randomTag";
if (! -e "$tempdir"){
    print `mkdir -p $tempdir`;
} else {
    print `rm -rf $tempdir/*`;
}
chdir($tempdir);

print `cp $outdir/$sourcefasta/$target.fasta $tempdir/$target.fasta`;
#makeHHR threads a query sequence through DIMERDB using HHsearch
if (! -e "$hhdir/$target.hhr.gz"){
    print `$hhsuitedir/bin/hhblits -i $tempdir/$target.fasta -oa3m $tempdir/$target.a3m -d $uniprotdb -n 2 -e 0.001`;
    print `$hhsuitedir/scripts/addss.pl $tempdir/$target.a3m`;
    print `$hhsuitedir/bin/hhmake -i $tempdir/$target.a3m -id 90 -diff 100 -cov 0 -qid 0`;
    print `$hhsuitedir/bin/hhsearch -i $tempdir/$target.hhm -d $dimerdb -id 90 -diff 100 -cov 0 -qid 0 -e 0.001 -p 20 -E 0.01 -Z 30000 -z 20000 -B 30000 -b 20000`;
    print `cp $tempdir/$target.hhr $hhdir/$target.hhr`;
    print `gzip $hhdir/$target.hhr`;
}

if (! -e "$hhdir/$target.hhr.gz"){
    print "HHsearch threading was not successful\n";
    exit(1);
}

print `cp $hhdir/$target.hhr.gz $tempdir/$target.hhr.gz`;
print `gzip -f -d $tempdir/$target.hhr.gz`;
if ($domaindiv){
    my $domainbound=detectDomains($target,$benchflag,$zthresh);
    if ($domainbound > 0){
	print "Domain split at $domainbound\n";
	splitDomains($target,$domainbound,$benchflag,$outdir);
    } else {
	print "Protein is single domain!\n";
    }
}

print `sync`;
print `rm -rf $tempdir`;
print `date`;

sub splitDomains{
    my $target=$_[0];
    my $boundary=$_[1];
    my $benchflag=$_[2];
    my $outdir=$_[3];

    open(my $unsplitf,"<","$tempdir/$target.fasta");
    my $throwaway=<$unsplitf>;
    my $sequence="";
    while (my $line=<$unsplitf>){
	chomp($line);
	$sequence="$sequence$line";
    }

    open(my $fastA,">","$outdir/$sourcefasta/$target\_A.fasta");
    print $fastA ">$target\_A\n";
    print $fastA substr($sequence,0,$boundary-1);
    close($fastA);

    open(my $fastB,">","$outdir/$sourcefasta/$target\_B.fasta");
    print $fastB ">$target\_B\n";
    print $fastB substr($sequence,$boundary-1);
    close($fastB);

    while (`squeue -u $user | wc -l`-1 >= $maxjobs){
        sleep(60);
    }
    my $args="-o $outdir -d";
    $args="$args --benchmark" if ($benchflag);

    print `sbatch -o $outdir/$sourcefasta/out_makeHHR$target\_A.log $bindir/makeHHR.pl -t $target\_A $args`;
    print `sbatch -o $outdir/$sourcefasta/out_makeHHR$target\_B.log $bindir/makeHHR.pl -t $target\_B $args`;
    print `sync`;
    print `rm -rf $tempdir`;
    print `date`;
    chdir("~");
    exit(0);
}

sub getSeqID{
    my $fname1=$_[0];
    my $fname2=$_[1];
    return 0.0 if (! -f $fname1 || ! -f $fname2);
    my $NWresult;
    if ($fname2=~/\.fasta/){
        $NWresult=`$bindir/NWalign $fname1 $fname2`;
    } elsif ($fname2=~/\.pdb/){
        $NWresult=`$bindir/NWalign $fname1 $fname2 2`;
    } else {
        return 0.0;
    }
    $NWresult=~/Identical length:\s+(\d+)/;
    my $idcount=$1;
    $NWresult=~/Length of sequence 1:\s+(\d+).*\nLength of sequence 2:\s+(\d+)/;
    my $seq1len=$1;
    my $seq2len=$2;
    $idcount=1 if ($idcount==0);
    return min($idcount/$seq1len,$idcount/$seq2len) if ($fname2=~/\.fasta/);
    return $idcount/$seq1len if ($fname2=~/\.pdb/);
    return 0.0;
}

sub detectDomains{
    my $prot=$_[0];
    my $benchmark=$_[1];
    my $zthresh=$_[2];

    open(my $hhrfile,"<","$tempdir/$prot.hhr");
    my %templates=();
    my @scores=();
    while (my $line=<$hhrfile>){
        next if (!($line=~/^>/));
        chomp($line);
        $line=~/^>(.*)/;
        my $tempname=$1;

        my $scoreline=<$hhrfile>;
	$scoreline=~/Score=(\d+\.\d+)/;
        my $sumprobs=$1;
	push(@scores,$sumprobs);
        $templates{$tempname}=$sumprobs if (! exists($templates{$tempname}) || $sumprobs > $templates{$tempname});
    }
    seek($hhrfile,0,0);

    my $mu=sum0(@scores)/scalar(@scores);
    my $sd=0.0;
    for my $score (@scores){
        $sd +=(($score-$mu)**2)/(scalar(@scores)-1);
    }
    $sd=$sd**0.5;

    for my $key (keys(%templates)){
        my $z=($templates{$key}-$mu)/$sd;
        $templates{$key}=$z;
    }
    
    my $Lgap=80;
    my $nstrong=0;
    my $ngap=0;
    my $nfull=0;
    my @ngaps=();
    my @cgaps=();

    #It's hacky, I know
    my $matchline=<$hhrfile>;
    $matchline=<$hhrfile>;
    $matchline=~/ (\d+)/;
    my $Lch=$1;
    for my $i (0..6){
	$matchline=<$hhrfile>;
    }

    my @usedTemplates=();
    while (my $line=<$hhrfile>){
	last if (!($line=~/^\s*\d+\s+/));
	my @parts=split(" ",$line);
	my $z=$templates{$parts[1]};
	next if ($z < $zthresh);
	next if ($benchmark && getSeqID("$prot.fasta","$springdb/monomers/$parts[1].pdb") > $homologthresh);

	my $used=0;
	for my $t (@usedTemplates){
	    if ($t eq $parts[1]){
		$used=1;
		last;
	    }
	}
	if ($used){
	    next;
	} else {
	    my $clustline=`grep "$parts[1]" $springdb/monomers.aliases`;
	    chomp($clustline);
	    my @usedClust=split(",",$clustline);
	    push(@usedTemplates,@usedClust);
	}

	print "$parts[1]:$z\n";
	$nstrong++;
	my @bounds=split("-",$parts[8]);
	if ($bounds[0]>$Lgap || $Lch-$bounds[1]>$Lgap){
	    $ngap++;
	    if ($bounds[0]>$Lgap){
		push(@ngaps,$bounds[0]);
	    } elsif ($Lch-$bounds[1]>$Lgap){
		push(@cgaps,$bounds[1]);
	    }
	} else {
	    $nfull++;
	}
    }
    close($hhrfile);

    if ($nstrong >= 4 && $nfull < $nstrong/2){
	if (scalar(@ngaps) > 0 && scalar(@cgaps) > 0){
	    my @mids=();
	    for my $n (@ngaps){
		for my $c (@cgaps){
		    my $mid=ceil(($n+$c)/2);
		    push(@mids,$mid);
		}
	    }
	    return int(sum0(@mids)/scalar(@mids)+0.5);
	} elsif (scalar(@cgaps)>0){
	    return int(sum0(@cgaps)/scalar(@cgaps)+0.5);
	} else {
	    return int(sum0(@ngaps)/scalar(@ngaps)+0.5);
	}
    } else {
	return 0;
    }
}
 
