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

my $user="$ENV{USER}"; # user name, please change it to your own name, i.e. 'jsmith'
my $outdir="";
######### Needed changes ended #################################

my $target="";
my $bindir="/nfs/amino-home/ewbell/PEPPI/bin";
my $domaindiv=0;
my $benchflag=0;
my $maxjobs=300;

GetOptions(
    "benchmark" => \$benchflag,
    "domains" => \$domaindiv,
    "outdir=s" => \$outdir,
    "target=s" => \$target
    ) or die "Invalid arguments were passed into makeHHR\n";

#User-set parameters
my $uniprotdb="/nfs/amino-library/local/hhsuite/uniprot20_2016_02/uniprot20_2016_02"; #location of Uniprot database for HHblits search
my $springdb="/nfs/amino-home/ewbell/SPRINGDB/";
my $dimerdb="$springdb/70CDHITstruct.db";
(my $sourcefasta=$target)=~s/\_[AB]//g;
my $hhdir="$outdir/$sourcefasta";
my $modeldir=$hhdir;
my $zthresh=3.8;
my $homologthresh=0.3;
my $toptm=1000;

print "benchmark: $benchflag\n";
#DO NOT CHANGE BENEATH THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING
#Processed parameters
$ENV{'HHLIB'}="$bindir/hhsuite/"; #necessary for proper function of HHsearch

my $randomTag=int(rand(1000000)); #This is to prevent multiple instances from deleting eachother's directories
print "Target: $target\n";
my $tempdir="/tmp/$user/makeHHR\_$target\_$randomTag";
if (! -e "$tempdir"){
    print `mkdir $tempdir`;
} else {
    print `rm -rf $tempdir/*`;
}
chdir($tempdir);

print `cp $outdir/$sourcefasta/$target.fasta $tempdir/$target.fasta`;
#makeHHR threads a query sequence through DIMERDB using HHsearch
if (! -e "$hhdir/$target.hhr.gz"){
    print `$bindir/hhsuite/bin/hhblits -i $tempdir/$target.fasta -oa3m $tempdir/$target.a3m -d $uniprotdb -n 2 -e 0.001`;
    print `$bindir/hhsuite/scripts/addss.pl $tempdir/$target.a3m`;
    print `$bindir/hhsuite/bin/hhmake -i $tempdir/$target.a3m -id 90 -diff 100 -cov 0 -qid 0`;
    print `$bindir/hhsuite/bin/hhsearch -i $tempdir/$target.hhm -d $dimerdb -id 90 -diff 100 -cov 0 -qid 0 -e 0.001 -p 20 -E 0.01 -Z 30000 -z 20000 -B 30000 -b 20000`;
    print `cp $tempdir/$target.hhr $hhdir/$target.hhr`;
    print `gzip $hhdir/$target.hhr`;
}

if (! -e "$hhdir/$target.hhr.gz"){
    print "HHsearch threading was not successful, creating trRosetta model\n";
    submitMakeModel($outdir,$target);
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

my $template=detectTemplate($target,$benchflag,$zthresh);

if ($template ne ""){
    homologyModel($target,$template);
} else {
    submitMakeModel($outdir,$target);
}

if (! -e "$modeldir/$target.pdb"){
    print "Homology modeling was not successful, creating trRosetta model\n";
    submitMakeModel($outdir,$target);
}

tmSearch($target);

print `sync`;
print `rm -rf $tempdir`;
print `date`;

sub submitMakeModel{
    my $outdir=$_[0];
    my $target=$_[1];

    my $args="-o $outdir -t $target\n";

    print `sbatch -o $outdir/$sourcefasta/out_makeModel_$target.log $bindir/makeModel.pl $args`;
    print `sync`;
    print `rm -rf $tempdir`;
    print `date`;
    chdir("~");
    exit(0);
}

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
    my $args="-o $outdir -d\n";
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

sub detectTemplate{
    my $prot=$_[0];
    my $benchmark=$_[1];
    my $zthresh=$_[2];

    open(my $hhrfile,"<","$tempdir/$prot.hhr");
    my @templates=();
    my @scores=();
    while (my $line=<$hhrfile>){
        next if (!($line=~/^>/));
        chomp($line);
        $line=~/^>(.*)/;
        my $tempname=$1;

        my $scoreline=<$hhrfile>;
        $scoreline=~/Sum_probs=(.*)/;
        my $sumprobs=$1;
        push(@scores,$sumprobs);
        my @pair=($tempname,$sumprobs);
        #print "$tempname,$sumprobs\n";
        push(@templates,\@pair);
    }
    close($hhrfile);

    my $mu=sum0(@scores)/scalar(@scores);
    my $sd=0.0;
    for my $score (@scores){
        $sd +=(($score-$mu)**2)/(scalar(@scores)-1);
    }
    $sd=$sd**0.5;

    my $returntemp="";

    @templates=sort{$b->[1]<=>$a->[1]} @templates;
    for my $i (0..scalar(@templates)-1){
        my $z=($templates[$i][1]-$mu)/$sd;
        last if ($z < $zthresh);
        next if ($benchmark && getSeqID("$prot.fasta","$springdb/monomers/$templates[$i][0].pdb") > $homologthresh);
        $returntemp=$templates[$i][0];
        print "Using template $templates[$i][0] with Z-score of $z for modeling\n";
        last;
    }

    return $returntemp;
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
        $scoreline=~/Sum_probs=(.*)/;
        my $sumprobs=$1;
	push(@scores,$sumprobs);
        $templates{$tempname}=$sumprobs
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

    #It's jank, I know
    my $matchline=<$hhrfile>;
    $matchline=<$hhrfile>;
    $matchline=~/ (\d+)/;
    my $Lch=$1;
    for my $i (0..6){
	$matchline=<$hhrfile>;
    }

    while (my $line=<$hhrfile>){
	last if (!($line=~/^\s*\d+\s+/));
	my @parts=split(" ",$line);
	my $z=$templates{$parts[1]};
	next if ($z < $zthresh || $parts[3] > 0.01);
	next if ($benchmark && getSeqID("$prot.fasta","$springdb/monomers/$parts[1].pdb") > $homologthresh);
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

    if ($nfull < $nstrong/2){
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

sub homologyModel{
    my $query=$_[0];
    my $template=$_[1];

    my %onetothree=('A'=>"ALA",'C'=>"CYS",'D'=>"ASP",'E'=>"GLU",'F'=>"PHE",
                    'G'=>"GLY",'H'=>"HIS",'I'=>"ILE",'K'=>"LYS",'L'=>"LEU",
                    'M'=>"MET",'N'=>"ASN",'P'=>"PRO",'Q'=>"GLN",'R'=>"ARG",
                    'S'=>"SER",'T'=>"THR",'V'=>"VAL",'W'=>"TRP",'Y'=>"TYR",
                    'B'=>"BBB",'Z'=>"ZZZ",'X'=>"XYZ",'U'=>"SEC",'O'=>"PYL");

    my @alignment=();
    my @qaa=();
    open(my $hhresultfile,"<","$tempdir/$query.hhr");
    my $readflag=0;
    while (my $line=<$hhresultfile>){
        if ($line=~/>$template/){
            while(1){
                for my $i (0..3){
                    $line=<$hhresultfile>;
                    last if ($line=~/Done!/);
                }
                last if (!($line=~/^Q/));
                #print "Query line:\n";
                #print $line;
                $line=~/Q .*\s(\d+) (\S+)\s+\d+ \(/;
                my $startq=$1;
                my $qseq=$2;
                #print "$startq,$qseq,test\n";
                for my $i (0..3){
                    $line=<$hhresultfile>;
                }
                #print "Template line:\n";
                #print $line;
                $line=~/T .*\s(\d+) (\S+)\s+\d+ \(/;
                my $startt=$1;
                my $tseq=$2;
            #print "$startt,$tseq,test\n";
                while (scalar(@qaa)<$startq){
                    push(@qaa,"XYZ");
                }
                while (scalar(@alignment)<$startt){
                    push(@alignment,-1);
                }
                print "Error: different sequnce lengths\n" if (length($qseq) != length($tseq));
                for my $i (0..length($qseq)-1){
                    my $qchar=substr($qseq,$i,1);
                    my $tchar=substr($tseq,$i,1);
                    if ($qchar eq "-"){
                        push(@alignment,-1);
                    } elsif ($tchar eq "-") {
                        if ($onetothree{$qchar} ne ""){
                            push(@qaa,$onetothree{$qchar});
                        } else {
                            push(@qaa,"UNK");
                        }
                        $startq++;
                    } else {
                        if ($onetothree{$qchar} ne ""){
                            push(@qaa,$onetothree{$qchar});
                        } else {
                            push(@qaa,"UNK");
                        }
                        push(@alignment,$startq);
                        $startq++;
                    }
                }
                for my $i (0..1){
                    $line=<$hhresultfile>;
                }
            }
            last;
        }
	
    }
    close($hhresultfile);
    
    open(my $modelout,">","$tempdir/$query.pdb");
    open(my $tempin,"<","$springdb/monomers/$template.pdb");
    my $i=1;
    while (my $line=<$tempin>){
        next if (substr($line,0,4) ne "ATOM" || substr($line,12,4) ne " CA ");
        my $resnum=substr($line,22,4);
        next if ($resnum >= scalar(@alignment) || $alignment[$resnum] < 0);
        my $resname=substr($line,17,3);
        chomp($line);
        substr($line,17,3)=$qaa[$alignment[$resnum]];
        substr($line,6,5)=sprintf("%5s",$resnum);
        substr($line,22,4)=sprintf("%4s",$alignment[$resnum]);
        $line=$line.sprintf("%5s",$resnum).sprintf(" %s",$resname);
        print $modelout "$line\n";
        $i++;
    }
    for my $j (1..$i-1){
        my $connection=sprintf("CONECT%5s%5s\n",$j,$j+1);
        #print $modelout $connection;
    }
    close($modelout);
    close($tempin);
    print `cp $tempdir/$query.pdb $modeldir/$query.pdb`;
}

sub TMalign{
    my $template=$_[0];
    my $query=$_[1];
    my $tmres=`$bindir/TMalign $template $query`;
    $tmres=~/TM-score= (.*) \(if normalized by length of Chain_1.*\nTM-score= (.*) \(if normalized by length of Chain_2.*/;
    my $tm1=$1;
    my $tm2=$2;
    #print "$tm1,$tm2\n";
    return ($tm1+$tm2)/2;
}

sub fTMalign{
    my $template=$_[0];
    my $query=$_[1];
    my $tmres=`$bindir/fTMalign $template $query -fast`;
    $tmres=~/TM-score= (.*) \(if normalized by length of Chain_1.*\nTM-score= (.*) \(if normalized by length of Chain_2.*/;
    my $tm1=$1;
    my $tm2=$2;
    #print "$tm1,$tm2\n";
    return ($tm1+$tm2)/2;
}

sub tmSearch{
    my $prot=$_[0];

    (my $monolist=$dimerdb)=~s/\.db/\.mono/;
    open(my $monofile,"<",$monolist);
    my @tmscores=();
    while (my $line=<$monofile>){
	chomp($line);
	my @pair;
	if ($benchflag && getSeqID("$tempdir/$prot.fasta","$springdb/monomers/$line.pdb") > $homologthresh){
	    print "Template $line given low TM-score because of high homology.\n";
	    @pair=($line,0.00001);
	} else {
	    @pair=($line,fTMalign("$springdb/monomers/$line.pdb","$tempdir/$prot.pdb"));
	}
	push(@tmscores,\@pair);
    }
    close($monofile);

    @tmscores=sort{$b->[1]<=>$a->[1]} @tmscores;
    for my $i (0..$toptm-1){
	$tmscores[$i][1]=TMalign("$springdb/monomers/$tmscores[$i][0].pdb","$tempdir/$prot.pdb");
    }
    @tmscores=sort{$b->[1]<=>$a->[1]} @tmscores;
    
    open(my $tmres,">","$modeldir/$prot.tm");
    for my $i (0..scalar(@tmscores)-1){
	print $tmres "$tmscores[$i][0] $tmscores[$i][1]\n";
    }
}
 
