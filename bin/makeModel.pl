#!/usr/bin/perl
#SBATCH -t 24:00:00
#SBATCH --mem=5G
#SBATCH -J makeModel.pl

use Getopt::Long qw(GetOptions);
use List::Util qw(sum0);

print `date`;
$ENV{'PATH'}="/nfs/amino-home/zhanglabs/bin:$ENV{'PATH'}";

$user="$ENV{USER}"; # user name, please change it to your own name, i.e. 'jsmith'
$outdir="";
######### Needed changes ended #################################

my $s="";
my $bindir="/nfs/amino-home/ewbell/PEPPI/bin";
my $benchflag=0;
my $domaindiv=0;
my $springdb="/nfs/amino-home/ewbell/SPRINGDB";
my $uniprotdb="/nfs/amino-library/local/hhsuite/uniprot20_2016_02/uniprot20_2016_02";
my $restriplet="/nfs/amino-library/contact/NEW/ResTriplet3_04132020";
my $trrosetta="/nfs/amino-home/zhng/local_library/anaconda3-tf/bin/python /nfs/amino-library/contact/NEW/trRosetta/trRosetta/trRosetta.py";
my $zthresh=10.0;
my $homologthresh=0.3;

GetOptions(
    "benchmark" => \$benchflag,
    "domains" => \$domaindiv,
    "outdir=s" => \$outdir,
    "target=s" => \$s
    ) or die "Invalid arguments were passed into C-I-TASSER";

my $randomTag=int(rand(1000000));
my $tempdir="/tmp/$user/PEPPI_MODEL_$s\_$randomTag";
if (! -e "$tempdir"){
    print `mkdir -p $tempdir`;
} else {
    print `rm -rf $tempdir/*`;
}
chdir($tempdir);

print `cp $outdir/$s/seq.fasta $tempdir/$s.fasta`;

if ($domaindiv){
    print "Domain division\n";
}

my $modeldir="$outdir/../model";

my $template="";
if (-f "$outdir/../hhr/$s.hhr"){
    print `cp $outdir/../hhr/$s.hhr $tempdir/$s.hhr`;
    $template=detectTemplate($s,$benchflag,$zthresh);
} else {
    print "HHsearch threading did not complete for protein $s, using trRosetta modeling\n";
}
print "$s\n";
if ($template ne ""){
    homologyModel($s,$template);
} else {
    modelSequence($s);
}

print `sync`;
print `rm -rf $tempdir`;
print `date`;

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

sub modelSequence{
    my $prot=$_[0];
    print "Making MSA for $prot...\n";
    print `$bindir/hhsuite/bin/hhblits -i $prot.fasta -oa3m $prot.a3m -d $uniprotdb -n 2 -e 0.001`;
    print "Predicting maps...\n";
    print `$restriplet/creat_npz.py $prot.a3m $prot`;
    print "Model generating...\n";
    print `$trrosetta $prot\_20.npz $prot.fasta final_1.pdb >> buildtrr.log`;
    print `mv final_1.pdb $prot.pdb`;
    print `cp $tempdir/$prot.pdb $modeldir/$prot.pdb`;
    print "Done!\n";
}
