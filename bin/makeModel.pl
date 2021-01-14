#!/usr/bin/perl
#SBATCH -t 24:00:00
#SBATCH --mem=10G
#SBATCH -J makeModel.pl

use Getopt::Long qw(GetOptions);

print `date`;
$ENV{'PATH'}="/nfs/amino-home/zhanglabs/bin:$ENV{'PATH'}";

$user="$ENV{USER}"; # user name, please change it to your own name, i.e. 'jsmith'
$outdir="";
######### Needed changes ended #################################

my $s="";
my $bindir="/nfs/amino-home/ewbell/PEPPI/bin";
my $springdb="/nfs/amino-home/ewbell/SPRINGDB/";
my $dimerdb="$springdb/70CDHITstruct.db";
my $uniprotdb="/nfs/amino-library/local/hhsuite/uniprot20_2016_02/uniprot20_2016_02";
my $restriplet="/nfs/amino-library/contact/NEW/ResTriplet3_04132020";
my $trrosetta="/nfs/amino-home/zhng/local_library/anaconda3-tf/bin/python /nfs/amino-library/contact/NEW/trRosetta/trRosetta/trRosetta.py";

GetOptions(
    "outdir=s" => \$outdir,
    "target=s" => \$s
    ) or die "Invalid arguments were passed into makeModel";

my $randomTag=int(rand(1000000));
my $tempdir="/tmp/$user/PEPPI_MODEL_$s\_$randomTag";
if (! -e "$tempdir"){
    print `mkdir -p $tempdir`;
} else {
    print `rm -rf $tempdir/*`;
}
chdir($tempdir);

(my $sourcefasta=$s)=~s/\_[AB]//g;
print `cp $outdir/$sourcefasta/$s.fasta $tempdir/$s.fasta`;

my $modeldir="$outdir/$sourcefasta";
my $toptm=1000;

modelSequence($s);
tmSearch($s);

print `sync`;
print `rm -rf $tempdir`;
print `date`;

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
        if ($benchflag && getSeqID("$tempdir/$prot.fasta","$springdb/monomers/$line.pdb") > $\
	    homologthresh){
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
        $tmscores[$i][1]=TMalign("$springdb/monomers/$tmscores[$i][0].pdb","$tempdir/$prot.pd\
b");
    }
    @tmscores=sort{$b->[1]<=>$a->[1]} @tmscores;

    open(my $tmres,">","$modeldir/$prot.tm");
    for my $i (0..scalar(@tmscores)-1){
        print $tmres "$tmscores[$i][0] $tmscores[$i][1]\n";
    }
}
