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

modelSequence($s);

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
