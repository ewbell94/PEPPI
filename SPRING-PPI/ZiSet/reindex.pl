use strict;
use warnings;

my $target=$ARGV[0];
my $indir="/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/ZiSet/$target/SPRING";
my $bindir="/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/bin";

exit(0) if (! -e "$indir/model0.pdb");
my @chain=split("-",$target);
open(my $newpdb,">","$indir/cpx.pdb");
open(my $oldpdb,"<","$indir/model0.pdb");

open(my $chainA,"<","$indir/$chain[0].pdb");
my $chainletter=substr($chain[0],-1);
while (my $line=<$chainA>){
    my $writeline=<$oldpdb>;
    if (substr($line,17,3) ne substr($writeline,17,3)){
	print "Residues are not equal\n";
	exit(1);
    }
    my $newindex=substr($line,22,4);
    substr($writeline,22,4)=$newindex;
    substr($writeline,21,1)=$chainletter;
    print $newpdb $writeline;
}
close($chainA);
my $throwaway=<$oldpdb>;
print $newpdb "TER\n";

open(my $chainB,"<","$indir/$chain[1].pdb");
$chainletter=substr($chain[1],-1);
while (my $line=<$chainB>){
    my $writeline=<$oldpdb>;
    if (substr($line,17,3) ne substr($writeline,17,3)){
	print "Residues are not equal\n";
	exit(1);
    }
    my $newindex=substr($line,22,4);
    substr($writeline,22,4)=$newindex;
    substr($writeline,21,1)=$chainletter;
    print $newpdb $writeline;
}
close($chainB);
print $newpdb "TER\n";

close($newpdb);
close($oldpdb);
