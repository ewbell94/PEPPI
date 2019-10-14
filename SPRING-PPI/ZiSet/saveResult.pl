use strict;
use warnings;

my $indir="/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/ZiSet";
my $outdir="$indir/$ARGV[0]";
if (! -e $outdir){
    print `mkdir $outdir`;
} else {
    print `rm -rf $outdir/*`;
}

open(my $list,"<","$indir/target_list");
while (my $target=<$list>){
    chomp($target);
    my @chains=split("-",$target);
    print `mkdir $outdir/$target`;
    print `cp $indir/$target/SPRING/TemplateSummary.txt $outdir/$target`;
    print `cp $indir/$target/SPRING/cpx.pdb $outdir/$target`;
    print `cp $indir/$target/SPRING/model0.pdb $outdir/$target`;
    #print `cp $indir/$target/SPRING/tms.txt $outdir/$target`;
    #print `cp $indir/$target/SPRING/pdb.tar.gz $outdir/$target`;
    print `cp $indir/$target/SPRING/$chains[0].pdb $outdir/$target`;
    print `cp $indir/$target/SPRING/$chains[1].pdb $outdir/$target`;
}
close($list);

print `cp $indir/res.csv $outdir/res.csv`;
