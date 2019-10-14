use strict;
use warnings;

my $benchmarkdir="/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/ZiSet";
my $indir="$benchmarkdir/$ARGV[0]";

open(my $list,"<","$benchmarkdir/target_list");
my @springdirs=glob("$benchmarkdir/*/SPRING");
if (scalar(@springdirs) > 0){
    print "SPRING results are still remaining, save and remove these results first.\n";
    exit(1);
}

while (my $target=<$list>){
    chomp($target);
    my @chains=split("-",$target);
    next if (! -e "$indir/$target");
    print `cp -r $indir/$target $benchmarkdir/$target/SPRING`;
}
close($list);
