#!/usr/bin/env perl

use strict;
use warnings;

my $maxjobs=300;
my $benchmarkdir="/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/ZiSet";
open(my $targets,"<","$benchmarkdir/target_list");
#open(my $targets,"<","$benchmarkdir/shortlist");
while (my $target=<$targets>){
    print "$target";
    chomp($target);
    my $targetdir="$benchmarkdir/$target";
    next if (-e "$targetdir/SPRING/TemplateSummary.txt");
    my $templatetext=`cat "$benchmarkdir/template.pl"`;
    $templatetext=~s/OUTDIR/$targetdir\/SPRING/;
    open(my $scriptfile,">","$targetdir/$target.pl");
    print $scriptfile $templatetext;
    close($scriptfile);
    print `chmod +x $targetdir/$target.pl`;
    my @chains=split("-",$target);
    while(`qstat -u ewbell | wc -l`-5>=$maxjobs){
	sleep(60);
    }
    print `qsub -N "SPRING_$target" -o $targetdir/out.log -e $targetdir/err.log -l nodes=1:ppn=1 -l walltime=10:00:00 $benchmarkdir/$target/$target.pl -F "$targetdir/$chains[0].seq $targetdir/$chains[1].seq"`;
    #print `$benchmarkdir/$target/$target.pl $targetdir/$chains[0].seq $targetdir/$chains[1].seq`;
}
