use strict;
use warnings;
use List::Util qw(min);

my $benchdir="/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/ZiSet";
open(my $targetlist,"<","$benchdir/realSPRING/target_list");
open(my $fixlist,">","$benchdir/fix_list");
while (my $target=<$targetlist>){
    chomp($target);
    if (! -e "$benchdir/realSPRING/$target/SPRING/TemplateSummary.txt" || ! -e "$benchdir/$target/SPRING/TemplateSummary.txt"){
	print "TemplateSummary does not exist for target $target\n";
	print $fixlist "$target\n";
	next;
    }
    open(my $realfile,"<","$benchdir/realSPRING/$target/SPRING/TemplateSummary.txt");
    my $throwaway=<$realfile>;
    open(my $clonefile,"<","$benchdir/$target/SPRING/TemplateSummary.txt");
    for my $i (0..4){
	my $realline=<$realfile>;
	my $cloneline=<$clonefile>;
	last if (!defined($realline) || !defined($cloneline) || $realline eq "DONE");
	my @realparts=split(' ',$realline);
	my @cloneparts=split(' ',$cloneline);
	my @chains=split('-',$realparts[1]);
	if ($chains[0] ne $cloneparts[0] || $chains[1] ne $cloneparts[1] || abs(min($realparts[8],$realparts[13])-$cloneparts[3]) > 0.1 || abs(min($realparts[7],$realparts[12])-$cloneparts[4]) > 0.0011 || abs($realparts[4]-$cloneparts[5]) > 0.01){
	    print $fixlist "$target\n";
	    print "Discrepancy found in $target: ";
	    print "Chains unequal " if ($chains[0] ne $cloneparts[0] || $chains[1] ne $cloneparts[1]);
	    print "Z-score " if (abs(min($realparts[8],$realparts[13])-$cloneparts[3]) > 0.1);
	    print "TM-score " if (abs(min($realparts[7],$realparts[12])-$cloneparts[4]) > 0.0011);
	    print "Dfire " if (abs($realparts[4]-$cloneparts[5]) > 0.01);
	    print "\n";
	    last;
	}
    }
    close($realfile);
    close($clonefile);
}
