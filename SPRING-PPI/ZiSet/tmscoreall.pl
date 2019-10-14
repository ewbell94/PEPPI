use strict;
use warnings;

open(my $targetlist,"<","target_list");
while (my $target=<$targetlist>){
    chomp($target);
    print "$target\n";
    next if (-e "$target/SPRING/tms.txt" || ! -e "$target/SPRING/TemplateSummary.txt");
    print `tar -zxf $target/SPRING/pdb.tar.gz`;
    print `mv ./tmp/*/*/*.pdb $target/SPRING/`;
    print `rm -rf ./tmp`;
    open(my $templatesummary,"<","$target/SPRING/TemplateSummary.txt");
    open(my $tmscorefile,">","$target/SPRING/tms.txt");
    while (my $line=<$templatesummary>){
	my @parts=split(" ",$line);
	my $template="$parts[0]-$parts[1].pdb";
	$template=~s/\//_/g;
	if (! -e "$target/SPRING/$template"){
	    print "Template file not found for $target: $template\n";
	    next;
	}
	print `mv $target/SPRING/$template $target/SPRING/model0.pdb`;
	print `perl reindex.pl $target`;
	my $tmresult=`./TMscore -c $target/SPRING/cpx.pdb $target/cpx.pdb`;
	$tmresult=~/TM-score    = ([01]\.\d+)/;
	my $tmscore=$1;
	print $tmscorefile "$parts[0]-$parts[1],$tmscore\n";
    }
}
