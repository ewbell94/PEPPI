use strict;
use warnings;

open(my $targets,"<","target_list");
while(my $target=<$targets>){
    chomp($target);
    my @chains=split('-',$target);
    print `cat $target/$chains[0].pdb > $target/cpx.pdb`;
    print `cat $target/$chains[1].pdb >> $target/cpx.pdb`;
}
