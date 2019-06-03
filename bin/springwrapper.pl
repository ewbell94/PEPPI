#!/usr/bin/env perl

my $springdir=$ARGV[0];
my $resdir=$ARGV[1];

for my $i (2..scalar(@ARGV)-1){
    my $query=$ARGV[$i];
    my $hhargs="-hhsearch";
    my @chains=split("-",$query);
    $hhargs="-hhr1 $resdir/../hhr/$chains[0].hhr -hhr2 $resdir/../hhr/$chains[1].hhr" if (-e "$resdir/../hhr/$chains[0].hhr" && -e "$resdir/../hhr/$chains[1].hhr");
    print `$springdir/spring.py -q $query -iDir $resdir -dockMono -scut 1.1 $hhargs`;
    for my $chain (@chains){
	print `mv $resdir/$query/SPRING/$chain.hhr $resdir/../hhr/` if (-e "$resdir/$query/SPRING/$chain.hhr");
    }
}
