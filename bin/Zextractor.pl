use strict;
use warnings;

my $hhrfilename=$ARGV[0];
my @hits=fetchHits($hhrfilename);
for my $hit (@hits){
    my @pair=@{$hit};
    print "$pair[0] $pair[1]\n";
}

sub fetchHits{
    my $hhrfilename=$_[0];
    my @templates=();
    my @scores=();
    open(my $hhrfile,"<",$hhrfilename);
    while (my $line=<$hhrfile>){
        if (substr($line,0,1) eq ">"){
            chomp($line);
            (my $target=$line)=~s/>//g;
            next if (grep(/^$target$/,@templates));
            push(@templates,$target);
            my $scoreline=<$hhrfile>;
            $scoreline=~/Sum_probs=(\S+)/;
            my $score=$1;
            #print "$scoreline$score\n";
            push(@scores,$score);
        }
    }
    close($hhrfile);

    print scalar(@templates)."\n";
    my $meanval=0.0;
    for my $score (@scores){
        $meanval+=$score/scalar(@scores);
    }
    print "$meanval\n";
    my $std=0.0;
    for my $score (@scores){
        $std+=($score-$meanval)**2/scalar(@scores);
    }
    $std=$std**(0.5);

    my @pairs=();
    for my $i (0..scalar(@templates)-1){
        my @pair=($templates[$i],($scores[$i]-$meanval)/$std);
        push(@pairs,\@pair);
    }

    @pairs=sort{$b->[1]<=>$a->[1]} @pairs;
    return @pairs
}
