#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use List::Util qw(min);
use POSIX qw(floor);

#User-set parameters
my $bindir="/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/bin";
my $outputdir="/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/12asA-12asB";
my $dbdir="/nfs/amino-home/liuzi/lz_program/TACOS/database";
my $maxtemplates=1000000;
my $writeAll=1;
my $scut=1.1;
#my $uniprotdb="/nfs/amino-library/local/hhsuite/uniprot20_2015_06/uniprot20_2015_06";
my $uniprotdb="/nfs/amino-library/local/hhsuite/uniprot20_2016_02/uniprot20_2016_02";
my $dimerdb="/nfs/amino-library/DIMERDB/HHsearch/hhm.db";
my $zmin=2.0;

#Processed parameters
my $user=`whoami`;
chomp($user);
$ENV{'HHLIB'}="$bindir/hhsuite/";
my @weights=(12.0,1.4);
my $homothresh=0.9;
my $seqmax=1500;
if (scalar(@ARGV) < 2){
    print "Not enough arguments were supplied\n";
    exit(1);
}
my $currdir=`pwd`;
chomp($currdir);
my $prot1file=$ARGV[0];
$prot1file="$currdir/$prot1file" if (substr($prot1file,0,1) ne "/");
my $prot2file=$ARGV[1];
$prot2file="$currdir/$prot2file" if (substr($prot2file,0,1) ne "/");
(my $preprot1=$prot1file)=~s/.*\///g;
(my $preprot2=$prot2file)=~s/.*\///g;
my @parts1=split('\.',$preprot1);
my @parts2=split('\.',$preprot2);
my $prot1=$parts1[0];
my $prot1ext=$parts1[1];
my $prot2=$parts2[0];
my $prot2ext=$parts2[1];
my $tempdir="/tmp/$user/PPI_$prot1-$prot2";
my $hhdir=$tempdir;

if (! -e "$prot1file" || ! -e "$prot2file"){
    print "Protein sequence files were not found!\n";
    exit(2);
}

if (-e "$outputdir/SPRING/TemplateSummary.txt"){
    print "SPRING has already been run!\n";
    exit(3);
}

#Make working directory
if (! -e "$tempdir"){
    print `mkdir $tempdir`;
} else {
    #print `rm -rf $tempdir/*`;
}

chdir("$tempdir");

print `cp $prot1file $tempdir/$prot1.fasta`;
print `cp $prot2file $tempdir/$prot2.fasta`;

#Copy HHR files or run HHsearch
print "Running HHsearch\n";
if (-e "$hhdir/$prot1.hhr"){
    print `cp $hhdir/$prot1.hhr $tempdir`;
} else {
    makeHHR($prot1);
    if (! -e "$tempdir/$prot1.hhr"){
	print "HHsearch failed for $prot1.\n";
	exit(4);
    }
}
if (-e "$hhdir/$prot2.hhr"){
    print `cp $hhdir/$prot2.hhr $tempdir`;
} else {
    makeHHR($prot2);
    if (! -e "$tempdir/$prot2.hhr"){
	print "HHsearch failed for $prot2.\n";
	exit(4);
    }
}

#Fetch HHsearch hits from HHR files
print "Fetching HHsearch hits\n";
my @prot1hits=fetchHits($prot1);
my @prot2hits=fetchHits($prot2);

#for my $i (0..scalar(@prot1hits)-1){
#    print "$prot1hits[$i][0],$prot1hits[$i][1]\n";
#}
#print "\n";

#for my $i (0..scalar(@prot2hits)-1){
#    print "$prot2hits[$i][0],$prot2hits[$i][1]\n";
#}
#print "\n";

#Store dfire and complex list
open(my $complexfile,"<","$dbdir/pdb/ComplexList.txt");
my @complexlines=<$complexfile>;
chomp(@complexlines);
close($complexfile);

open(my $dfirefile,"<","$bindir/dfire.txt");
my @dfire=<$dfirefile>;
chomp(@dfire);
close($dfirefile);

#Build index match dictionary
print "Building index\n";

open(my $indexfile,"<","$dbdir/SPRING/indexAll.txt");
my @indexlines=<$indexfile>;
#chomp(@indexlines);
close($indexfile);

my @prothits=();
for my $i (0..min($maxtemplates-1,scalar(@prot1hits)-1)){
    push(@prothits,$prot1hits[$i][0]);
}
for my $i (0..min($maxtemplates-1,scalar(@prot2hits)-1)){
    push(@prothits,$prot2hits[$i][0]);
}
my %index=();
for my $hit (@prothits){
    if (!exists($index{$hit})){
	my @grephits=grep(/ $hit/,@indexlines);
	#print "$hit,".scalar(@grephits)."\n";
	$index{$hit}=\@grephits;
    }
}

#Search for dimer templates given monomeric hits
my @dimerTemplates=fetchDimers(\@prot1hits,\@prot2hits);

if (getSeqID("$tempdir/$prot1.fasta","$tempdir/$prot2.fasta") < $homothresh){
    my @flippedTemplates=fetchDimers(\@prot2hits,\@prot1hits);
    for my $i (0..scalar(@flippedTemplates)-1){
	my $dupflag=0;
	for my $j (0..scalar(@dimerTemplates)-1){
	    if ($dimerTemplates[$j][0] eq $flippedTemplates[$i][0] && $dimerTemplates[$j][1] eq $flippedTemplates[$i][1]){
		$dupflag=1;
		last;
	    }
	}
	if (!$dupflag){
	    my @flipped=($flippedTemplates[$i][1],$flippedTemplates[$i][0],$flippedTemplates[$i][2]);
	    push(@dimerTemplates,\@flipped);
	}
    }
}

#Filter out templates based on cutoff seq identity
#if ($scut<1.0){
#    my @temp=@dimerTemplates;
#    @dimerTemplates=();
#    for my $i (0..scalar(@temp)-1){
#	my $seqidA=getSeqID("$tempdir/$prot1.fasta","$dbdir/pdb/PDBall/".substr($temp[$i][0],1,2)."/$temp[$i][0].pdb");
#	next if ($seqidA >= $scut);
#	my $seqidB=getSeqID("$tempdir/$prot2.fasta","$dbdir/pdb/PDBall/".substr($temp[$i][1],1,2)."/$temp[$i][1].pdb");
#	next if ($seqidB >= $scut);
#	push(@dimerTemplates,$temp[$i]);
#    }
#}

#Create and score models from selected dimer templates
if (scalar(@dimerTemplates)==0){
    print "No dimer templates found!\n";
    exit(0);
}

constructMonomer($prot1,$prot1hits[0][0]);
constructMonomer($prot2,$prot2hits[0][0]);

print "Constructing models\n";
my $seq1len=`grep -v "^>" $tempdir/$prot1.fasta | wc -c`-`grep -v "^>" $tempdir/$prot1.fasta | wc -l`;
my $seq2len=`grep -v "^>" $tempdir/$prot2.fasta | wc -c`-`grep -v "^>" $tempdir/$prot2.fasta | wc -l`;
for my $i (0..scalar(@dimerTemplates)-1){
    print "$dimerTemplates[$i][0]-$dimerTemplates[$i][1]\n";
    my $springscore=constructModel($prot1,$prot2,$dimerTemplates[$i][0],$dimerTemplates[$i][1],$dimerTemplates[$i][2]);
    $dimerTemplates[$i][2]=$springscore;
}

@dimerTemplates=sort{$b->[2]<=>$a->[2]} @dimerTemplates;
for my $i (0..scalar(@dimerTemplates)-1){
    print "$dimerTemplates[$i][0],$dimerTemplates[$i][1],$dimerTemplates[$i][2]\n";
}

sub getSeqID{
    my $fname1=$_[0];
    my $fname2=$_[1];
    return 0.0 if (! -f $fname1 || ! -f $fname2);
    my $NWresult;
    if ($fname2=~/\.fasta/){
	$NWresult=`$bindir/NWalign $fname1 $fname2`;
    } elsif ($fname2=~/\.pdb/){
	$NWresult=`$bindir/NWalign $fname1 $fname2 2`;
    } else {
	return 0.0;
    }
    $NWresult=~/Sequence identity: (.*)\(/;
    return $1*1.0;
}

sub makeHHR{
    my $prot=$_[0];
    print `$bindir/hhsuite/bin/hhblits -i $tempdir/$prot.fasta -oa3m $tempdir/$prot.a3m -d $uniprotdb -n 2 -e 0.001`;
    print `$bindir/hhsuite/scripts/addss.pl $tempdir/$prot.a3m`;
    print `$bindir/hhsuite/bin/hhmake -i $tempdir/$prot.a3m -id 90 -diff 100 -cov 0 -qid 0`;
    print `$bindir/hhsuite/bin/hhsearch -i $tempdir/$prot.hhm -d $dimerdb -id 90 -diff 100 -cov 0 -qid 0 -e 0.001 -p 20 -E 0.01 -Z 30000 -z 20000 -B 30000 -b 20000`;
}

sub fetchHits{
    my $prot=$_[0];
    my @templates=();
    my @scores=();
    open(my $hhrfile,"<","$tempdir/$prot.hhr");
    for my $i (0..8){
	my $throwaway=<$hhrfile>;
    }
    while (my $line=<$hhrfile>){
	my @parts=split(' ',$line);
	last if (!looks_like_number($parts[0]));
	push(@templates,$parts[1]);
	#push(@scores,$parts[5]/$parts[7]);
	push(@scores,$parts[5]);
	#push(@scores,-1*log($parts[3])/log(10));
    }
    close($hhrfile);
    
    my $meanval=0.0;
    for my $score (@scores){
	$meanval+=$score/scalar(@scores);
    }
    my $std=0.0;
    for my $score (@scores){
	$std+=($score-$meanval)**2/scalar(@scores);
    }
    $std=$std**(0.5);
    my @outlist=();
    for my $i (0..scalar(@templates)-1){
	my @pair=($templates[$i],($scores[$i]-$meanval)/$std);
	#my @pair=($templates[$i],$scores[$i]);
	push(@outlist,\@pair) if ($pair[1] >= $zmin && getSeqID("$tempdir/$prot.fasta","$dbdir/pdb/chains/".substr($templates[$i],1,2)."/$templates[$i].pdb") < $scut);
    }
    @outlist=sort{$b->[1]<=>$a->[1]} @outlist;
    return @outlist;
}

sub fetchDimers{
    print "Fetching dimers\n";
    my @prot1list=@{$_[0]};
    my @prot2list=@{$_[1]};
    my @dimerlist=();
    for my $i (0..min($maxtemplates-1,scalar(@prot2hits)-1)){
	my @prot1hits=@{$index{$prot1list[$i][0]}};
	my $hitcount=scalar(@prot1hits);
	for my $prot1hitline (@prot1hits){
	    my @prot1parts=split(' ',$prot1hitline);
	    my $biomol1="$prot1parts[0]/$prot1parts[1]";
	    next if ($biomol1=~/_1_/);
	    for my $j (0..min($maxtemplates-1,scalar(@prot2hits)-1)){
		my $zscore=min($prot1list[$i][1],$prot2list[$j][1]);
		my @prot2hits=@{$index{$prot2list[$j][0]}};
		my $hit2count=scalar(@prot2hits);
		for my $prot2hitline (@prot2hits){
		    my @prot2parts=split(' ',$prot2hitline);
		    my $biomol2="$prot2parts[0]/$prot2parts[1]";
		    next if ($biomol2=~/_0_/ || $biomol1 ge $biomol2);
		    if (grep(/$biomol1-$biomol2/,@complexlines)){
			#print "$biomol1,$biomol2,$zscore\n";
			my @dimerpair=($biomol1,$biomol2,$zscore);
			push(@dimerlist,\@dimerpair);
		    }
		}
	    }
	}
    }
    return @dimerlist;
}

sub constructMonomer{
    my $query=$_[0];
    my $template=$_[1];
    my $templatesub=substr($template,1,2);
    my %onetothree=('A'=>"ALA",'C'=>"CYS",'D'=>"ASP",'E'=>"GLU",'F'=>"PHE",
		    'G'=>"GLY",'H'=>"HIS",'I'=>"ILE",'K'=>"LYS",'L'=>"LEU",
		    'M'=>"MET",'N'=>"ASN",'P'=>"PRO",'Q'=>"GLN",'R'=>"ARG",
		    'S'=>"SER",'T'=>"THR",'V'=>"VAL",'W'=>"TRP",'Y'=>"TYR",
		    'B'=>"BBB",'Z'=>"ZZZ",'X'=>"XYZ");
    
    my @alignment=();
    my @qaa=();
    open(my $hhresultfile,"<","$tempdir/$query.hhr");
    my $readflag=0;
    while (my $line=<$hhresultfile>){
	if ($line=~/>$template/){
	    while(1){
		for my $i (0..3){
		    $line=<$hhresultfile>;
		    last if ($line=~/Done!/);
		}
		last if (!($line=~/^Q/));
		#print "Query line:\n";
		#print $line;
		$line=~/Q .*\s(\d+) (\S+)\s+\d+ \(/;
		my $startq=$1;
		my $qseq=$2;
		#print "$startq,$qseq,test\n";
		for my $i (0..3){
		    $line=<$hhresultfile>;
		}
		#print "Template line:\n";
		#print $line;
		$line=~/T .*\s(\d+) (\S+)\s+\d+ \(/;
		my $startt=$1;
		my $tseq=$2;
		#print "$startt,$tseq,test\n";
		while (scalar(@qaa)<$startq){
		    push(@qaa,"XYZ");
		}
		while (scalar(@alignment)<$startt){
		    push(@alignment,-1);
		}
		print "Error: different sequnce lengths\n" if (length($qseq) != length($tseq));
		for my $i (0..length($qseq)-1){
		    my $qchar=substr($qseq,$i,1);
		    my $tchar=substr($tseq,$i,1);
		    if ($qchar eq "-"){
			push(@alignment,-1);
		    } elsif ($tchar eq "-") {
			push(@qaa,$onetothree{$qchar});
			$startq++;
		    } else {
			push(@qaa,$onetothree{$qchar});
			push(@alignment,$startq);
			$startq++;
		    }
		}
		for my $i (0..1){
		    $line=<$hhresultfile>;
		}
	    }
	    last;
	}
	
    }
    close($hhresultfile);
    
    open(my $modelout,">","$tempdir/$query.pdb");
    open(my $tempin,"<","$dbdir/pdb/chains/$templatesub/$template.pdb");
    my $i=1;
    while (my $line=<$tempin>){
	next if (substr($line,0,4) ne "ATOM" || substr($line,12,4) ne " CA ");
	my $resnum=substr($line,22,4);
	next if ($resnum >= scalar(@alignment) || $alignment[$resnum] < 0);
	my $resname=substr($line,17,3);
	chomp($line);
	substr($line,17,3)=$qaa[$alignment[$resnum]];
	substr($line,6,5)=sprintf("%5s",$resnum);
	substr($line,22,4)=sprintf("%4s",$alignment[$resnum]);
	$line=$line.sprintf("%5s",$resnum).sprintf(" %s",$resname);
	print $modelout "$line\n";
	$i++;
    }
    for my $j (1..$i-1){
	my $connection=sprintf("CONECT%5s%5s\n",$j,$j+1);
	#print $modelout $connection;
    }
    close($modelout);
    close($tempin);
}

sub constructModel{
    my $prot1=$_[0];
    my $prot2=$_[1];
    my $dimer1temp=$_[2];
    my $dimer2temp=$_[3];
    my $zscore=$_[4];

    my $dimer1sub=substr($dimer1temp,1,2);
    my $dimer2sub=substr($dimer2temp,1,2);

    (my $dimer1name=$dimer1temp)=~s/\//_/g;
    (my $dimer2name=$dimer2temp)=~s/\//_/g;
    open(my $modelfile,">","$tempdir/$dimer1name-$dimer2name.pdb");
    
    my $TM1result=`$bindir/TMalign "$tempdir/$prot1.pdb" "$dbdir/pdb/PDBall/$dimer1sub/$dimer1temp.pdb" -L $seq1len -o $tempdir/out`;
    #print "$TM1result\n";
    $TM1result=~/TM-score= (.*) \(if scaled/;
    my $tm1score=$1;
    #print "$tm1score\n";
    my $end1ind=0;
    open(my $supfile,"<","$tempdir/out_all");
    while (my $line=<$supfile>){
	last if ($line=~/^TER/);
	if ($line=~/^ATOM/){
	    $end1ind=substr($line,6,5);
	    print $modelfile $line;
	}
    }
    print $modelfile "TER\n";
    close($supfile);

    my $TM2result=`$bindir/TMalign "$tempdir/$prot2.pdb" "$dbdir/pdb/PDBall/$dimer2sub/$dimer2temp.pdb" -L $seq2len -o $tempdir/out`;
    #print "$TM2result\n";
    $TM2result=~/TM-score= (.*) \(if scaled/;
    my $tm2score=$1;
    #print "$tm2score\n";
    my $end2ind=$end1ind;
    open($supfile,"<","$tempdir/out_all");
    while (my $line=<$supfile>){
	last if ($line=~/^TER/);
	if ($line=~/^ATOM/){
	    substr($line,21,1)="B";
	    $end2ind=$end1ind+substr($line,6,5);
	    substr($line,6,5)=sprintf("%5s",$end2ind);
	    print $modelfile $line;
	}
    }
    print $modelfile "TER\n";
    close($supfile);

    for my $i (1..$end1ind-1){
	my $connectline=sprintf("CONECT%5s%5s\n",$i,$i+1);
	print $modelfile $connectline;
    }
    for my $i ($end1ind+1..$end2ind-1){
	my $connectline=sprintf("CONECT%5s%5s\n",$i,$i+1);
	print $modelfile $connectline;
    }
    close($modelfile);

    my %aminocode=('ALA'=>0,'CYS'=>1,'ASP'=>2,'GLU'=>3,'PHE'=>4,
		   'GLY'=>5,'HIS'=>6,'ILE'=>7,'LYS'=>8,'LEU'=>9,
		   'MET'=>10,'ASN'=>11,'PRO'=>12,'GLN'=>13,'ARG'=>14,
		   'SER'=>15,'THR'=>16,'VAL'=>17,'TRP'=>18,'TYR'=>19,
		   'BBB'=>20,'ZZZ'=>20,'XYZ'=>20);
    my @Acoord=();
    my @Aseq=();
    my @Bcoord=();
    my @Bseq=();
    open($modelfile,"<","$tempdir/$dimer1name-$dimer2name.pdb");
    while (my $line=<$modelfile>){
	if ($line=~/^ATOM/){
	    my @coord=(substr($line,30,8),substr($line,38,8),substr($line,46,8));
	    if (substr($line,21,1) eq "A"){
		push(@Acoord,\@coord);
		push(@Aseq,substr($line,17,3));
	    } elsif (substr($line,21,1) eq "B"){
		push(@Bcoord,\@coord); 
		push(@Bseq,substr($line,17,3));
	    }
	}
    }
    close($modelfile);
    
    my $dfire=0.0;
    for my $i (0..scalar(@Aseq)-1){
	for my $j (0..scalar(@Bseq)-1){
	    my $dist=0.0;
	    for my $n (0..2){
		$dist+=($Acoord[$i][$n]-$Bcoord[$j][$n])**2;
	    }
	    $dist=$dist**(0.5);
	    if ($dist < 10.0){
		#print "$dist,$Aseq[$i] $i,$Bseq[$j] $j\n";
		my $index=$aminocode{$Aseq[$i]}*21*20+$aminocode{$Bseq[$j]}*20+floor($dist*2.0);
		$dfire+=$dfire[$index];
	    }
	}
    }
    
    my $tmscore=min($tm1score,$tm2score);
    my $springscore=$zscore*($weights[0]*$tmscore+$weights[1]*$dfire);
    print "$zscore,$tmscore,$dfire,$springscore\n\n";
    return $springscore;
}
