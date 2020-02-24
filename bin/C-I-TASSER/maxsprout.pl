#!/usr/bin/perl

# usage: maxsprout.pl CA atom.pdb
# Please change dirs at fort.71, dglp.list

if (@ARGV != 2) {
    usage();
    exit 1;
}

$lib="/nfs/amino-library";

#############
#############
#############
`pwd`=~/(\S+)/;
$outdir=$1;
$bindir="$lib/bin/maxsprout";
$tag=int(rand()*10000000000);
$WORKDIR = "/tmp/$tag";
#############
#############
#############
`mkdir -p $WORKDIR`;
chdir "$WORKDIR";

my $INFILE = $ARGV[0];
my $OUTFILE = $ARGV[1];
my $READBRK = "$bindir/readbrk";
my $BUILDBACKBONE = "$bindir/buildbackbone";
my $DGLPLIST = "$bindir/dglp.list";
my $TORSO = "$bindir/torso2";
my $File_torso_fort71 = "$bindir/fort.71";
my $File_torso_fort72 = "$bindir/fort.72";
my $File_torso_fort81 = "$bindir/fort.81";

$INFILE = $ARGV[0];

# run MaxSprout
`cp $outdir/$INFILE $WORKDIR`;
runMaxSprout($INFILE, $OUTFILE);

########## remove crap files ################
chdir "$outdir";
open(a,"$WORKDIR/$OUTFILE");
open(b,">$outdir/$OUTFILE");
while($line=<a>){
    if($line=~/^ATOM/ && substr($line,12,4)=~/\w/){
	printf b "%54s\n",substr($line,0,54);
	printf "%54s\n",substr($line,0,54);
    }
}
close(a);
close(b);
#`mv $WORKDIR/$OUTFILE $outdir`;
`rm -fr $WORKDIR`;

### run pulchra.pl if maxsprout not work ########
$Lch1=-1;
$Lch2=-2;
$rst1=`cat $outdir/$INFILE |grep CA |wc`;
if($rst1=~/(\d+)/){
    $Lch1=$1;
}
$rst2=`cat $outdir/$OUTFILE |grep CA |wc`;
if($rst2=~/(\d+)/){
    $Lch2=$1;
}
#printf "$Lch1 -------- $Lch2\n"; 
if($Lch1 != $Lch2){
    `$lib/bin/pulchra.pl $INFILE`;
    `mv pul_$INFILE $OUTFILE`;
}

exit();


sub runMaxSprout($$) {
    my ($file_pdb, $file_pdb_mod) = @_;
    my ($PDB, @MAXSPROUT, @TORSO, $use_maxsprout, $is_empty);
    my ($file_pdb_mod_tmp, $chain_name, $file_maxsprout);
    my ($file_pdb_chain, $file_torso, $file_cmd_torso);
    my ($res, $msg);

    my ($pdb_basename) = $file_pdb;
    $pdb_basename =~ s/\..*//;

    #----------------------------------------------------------------------
    # This subroutine runs maxsprout modelling on a pdb file.
    # The structure is first split into different chains (if there is more
    # than one) and then maxsprout and torso are run for each chain. 
    # The result file is a merge of the model files for each chain.
    #----------------------------------------------------------------------

    # initialise flags
    $use_maxsprout = 0;
    
    # define tmp result file name
    $file_pdb_mod_tmp =  "$file_pdb_mod"  . ".tmp";
    
    # remove the result file if it exists
    if (-e $file_pdb_mod) { 
	unlink $file_pdb_mod;
    }
     
    # copy some mandatory files in request dir
    #open(a1,"$bindir/fort.71");
    #open(a2,">$WORKDIR/fort.71");
    #while($line=<a1>){
    #$line=~s/\.\.\/rotlib\//$bindir\/rotlib\//mg;
    #print a2 "$line";
    #}
    #close(a1);
    #close(a2);
    system("cd $WORKDIR;cp $File_torso_fort71 ./fort.71");
    system("cd $WORKDIR;cp $File_torso_fort72 ./fort.72");
    system("cd $WORKDIR;cp $File_torso_fort81 ./fort.81");

    # extract chain sub files

    $PDB = pdbExtractChain($file_pdb);

    # LOOP FOR EACH CHAIN EXTRACTED
    # =============================
    foreach $file_pdb_chain (@$PDB) {
	# build the file name
	$chain_name = $file_pdb_chain;
	$chain_name =~ s/.*\///;

        print STDERR "$0: processing chain $chain_name\n";

	# define working file names
	$file_maxsprout = $file_pdb_chain . ".maxsprout";
	$file_torso     = $file_pdb_chain . "_sc.brk";
	$file_cmd_torso = $file_pdb_chain . "_sc.cmd";
     
	# RUN MAXSPROUT PROGRAM (Now all on unix - CD - Dec 96)
	# =====================

        print STDERR "$READBRK -rd $WORKDIR -wd $WORKDIR -pdb $chain_name.brk\n";

	# run readbrk
	($res, $msg) = systemCall("$READBRK -rd $WORKDIR -wd $WORKDIR -pdb $chain_name.brk");
        die "readbrk failed: $msg\n" unless $res;

        print STDERR "echo 'y' | $BUILDBACKBONE -pl $DGLPLIST -pdb $WORKDIR/$chain_name -d -y\n";

	# run buildbackbone
	($res, $msg) = systemCall("echo 'y' | $BUILDBACKBONE -pl $DGLPLIST -pdb $WORKDIR/$chain_name -d -y ");
        die "buildbackbone failed: $msg\n" unless $res;

	# get back the model file
	($res, $msg) = systemCall("cp $WORKDIR/$chain_name.brk_mod $file_maxsprout");
        die "could not copy model file: $msg\n" unless $res;

	# clean up
	unlink $file_pdb_chain;

	# next chain if the maxsprout file is empty
        if (!containsAtom($file_maxsprout)) {
            next;
        }

	# build MAXSPROUT list
	push (@MAXSPROUT, $file_maxsprout);

	# RUN TORSO PROGRAM
	# =================

        # TODO: do proper file handling here
	# build a command file
	system "echo $file_maxsprout  >  $file_cmd_torso";
	system "echo N               >>  $file_cmd_torso";
	system "echo N               >>  $file_cmd_torso";
	system "echo N               >>  $file_cmd_torso";

        print STDERR "cd $WORKDIR; $TORSO < $file_cmd_torso\n";

	# run torso
	($res, $msg) = systemCall("cd $WORKDIR; $TORSO < $file_cmd_torso");
        die "torso failed: $msg\n" unless $res;

	# next chain if the file was not created
	# set a flag to use maxsprout file instead of torso
	unless (-e $file_torso) {
	    $use_maxsprout = 1;
	    next;
	} else {
	    # append the file to the TORSO list
	    push (@TORSO, $file_torso);
	}

    }	# END OF FOREACH CHAIN

    # BUILD THE RESULT FILE
    # =====================

    # put the header of the original pdb file
    open(PDB, "< $file_pdb") or die "could not read $file_pdb: $!\n";
    open(MOD, "> $file_pdb_mod_tmp") or die "could not write $file_pdb_mod_tmp: $!\n";
    while (<PDB>) {
	if (/^HEADER/ || /^COMPND/ || /^SOURCE/ || /^AUTHOR/) {
	    print MOD "$_";
	} else {
	    last;
	}
    }
    close PDB;
    close MOD;

    if (!$use_maxsprout) {
	# append the results of torso
	foreach $file_torso (@TORSO) {
	    if (containsAtom($file_torso)) {
		$is_empty = 0;
		system("cat $file_torso >> $file_pdb_mod_tmp");
	    }
	}
    } else {
	# append the results of maxsprout
	foreach $file_maxsprout (@MAXSPROUT) {
	    $is_empty = 0;
	    system("cat $file_maxsprout >> $file_pdb_mod_tmp");
	}
    }

    # if the result is empty: do not create it
    if ($is_empty) {
	print "DEBUG : WARNING : no model file for $file_pdb\n";
	return 0;
    }	

    # move the file to final result file
    system("mv $file_pdb_mod_tmp $file_pdb_mod");

    # remove unwanted temp files from working dir
    system("rm -f fort* *.brk_mod *.dist *.main *.maxsprout *.cmd *.side *.stat $pdb_basename" . "_*.brk");

    return 1;
}


# =================================================================
# sub to extract pdb chain files
# =================================================================
sub pdbExtractChain {
    my ($pdb_file) = @_;
    my ($basename, $chain_num, $chain, @CHAIN, $file_chain);
    my ($line, $model);

    # extract basename from $pdb_file
    $basename = $pdb_file;
    $basename =~ s/\..*//;

    print STDERR "$0: pdb_file '$pdb_file' basename '$basename'\n";

    open(PDB, "< $pdb_file") or die "could not open $pdb_file: $!\n";

    while (<PDB>) {
	chop;
	$line = $_;
	# for models (RMN), keep only the first model: exit after end of model
	if (/^MODEL / && ! $model) {
	    $model = 1;
	}
	if (/^ENDMDL / && $model) {
	    $model = 0;
	    last;
	}
	if (/^ATOM/) {
	    if (substr($_,21,1) ne $chain) {
		$chain_num++;
		$file_chain = "$basename" . "_$chain_num";
		if ($chain) {
		    close CHAIN;
		}
		open(CHAIN, "> $file_chain.brk") or die "could not open $file_chain.brk: $!";
		push(@CHAIN, $file_chain);
		$chain = substr($_,21,1);
	    }
	    print CHAIN "$line\n";
	}
    }
    if ($chain) {
	close CHAIN;
    }

    print STDERR "$0: got chain\n";

    return (\@CHAIN);
}


# =================================================================
# sub to test if a file contains line starting with ATOM 
# =================================================================
sub containsAtom{
    my ($file) = @_;

    open(ATOM, "< $file") or die "could not open $file: $!\n";

    while (<ATOM>) {
	if (/^ATOM/) {
	    close ATOM;
	    return 1;
	}
    }

    close ATOM;
    return 0;
}

# =====================================================================
#  Description: Execute a command as a system call.
#
#  Arguments:   $cmd    command (with args) to execute
#
#  Returns:     1, ''  on success
#               0, msg on error
# =====================================================================

sub systemCall ($)
{
    my $cmd = shift;

    my $rc = 0xffff & system($cmd);
    my $msg;

    if ($rc == 0) {
        #$self->logDebug(IDENT . "system: $cmd exited normally");
    }
    elsif ($rc == 0xff00) {
        $msg = "$cmd: $ERRNO";
    } 
    elsif (($rc & 0xff) == 0) {
        $rc >>= 8;
        $msg = "$cmd: exit status $rc";
    }  
    else {
        $msg = "$cmd: ";
        if ($rc &   0x80) {
            $rc &= ~0x80;
            $msg .= "coredump from ";
        } 
        $msg .= "signal $rc";
    }
    return ($rc == 0, $msg);
}

sub usage ()
{
    print STDERR "usage: $PROGRAM_NAME pdbfile outfile\n";
}

