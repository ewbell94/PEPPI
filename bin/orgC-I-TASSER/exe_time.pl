#!/usr/bin/perl

#############################################
# usage exe_time.pl time command
#
# this program is to implement a unix command in limited time
# if finished normal, print "finished"
# if not normal, print "unfinished"
#############################################

$timeout="$ARGV[0]"*60; # in minite timeout
$exe="$ARGV[1]";
for($i=2;$i<=10;$i++){
    goto pos1 if(length $ARGV[$i] <1);
    $exe .=" $ARGV[$i]";
}
 pos1:;
printf "---$exe=====\n";
#exit();


eval{
    local $SIG{ALRM} = sub{ die "alarm\n" };
    alarm $timeout;
    $rst=qx($exe);
    printf "$rst\n";
    alarm 0;
};

if($@){
    die unless $@ eq "alarm\n";
    print "unfinished\n";
    qx(ps -ef | grep "$exe" | egrep -v " grep | vi | more | cat | pg | egrxep " | awk '{print \$2}' | xargs kill > /dev/null 2>&1);
}else{
    print "finished\n";
}

exit();
