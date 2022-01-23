#!/bin/bash

#ARGS (1) file 1 (2) file 2
if [ $# -ne 1 ]; then
   echo "usage: ARGS (1) std output from risk calculation"
   exit 1
fi 

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    if [ -d "$GPAC_SCRIPT_PATH" ];then
	echo "$0 using GPAC_SCRIPTH_PATH $GPAC_SCRIPT_PATH"
    else
	export GPAC_SCRIPT_PATH=/usr/gapps/GEOS/scripts
	if [ -d "$GPAC_SCRIPT_PATH" ];then
	    echo "$0 using GPAC_SCRIPTH_PATH $GPAC_SCRIPT_PATH"
	else
	    echo "$0 cannot find GPAC_SCRIPT_PATH $GPAC_SCRIPT_PATH"
	    exit 1
	fi
    fi
fi

#extract the hazard headings
grep ^HAZARD $1 | awk '{print $2}' | uniq > uniq_times.tmp
echo "set xlabel 'log(acc, cm/s/s)'; set ylabel 'frequency (1/yr)'; plot \\" > hazard.gp
echo "set xlabel 'nuisance(% of people)'; set ylabel 'frequency (1/yr)'; plot \\" > risk.gp

export i=0
for t in `cat uniq_times.tmp`
do
    grep ^HAZARD $1 | awk '{print $2" "$3" "$4" "$5" "$6}' | grep ^$t > hazard_$t
    perl -e 'use strict;
my $file = shift;
my $index = shift;
my $sep = ",";
if($index==0)
{
  $sep="";
}
open(IN,"<$file");
my $line = <IN>;
my @arr = split(/\s+/,$line);
print "$sep '\''$file'\'' u 3:5 w lines ti \"$arr[0]-$arr[1]\"\\\n";
close IN;' -f hazard_$t $i >> hazard.gp

    perl -e 'use strict;
my $file = shift;
my $index = shift;
my $sep = ",";
if($index==0)
{
  $sep = "";
}
open(IN,"<$file");
my $line = <IN>;
my @arr = split(/\s+/,$line);
print "$sep '\''$file'\'' u (100*\$4):5 w lines ti \"$arr[0]-$arr[1]\"\\\n";
close IN;' -f hazard_$t $i >> risk.gp

    (( i++ ))
done
