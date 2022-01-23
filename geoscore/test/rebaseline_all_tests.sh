#!/bin/bash
for f in `find ./ -name "*.ats"`
do
    echo "PROCESSING $f"
    echo `perl -e 'use File::Basename; my $f = shift; my @arr = fileparse($f); print "pushd $arr[1]";' -f $f`
    `perl -e 'use File::Basename; my $f = shift; my @arr = fileparse($f); print "pushd $arr[1]";' -f $f`
    pwd

    perl -e 'use File::Basename;
use strict;
my $filename = shift;
my $newdir = shift;

if($filename=~m/gpac\.ats/)
{
  exit 1;
}

print "-->(1)copying $filename definition to $newdir/restart.Linux\n";

my @arr = fileparse($filename);
my $dir = $arr[1];
die "cannot open $filename" if(!open(IN,"<$arr[0]"));
my $status = 1;
while(<IN>)
{
  if(/restartNow\=\"(\S+)\"/)
  {
    my $rnow = $1;
    system("mkdir -p $newdir");
    system("cp $rnow $newdir"."/restart.Linux");
    $status = 0;
  }
}
close IN;
exit $status;
' -f $f `ls -ltr | grep restart.Linux | awk '{print $11}' | sed s/restart\.Linux//`
#    if [ $? -eq "0" ]; then 
#	popd
#	exit
#    else
	popd
#    fi
done
