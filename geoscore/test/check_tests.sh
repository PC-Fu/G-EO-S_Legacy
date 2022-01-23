#!/bin/bash
for f in `find ./ -name "*.ats"`
do
    perl -e 'use File::Basename;
use strict;
my $filename = shift;
my @arr = fileparse($filename);
my $dir = $arr[1];
open(IN,"<$filename");
my $hasCurrent = 0;
my $hasRestart = 0;
my $isSource = 0;
while(<IN>)
{
  if(/^source/)
  {
    if(/gpac\.ats/)
    {
    }
    else
    {
      $isSource = 1;
    }
  }
  if(/restartNow/)
  {
    $hasCurrent = 1;
  }
  if(/restartTest/)
  {
    $hasRestart = 1;
  }
} 
if($isSource == 0 && ($hasCurrent == 0 || $hasRestart == 0))
{
  print "$filename: hasCurrent = $hasCurrent hasRestart = $hasRestart\n";
}
close IN;' -f $f
done
