#!/bin/perl
use File::Basename;
use strict;

my $filename = shift;
my @arr = fileparse($filename);
my $dir = $arr[1];

my $on = 1;
my $written = 0;
if($dir=~m/ScriptInc/)
{
    $on = 0;
}

die "cannot open $filename" if(!open(IN,"<$filename"));

#FIRST, READ IN THE NEW LICENSE
my $licfile = "bsd_notice.txt";
my $heading = "";
{
  die "cannot open $licfile" if(!open(LIC,"<$licfile"));
  while(<LIC>)
  {
    $heading .= $_;
  }
  close LIC;
}

#NEXT, PRINT THE NEW LICENSE INSTEAD OF THE OLD
while(<IN>)
{
  if($on == 1)
  {
    if(/^\/\//)
    {
        if($written == 0)
        {
          print $heading;
          $written = 1;
        }
    }
    else
    {
      $on = 0;
    }
  }
  if($on == 0)
  {
    print $_;
  }
} 
close IN;
