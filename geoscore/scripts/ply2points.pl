#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <VisIt PLY format filename>\n" if($narg != 1);

my $filename = shift;
die "Cannot open $filename" if(!open(IN,"<$filename"));

my $on = 0;
while(<IN>)
{
    if($on > 0)
    {
	print $_;
    }
    if(/end_header/)
    {
	$on = 1;
    }
}
