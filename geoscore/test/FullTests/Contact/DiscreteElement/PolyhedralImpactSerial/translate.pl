#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <filename> <nodes> <translate-x>\n" if($narg != 3);
my $fname = shift;
my $nnds = shift;
my $dx = shift;

my $node_on = 0;
my $ind = 0;
open(IN,"<$fname");
while(<IN>)
{
    if(/^\s*\*/)
    {
	if(/NODE/)
	{
	    $node_on = 1;	    
	}
	else
	{
	    $node_on = 0;
	}
	print $_;
    }
    else
    {
	if($node_on)
	{
	    my $line = $_;
	    chomp($line);
	    $line=~s/^\s+//;
	    my @arr = split(/\,\s*/,$line);
	    my $x = $arr[1];
	    if($ind < $nnds)
	    {
		$x += $dx;
	    }
	    $ind++;
	    print "$arr[0],$x,$arr[2],$arr[3]\n";
	}
	else
	{
	    print $_;
	}
    }
}
