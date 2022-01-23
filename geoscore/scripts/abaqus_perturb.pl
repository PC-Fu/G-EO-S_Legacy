#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <file> <x> <y> <z> <nx> <ny> <nz> <d> <rmin> <max perturb>\n" if ($narg != 10);

my $fname = shift;
die "Cannot open $fname\n" if(!open(IN,"<$fname"));

my @x = (0,0,0);
my @nx = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
    $x[$i] = shift;
}
for(my $i = 0; $i < 3; $i++)
{
    $nx[$i] = shift;
}
my $d = shift;
my $rmin = shift;
my $dx = shift;

#READ THROUGH ABAQUS FILE
{
    my $node_on = 0;
    my $icurr = 0;
    while(<IN>)
    {
	my $lineall = $_;
	if(/^\*(\S+)/)
	{
	    if(/NODE/)
	    {
		$node_on = 1;
	    }
	    else
	    {
		$node_on=0;
	    }
	    print $lineall;
	}
	else
	{
	    if($node_on)
	    {
		my $line = $lineall;
		$line=~s/^\s+//;
		$line=~s/\s+$//;
		my @arr = split(/\,\s*/,$line);
		my @xx = (0,0,0);
		my $dot = 0;
		for(my $i = 0; $i <= $#x; $i++)
		{
		    $xx[$i] = $arr[$i + 1];
		    $xx[$i] -= $x[$i];
		    $dot += $nx[$i] * $xx[$i];
		}
		my $d2 = sqrt($dot * $dot);
		if($d2 > $d)
		{
		    print $lineall;
		}
		else 
		{
		    my $r = 0;
		    my @xr = ();
		    for(my $i = 0; $i <= $#x; $i++)
		    {
			$xr[$i] = $xx[$i] - ($dot * $nx[$i]);
			$r += $xr[$i] * $xr[$i];
		    }
		    $r = sqrt($r);
		    #die "X, NX, D2<D, R>RMIN XX XR: $x[0] $x[1] $x[2] $nx[0] $nx[1] $nx[2] $d2 $d $r $rmin $xx[0] $xx[1] $xx[2] $xr[0] $xr[1] $xr[2]\n$arr[1] $arr[2] $arr[3]\n$lineall"; 
		    if($r < $rmin)
		    {
			print $lineall;
		    }
		    else
		    {
			print "$arr[0]";
			my $v = 0;
			for(my $i = 0; $i <= $#x; $i++)
			{
			    $v = $dx;
			    $v *= 0.5 - rand();
			    my $xcurr = $xr[$i];
			    $xcurr *= ($r + $v) / $r;
			    $xcurr += $dot * $nx[$i];
			    $xcurr += $x[$i];
			    print ", $xcurr";
			}
			print "\n";
		    }
		}
	    }
	    else
	    {
		print $lineall;
	    }
	}
    }
    close IN;
}
