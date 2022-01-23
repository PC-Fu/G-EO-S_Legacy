#!/bin/perl
use strict;
use Math::Trig; 

my $narg = $#ARGV + 1;

die "usage: $0 <points>\n" if($narg != 1);

my $fname = shift;
my @lastpt = ();
die "cannot open $fname" if(!open(IN,"<$fname"));
my $arclength = 0;
my $atansum = 0;
while(<IN>)
{
    my $line = $_;
    $line=~s/^\s+//;
    $line=~s/\s+$//;
    my @arr = split(/\s+/,$line);
    die "don't recognize line" if($#arr < 1);
    my @pt = ($arr[0], $arr[1]);
    if($#lastpt > -1)
    {
	my $dx = $pt[0] - $lastpt[0];
	my $dx2 = $dx * $dx;
	my $dy = $pt[1] - $lastpt[1];
	my $dy2 = $dy * $dy;
	$dx2 += $dy2;
	if($dx2 > 0)
	{
	    my $arc = sqrt($dx2);
	    $atansum += $arc * atan2($dx, $dy);
	    $arclength += $arc;
	}
    }
    $lastpt[0] = $pt[0];
    $lastpt[1] = $pt[1];
}
close IN;

die "arc length of the given path is 0!" if($arclength== 0);

$atansum *= 45 / (atan(1) * $arclength);
print "$atansum $arclength\n";
