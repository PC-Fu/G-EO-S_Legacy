#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <z values in increasing order>\n" if($narg < 1);

my $dH0dz = -25.56e3;
my $dH1dz = -13.06e3;
my $dvdz = -25.56e3;
my $dhdz = -13.18e3;
my $z0 = -1000;

my $zlast = -1e100;
for(my $i = 0; $i < $narg; $i++)
{
    my $z = shift;
    die "$z < $zlast <-- z must be in ascending order!" if($z < $zlast);
    my $sigH = $dH0dz * $z;
    if($z < $z0)
    {
	$sigH = $dH0dz * $z0 + $dH1dz * ($z - $z0);
    }
    my $sigh = $dhdz * $z;
    my $sigv = $dvdz * $z;
    print "$z $sigH $sigh $sigv\n";
    $zlast = $z;
}
