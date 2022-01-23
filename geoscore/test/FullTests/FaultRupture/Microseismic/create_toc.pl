#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <min toc> <max toc> <toc range> <dx>\n" if($narg != 4);

my $tocmin = shift;
my $tocmax = shift;
my $tocrange = shift;
my $dx = shift;
die "min toc ($tocmin) should be greater than max toc ($tocmax)" if($tocmin > $tocmax);
die "min toc ($tocmin) should be greater than 0" if($tocmin <= 0);
die "max toc ($tocmax) should be less than 1" if($tocmax >= 1);
die "toc range ($tocrange) should be >= 0" if($tocrange < 0);
die "dx ($dx) should be > 0" if($dx <= 0);

my @x = ();
my @y = ();
my @z = ();

my $n = int(2.0 / $dx);
if($n < 3)
{
    $n = 3;
}
$dx = 2.0/$n;

open(OUTX, ">x");
open(OUTY, ">y");
open(OUTZ, ">z");
for(my $i = 0; $i <= $n; $i++)
{
    my $xx = -1 + $i * $dx;
    print OUTX "$xx\n";
    print OUTY "$xx\n";
    print OUTZ "$xx\n";
    push(@z, $xx);
}
close(OUTX);
close(OUTY);
close(OUTZ);

open(OUT,">toc");
my $nn = $n + 1;
print OUT "$nn $nn $nn\n";
for(my $k = 0; $k <= $n; $k++)
{
    my $zz = $z[$k];
    my $fct = (1.0 - $zz)/2.0;
    my $tocmean = $tocmin + $fct * ($tocmax - $tocmin);
    for(my $j = 0; $j <= $n; $j++)
    {
	for(my $i = 0; $i <= $n; $i++)
	{
	    my $xx = $tocmean + $tocrange * (0.5 - rand());
	    if($xx < $tocmin)
	    {
		$xx = $tocmin;
	    }
	    if($xx > $tocmax)
	    {
		$xx = $tocmax;
	    }
	    print OUT "$xx\n";
	}
    }
}
close OUT;
