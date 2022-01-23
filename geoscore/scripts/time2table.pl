#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <file> <opt: nper - default: 1>\n" if($narg < 1 || $narg > 2);

my $filename = shift;
die "cannot open $filename\n" if(!open(IN,"<$filename"));
my $nper = 1;
if($narg > 1)
{
    $nper = shift;
}

my $tlast = -1;
my @ts = ();
my @dts = ();
while(<IN>)
{
    my $t = 1.0 * $_;
    if($tlast >= 0)
    {
	my $dt = $t - $tlast;
	$dt /= $nper;
	push(@dts, $dt);
    }
    push(@ts, $t);
    $tlast = $t;
}
push(@dts, $tlast);
close IN;

print "<Table1D name=\"ttable\" coord=\"";
for(my $i = 0; $i < $#ts; $i++)
{
    print "$ts[$i], ";
}
print "$ts[$#ts]\" value=\"";
for(my $i = 0; $i < $#dts; $i++)
{
    print "$dts[$i], ";
}
print "$dts[$#dts]\"/>\n";
