#!/usr/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <lower-x lower-y lower-z upper-x upper-y upper-z> <property> <value>\n" if($narg != 8);

my @x0 = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
    $x0[$i] = shift;
}

my @x1 = (0,0,0);
my @dx = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
    $x1[$i] = shift;
    $dx[$i] = ($x1[$i] - $x0[$i]) / 3.0;
}

my $pstr = shift;
my $pval = shift;


open(P,">$pstr");
print P "4 4 4\n";
for(my $k = 0; $k < 4; $k++)
{
    for(my $j = 0; $j < 4; $j++)
    {
	for(my $i = 0; $i < 4; $i++)
	{
	    print P "$pval\n";
	}
    }
}
close P;

open(X,">x");
for(my $i = 0; $i < 4; $i++)
{
    my $val = $i * $dx[0];
    print X "$val\n";
}
close X;

open(Y,">y");
for(my $i = 0; $i < 4; $i++)
{
    my $val = $i * $dx[1];
    print Y "$val\n";
}
close Y;

open(Z,">z");
for(my $i = 0; $i < 4; $i++)
{
    my $val = $i * $dx[2];
    print Z "$val\n";
}
close Z;


#print "         <Table3D name=\"$pstr\" x_file=\"x\" y_file=\"y\" z_file=\"z\" voxel_file=\"$pstr\"/>\n";
print "         <InitialConstitutiveValue object=\"Fault_ElementManager\" propertytype=\"$pstr\" tablename=\"$pstr\" toregion=\"EB1\"/>\n";
