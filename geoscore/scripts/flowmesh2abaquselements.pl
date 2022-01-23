#!/usr/bin/perl
use strict;
use Math::Trig;

my $narg = $#ARGV + 1;
die "usage: $0 <fault definition file: space-delimited tuples of [x,y,z,length,depth,dip(deg),strike(deg)]>\nNote: depth is the (positive) dimension of the fault along the dip direction\n" if($narg != 1);

my $ffile = shift;

#---------------------------
#  GET FAULTS
#---------------------------
die "unable to open $ffile" if(!open(IN,"<$ffile"));
my $nfaults = 0;
my $inode = 1;
{
	print "*HEADING
flowmesh2abaquselements.pl script from GPAC: SJ 11/19/2012
**
*NODE, NSET=ALLNODES
";
	my $pi = 4 * atan(1);
    while(<IN>)
    {
		my $line = $_;
		$line =~s/^\s+//;
		chomp($line);
		my @arr = split(/\s+/,$line);
		
		die "fault entry must be length 7: x y z length depth dip strike!!" if($#arr != 6);
		
		#get values
		my $x = $arr[0];
		my $y = $arr[1];
		my $z = $arr[2];
		my $lu = $arr[3];
		my $lv = $arr[4];
		my $dip = $arr[5];
		my $strike = $arr[6];
	
		#calculate normal
		my $nz = cos($dip*$pi/180.0);
		my $nxy = sqrt(1.0 - $nz * $nz);
		my $nx = $nxy * cos($strike*$pi/180.0);
		my $ny = -$nxy * sin($strike*$pi/180.0);
		
		#get ux uy uz (vector along strike)
		my $ux = sin($strike*$pi/180.0);
		my $uy = cos($strike*$pi/180.0);
		my $uz = 0;
		
		#get vx vy vz (vector along dip)
		my $vx = $ny*$uz - $nz*$uy;
		my $vy = $nz*$ux - $nx*$uz;
		my $vz = $nx*$uy - $ny*$ux;
		
		#get element loci (+u, +v)
		{
			my $xx = $x + 0.5 * ($ux * $lu + $vx * $lv);
			my $yy = $y + 0.5 * ($uy * $lu + $vy * $lv);
			my $zz = $y + 0.5 * ($uz * $lu + $vz * $lv);
			print "   $inode, $xx, $yy, $zz\n";
			$inode++;
		}
		#get element loci (+u, -v)
		{
			my $xx = $x + 0.5 * ($ux * $lu - $vx * $lv);
			my $yy = $y + 0.5 * ($uy * $lu - $vy * $lv);
			my $zz = $y + 0.5 * ($uz * $lu - $vz * $lv);
			print "   $inode, $xx, $yy, $zz\n";
			$inode++;
		}	
		#get element loci (-u, -v)
		{
			my $xx = $x - 0.5 * ($ux * $lu + $vx * $lv);
			my $yy = $y - 0.5 * ($uy * $lu + $vy * $lv);
			my $zz = $y - 0.5 * ($uz * $lu + $vz * $lv);
			print "   $inode, $xx, $yy, $zz\n";
			$inode++;
		}
		#get element loci (-u, +v)
		{
			my $xx = $x - 0.5 * ($ux * $lu - $vx * $lv);
			my $yy = $y - 0.5 * ($uy * $lu - $vy * $lv);
			my $zz = $y - 0.5 * ($uz * $lu - $vz * $lv);
			print "   $inode, $xx, $yy, $zz\n";
			$inode++;
		}
		$nfaults++;
    }
}
close IN;

print "*ELEMENT, TYPE=S4R, ELSET=EB1\n";
my $ielement = 1;
for(my $i = 0; $i < $nfaults; $i++)
{
	my $inode0 = 4*$i + 1;
	my $inode1 = 4*$i + 2;
	my $inode2 = 4*$i + 3;
	my $inode3 = 4*$i + 4;
	print "   $ielement, $inode0, $inode1, $inode2, $inode3\n";
	$ielement++;
}
