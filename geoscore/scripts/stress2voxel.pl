#!/bin/perl
use strict;
use Math::Trig;

my $narg = $#ARGV + 1;
die
"usage: $0 <profile file name with z sigH sigh sigv tuples> <sigH strike> <fault strike> <fault dip>\n"
  if ( $narg != 4 );

my $filename = shift;
die "cannot open $filename" if ( !open( IN, "<$filename" ) );

#get transform, where theta = angle between fault strike and sigH strike
my $ctst = 0;    #cos (theta) * sin(theta)
my $ct2  = 0;    #cos^2(theta)
my $st2  = 0;    #sin^2(theta)

my $cdsd = 0;    #cos(dip) * sin(dip)
my $cd2  = 0;    #cos^2(dip)
my $sd2  = 0;    #sin^2(dip)

my @vstrike = ( 0, 0, 0 );
my @vdip    = ( 0, 0, 0 );
{
	my $fct     = atan(1.0) / 45;
	my $sstrike = shift;
	$sstrike *= $fct;
	my $fstrike = shift;
	$fstrike *= $fct;
	my $fdip = shift;
	$fdip *= $fct;
	my $theta = $fstrike - $sstrike;

	my $ct = cos($theta);
	my $st = sin($theta);
	$ct2  = $ct * $ct;
	$st2  = $st * $st;
	$ctst = $ct * $st;

	my $cd = cos($fdip);
	my $sd = sin($fdip);
	$cd2  = $cd * $cd;
	$sd2  = $sd * $sd;
	$cdsd = $cd * $sd;

	$vstrike[0] = sin($fstrike);
	$vstrike[1] = cos($fstrike);
	$vdip[2]    = -$sd;
	$vdip[0]    = $vstrike[1] * $cd;
	$vdip[1]    = -$vstrike[0] * $cd;
}

##################################
#PRINT THE (UNUSED) COORDINATES IN X AND Y
##################################
{
	open( XOUT, ">xx" );
	print XOUT "0\n";
	close XOUT;

	open( YOUT, ">yy" );
	print YOUT "0\n";
	close YOUT;
}

##################################
#PRINT THE COORDINATES IN Z
##################################

open( T, ">tau" );
print T "1 1 ";
open( S, ">sigma" );
print S "1 1 ";
open( ZOUT, ">zz" );

#get sigma entries from file ... will clip to fit domain later
my @asigs = ();
my $z     = -1e10;
my $nsigs = 0;
while (<IN>) {
	my $line = $_;
	$line =~ s/^\s+//;
	chomp($line);
	my @arr = split( /\s+/, $line );
	die "$0 : cannot parse line: $_" if ( $#arr != 3 );
	die
	  "$0 : cannot have non-sequential points in the sigma file: $arr[0] <= $z"
	  if ( $arr[0] <= $z );
	$z = $arr[0];
	print ZOUT "$z\n";
	$nsigs++;
	die "$0 : cannot have negative sigmas\n"
	  if ( $arr[1] < 0 || $arr[2] < 0 || $arr[3] < 0 );
	push( @asigs, \@arr );
}
close IN;
close ZOUT;

die "must have at least two sigma entries in $filename\n"
  if ( $nsigs < 2 );

print T "$nsigs\n";
print S "$nsigs\n";

my $ir2 = 1.0 / sqrt(2);
for ( my $i = 0 ; $i < $nsigs ; $i++ ) {
	my $z    = $asigs[$i][0];
	my $sigH = $asigs[$i][1];
	my $sigh = $asigs[$i][2];
	my $sigv = $asigs[$i][3];

	#get the sigma normal in the strike frame
	my $sig_ns = $sigh * $ct2 + $sigH * $st2;

	#get the sigma normal in the plane frame
	my $sig = $sigv * $cd2 + $sig_ns * $sd2;
	print S "$sig\n";

	#get the shear stress in the strike direction
	my $taus = ( $sigH - $sigh ) * $ctst;

	#get the shear stress in the dip direction
	my $taud = ( $sigv - $sig_ns ) * $cdsd;

	#get the magnitude of the total shear stress
	my $tau = sqrt( $taud * $taud + $taus * $taus );

	#print the value
	print T "$tau\n";

	my $fctd = $ir2;
	my $fcts = $ir2;
	if ( $tau > 0 ) {
		$fctd = $taud / $tau;
		$fcts = $taus / $tau;
	}
	my @dislocation = ( 0, 0, 0 );
	my $sum = 0.0;
	for ( my $j = 0 ; $j < 3 ; $j++ ) {
		$dislocation[$j] = $fctd * $vdip[$j] + $fcts * $vstrike[$j];
		$sum += $dislocation[$j] * $dislocation[$j];
	}
	if ( $sum > 0 ) {
		$sum = 1.0 / sqrt($sum);
	}
	for ( my $j = 0 ; $j < 3 ; $j++ ) {
		$dislocation[$j] *= $sum;
	}
	print
"\t\t\t   rakeVector=\"$dislocation[0] $dislocation[1] $dislocation[2]\"\n";
}
