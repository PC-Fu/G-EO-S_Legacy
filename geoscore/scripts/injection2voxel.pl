#!/bin/perl

#Ported from single_fault.m (Josh White)

use strict;

my $narg = $#ARGV + 1;
die
"usage: $0 <x-file> <y-file> <z-file> <offset t in years> <x-coord> <y-coord> <reservoir thickness>"
  if ( $narg != 7 );

#############################
#READ X Y Z T0
my @xs = ();
my @ys = ();
my @zs = ();
for ( my $i = 0 ; $i < 3 ; $i++ ) {
	my $fn = shift;
	die "$0 : cannot open $fn\n" if ( !open( IN, "<$fn" ) );
	while (<IN>) {
		my $line = $_;
		$line =~ s/^\s+//;
		chomp($line);
		if ( $i == 0 ) {
			push( @xs, $line );
		}
		elsif ( $i == 1 ) {
			push( @ys, $line );
		}
		else {
			push( @zs, $line );
		}
	}
	close IN;
}
my $time0 = shift;

#convert time from years to hours
$time0 *= 24 * 365.25;

my $xinj = shift;
my $yinj = shift;
my $dz = shift;

#############################
#SET RESERVOIR PARAMS

my $reservoir_depth = 2000;     # meters
my $reservoir_k     = 200;      # perm, millidarcy
my $reservoir_mu    = 0.30;     # viscosity, centipoise (0.30 brine, 0.05 co2)
my $reservoir_phi   = 0.10;     # porosity
my $reservoir_c     = 0.025;    # total compressibility, 1/MPa
my $reservoir_b     = 1;        # formation volume factor, Rm3/Sm3
my $reservoir_rw    = 0.1;      # wellbore radius, meters
my $reservoir_h     = $dz;      # reservoir thickness, meters

my $q    = 250;                 # injection rate, m3/hr
my $time = 50 * 365.25 * 24;    # injection time, 50 years in hours

my @injectors_rate     = ( $q,    0 );         # inject at Q and then stop
my @injectors_duration = ( $time, 3 * $time ); # inject for 50, drawdown for 150

my @injectors_xy = ( $xinj, $yinj );     # one injector 

my $hydrostatic = 0.01;                        # MPa/m

#############################
#CREATE VOXEL FILE

my @times = ( 0 );
if($time0 > 0)
{
    push(@times, $time0);
}
my $ishutin = -1;
for ( my $id = 0 ; $id <= $#injectors_duration ; $id++ ) {
	my $tlast = $times[$#times];
	$ishutin = $#times;
	for ( my $it = 0 ; $it < 11 ; $it++ ) {
		my $tt = $injectors_duration[$id] * exp( -5 + 0.5 * $it ) + $tlast;
		push( @times, $tt );
	}
}

my $nt = $#times + 1;
my $nx = $#xs + 1;
my $ny = $#ys + 1;
my $nz = $#zs + 1;

open( OUT,  ">voxel" );
open( TOUT, ">t" );
print OUT "$nx $ny $nz $nt\n";
my @dps = ();
my $ndp = ($#ys + 1) * ($#xs + 1);
$dps[$ndp-1] = 0;
for ( my $it = 0 ; $it <= $#times ; $it++ ) {

	#convert time in hours to SI seconds
	my $tts = $times[$it] * 3600;
	print TOUT "$tts\n";

	#convert to the pseudo time
	my $qcurr = $injectors_rate[0];
	my $tt = $times[$it];
	if ( $it > 0 ) {
	        if ( $it > $ishutin ) {
			$qcurr = $injectors_rate[1] - $injectors_rate[0];
			$tt -= $times[$ishutin];
		}
		elsif ( $it > 0 ) {
			$tt -= $times[1];
		}
		print "QCURR: $qcurr TT: $tt IT: $it ($ishutin)\n";
	}

	#set the pressure at each point
	my $idp = 0;
	for ( my $k = 0 ; $k <= $#zs ; $k++ ) {
		my $zz = $zs[$k];
		my $phydro = $zz * ( -$hydrostatic );
		for ( my $j = 0 ; $j <= $#ys ; $j++ ) {
			my $yy = $ys[$j] - $injectors_xy[1];
			$yy *= $yy;
			for ( my $i = 0 ; $i <= $#xs ; $i++) {
				my $dp = 0;
				if ( $k == int( ($#zs + 1) / 2 ) ) {
				    if($tt > 0 && $it < $#times) {
					my $xx = $xs[$i] - $injectors_xy[0];
					$xx *= $xx;
					my $radius = sqrt( $xx + $yy );
					if ( $radius < $reservoir_rw ) {
						$radius = $reservoir_rw;
					}
					my $dpk = &PressureContribution(
					    $radius,       $tt,          $qcurr,
					    $reservoir_rw, $reservoir_k, $reservoir_phi,
					    $reservoir_mu, $reservoir_c );
					$dps[$idp] += $dpk;
					#print "DP: $xx $yy $dpk\n";
				    }
				    $dp += $dps[$idp];
				    $idp++;
				}
				$dp += $phydro;
				$dp *= 1e6;
				print OUT "$dp ";
			}
			print OUT "\n";
		}
	}
}
close OUT;
close TOUT;

sub PressureContribution {
	my $r = shift;
	my $t = shift;
	my $q = shift;

	my $reservoir_rw  = shift;
	my $reservoir_k   = shift;
	my $reservoir_phi = shift;
	my $reservoir_mu  = shift;
	my $reservoir_c   = shift;

	my $psi    = 0.00689475729;    # MPa
	my $feet   = 0.3048;           # m
	my $barrel = 0.158987295;      # m3

	my $A =
	  0.0002637 * $feet *
	  $feet /
	  $psi;    # goofy coefficients to make inconsistent units work
	my $B = 141.2 * 24 / $barrel * $feet * $psi;

	my $rd = $r / $reservoir_rw;

	my $tmp1 = $A * $reservoir_k * $t;
	my $tmp2 =
	  $reservoir_phi *
	  $reservoir_mu *
	  $reservoir_c *
	  $reservoir_rw *
	  $reservoir_rw;

	my $td = $tmp1 / $tmp2;

	my $x = -$rd * $rd / ( 4.0 * $td );
	my $pd = &Ei(-$x);
	$pd *= -0.5;

	#die "$tmp1 $tmp2 $td $x $pd\n" if($pd > 1e100);

	$tmp1 = $B * $q * $reservoir_b * $reservoir_mu;
	$tmp2 = $reservoir_k * $reservoir_h;

	my $delta_p = ( $tmp1 / $tmp2 ) * ($pd);

	return $delta_p;
}

#FROM: Numerical Recipes in C
sub Ei {
	my $x = shift;
	die "Bad argument: $x\n" if ( $x <= 0.0 );

	my $EULER = 0.57721566;
	my $MAXIT = 100;
	my $FPMIN = 1.0e-30;
	my $EPS   = 6.0e-8;

	#Special case
	#avoid failure of convergence test because of underflow
	if ( $x < $FPMIN ) {
		return log($x) + $EULER;
	}

	#Use power series
	if ( $x <= -log($EPS) ) {
		my $sum  = 0.0;
		my $fact = 1.0;
		my $k = 0;
		for ( $k = 1 ; $k <= $MAXIT ; $k++ ) {
			$fact *= $x / $k;
			my $term = $fact / $k;
			$sum += $term;
			last if ( $term < $EPS * $sum );
		}
		die "Series failed in ei" if ( $k > $MAXIT );
		return $sum + log($x) + $EULER;
	}

	#Use asymptotic series. Start with second term.
	else {
		my $sum  = 0.0;
		my $term = 1.0;
		for ( my $k = 1 ; $k <= $MAXIT ; $k++ ) {
			my $prev = $term;
			$term *= $k / $x;
			last if ( $term < $EPS );
			#Since final sum is greater than one, 
			#term itself approximates the relative error .
			if ( $term < $prev )
			{
				$sum += $term;
			} else {
				$sum -= $prev;
				last;
			}
		}
		return exp($x) * ( 1.0 + $sum ) / $x;
	}
}
