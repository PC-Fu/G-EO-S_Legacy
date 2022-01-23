#!/opt/local/bin/perl
use strict;

my $narg = $#ARGV + 1;

die "usage: $0 <psuade definition file> <psuade variable input file>\n"
  if ( $narg != 2 );
{
	my $filename = shift;
	die "cannot open $filename\n" if ( !open( IND, "<$filename" ) );
	$filename = shift;
	die "cannot open $filename\n" if ( !open( INV, "<$filename" ) );
}

#set the defaults

my %globals = ();

$globals{"G"}          = 3e10;
$globals{"K"}          = 5e10;
$globals{"SHEARRATE"}  = 3.1709792e-09;
$globals{"SIGHSTRIKE"} = 135;
$globals{"TBURN"}      = 20;

$globals{"FAULTSTRIKE"} = 47;
$globals{"FAULTLENGTH"} = 1690;
$globals{"FAULTDIP"}    = 90;
$globals{"FAULTHEIGHT"} = 1515;
$globals{"NDIP"}        = 18;
$globals{"NSTRIKE"}     = 20;

$globals{"XMIN"} = 0;
$globals{"YMIN"} = 0;
$globals{"ZMIN"} = -2490;

$globals{"XMAX"} = 4000;
$globals{"YMAX"} = 4000;
$globals{"ZMAX"} = -800;

$globals{"XORG"} = 787;
$globals{"YORG"} = 360;
$globals{"ZTOP"} = -885;

$globals{"NVOXX"} = 20;
$globals{"NVOXY"} = 20;
$globals{"NVOXZ"} = 20;

my %randoms = ();

$randoms{"Dc"}                  = [ 0.000015, 0.000035 ];
$randoms{"frictionCoefficient"} = [ 0.6,      0.9 ];
$randoms{"A"}                   = [ 0.005,    0.008 ];
$randoms{"B"}                   = [ 0.015,    0.015 ];

$randoms{"currentFrictionCoefficient"} = [ 0.6,  0.9 ];
$randoms{"alpha"}                      = [ 0.25, 0.25 ];
$randoms{"shearRateStar"}              = [ 1e-6, 1e-6 ];
$randoms{"shearSlipRateAB"}            = [ 1e-6, 1e-6 ];
$randoms{"stateStar"}                  = [ 0.1,  0.1 ];
$randoms{"state"}                      = [ 1e8,  1e8 ];
$randoms{"shearRateDrive"} =
  [ $globals{"SHEARRATE"}, $globals{"SHEARRATE"} ];
$randoms{"stressRate"} = [ -1e-18, -1e-18 ];


#get variable names
my @vars = ();
{
	while (<IND>) {
		if (/^\s*variable\s+(\S+)\s+(\S+)/) {
			push( @vars, $2 );
			my $ok = 0;
			my $varName = $vars[$#vars];
			if ( $varName =~ m/ABratio/ ) {
			    $ok = 1;
			}
			elsif ( exists( $globals{$vars[$#vars]} ) ) {
			    $ok = 1;
			}
			else
			{
			    $varName=~s/Max//;
			    $varName=~s/Min//;
			    if ( exists( $randoms{$vars[$#vars]} ) ) {
				$ok = 1;
			    }
			}
			die "Variable $vars[$#vars] is not an allowable variable name for GPAC! Exiting." if($ok == 0);
		}
		elsif (/END/) {
			last;
		}
	}
}
close IND;

#get variable values
{
	my $num = <INV>;
	for ( my $i = 0 ; $i < $num ; $i++ ) {

		#VALUE
		my $val = <INV>;
		chomp($val);
		$val =~ s/^\s+//;

		#VARIABLE NAME
		my $varName = $vars[$i];
		chomp($varName);
		$varName=~s/^\s+//;
		if ( $varName =~ m/ABratio/ ) {
			my $B = 0.015000;
			my $A = $val * $B;
			$randoms{"A"} = [ $A, $A ];
			$randoms{"B"} = [ $B, $B ];
		}
		elsif ( $varName =~ m/Min/ ) {
			$varName =~ s/Min//;
			$randoms{$varName}[0] = $val;
		}
		elsif ( $varName =~ m/Max/ ) {
			$varName =~ s/Max//;
			$randoms{$varName}[1] = $val;
		}
		else {
			if ( exists( $globals{$varName} ) ) {
			    $globals{$varName} = $val;
			}
			elsif ( exists( $randoms{$varName} ) ) {
			    $randoms{$varName} = [$val, $val];
			}
		}
	}
}
close INV;

#print the file

#2) go through all globals and print
my $key = "";
foreach $key ( sort(keys %globals) ) {
    print "GV $key $globals{$key}\n";
}
print "\n\n";

#3) go through all randoms and print
foreach $key ( sort(keys %randoms) ) {
    print "IC $key $randoms{$key}[0] $randoms{$key}[1]\n";
}
