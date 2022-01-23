#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;

#determine if the user is asking for help (either "help" or no arguments)
my $print_help = 0;
my $file       = "";
if ( $narg > 0 ) {
	$file = shift;
	if ( $file eq "help" ) { $print_help = 1; }
}
else {
	$print_help = 1;
}

if ($print_help) {
	print "usage: $0 <input_file> <scalar0> <scalar1> ... <scalarN>\n";
	print "input is space-delimited, carriage return terminated\n";
	print
	  "note: positions must be defined by the field names 'x', 'y', and 'z'\n";
	print "      also, output to a .okc file for visit to auto-recognize\n";
	exit();
}

my $nscalars = $narg - 1;

#first, assume 12 characters per scalar
my $nchar = 12;
my $ncols = $nscalars * $nchar - 1;

my @args      = ();
my @positions = ();
my %argorder  = ();

if ( $nscalars == 0 ) {

	#by default, assume the file only contains points
	push( @args, "x" );
	push( @args, "y" );
	push( @args, "z" );
}
else {
	my @args_tmp = ();
	for ( my $i = 0 ; $i < $nscalars ; $i++ ) {
		my $arg = shift;	
		push(@args_tmp, $arg);
	}	
	my $j = 0;
	for ( my $i = 0 ; $i < $nscalars ; $i++ ) {
		if ( $args_tmp[$i] eq "x" ) {
			$positions[$j] = $i;
			push( @args, $args_tmp[$i] );
			$j++;
		}
	}
	for ( my $i = 0 ; $i < $nscalars ; $i++ ) {
		if ( $args_tmp[$i] eq "y" ) {
			$positions[$j] = $i;
			push( @args, $args_tmp[$i] );
			$j++;
		}
	}
	for ( my $i = 0 ; $i < $nscalars ; $i++ ) {
		if ( $args_tmp[$i] eq "z" ) {
			$positions[$j] = $i;
			push( @args, $args_tmp[$i] );
			$j++;
		}
	}
	for ( my $i = 0 ; $i < $nscalars ; $i++ ) {
		my $arg = $args_tmp[$i];
		my $ispos = 0;
		if ( $arg eq "x" || $arg eq "y" || $arg eq "z" ) {
			$ispos = 1;
		}
		if ( $ispos == 0 ) {
			$positions[$j] = $i;
			push( @args, $arg );
			$j++;
		}
	}
}

#count points
my $npoint = 0;
{
	open( IN, "<$file" );
	while (<IN>) {
		$npoint++;
	}
	close IN;
}

#print the file header
#
# NOTE: XMDV FILE FORMAT: (N = the number of variables; R = the number of rows; C = the number of columns)
#
#N R 12
#varname0
#...
#varnameN-1
#var0minvalue   var0maxvalue    10
#...
#varN-1minvalue varN-1maxvalue  10
#data0,0   data1,0   data2,0  ...  dataC-1,0
#...
#data0,R-1 data1,R-1 data2,R-1 ... dataC-1,R-1

print "$nscalars $npoint $ncols\n";
for ( my $i = 0 ; $i < $nscalars ; $i++ ) {
	print "$args[$i]\n";
}

#get min max of each dimension
{
	my @mins = ();
	my @maxs = ();
	for ( my $ii = 0 ; $ii < $nscalars ; $ii++ ) {
		push( @mins, 1e100 );
		push( @maxs, -1e100 );
	}
	{
		die "cannot find $file\n" if ( !open( IN, "<$file" ) );
		while (<IN>) {
			my $line = $_;
			chomp($line);
			$line =~ s/^\s+//;
			my @values = split( /\s+/, $line );
			for ( my $i = 0 ; $i <= $#values ; $i++ ) {
				my $value = $values[$i];
				if ( $value < $mins[$i] ) {
					$mins[$i] = $value;
				}
				if ( $value > $maxs[$i] ) {
					$maxs[$i] = $value;
				}
			}
		}
		close IN;
	}
	for ( my $ii = 0 ; $ii < $nscalars ; $ii++ ) {
		print "$mins[$positions[$ii]] $maxs[$positions[$ii]] 10\n";
	}
}

#print points
{
	open( IN, "<$file" );
	while (<IN>) {
		my $line = $_;
		chomp($line);
		$line =~ s/^\s+//;
		my @values = split( /\s+/, $line );
		for ( my $i = 0 ; $i < $nscalars ; $i++ ) {
			my $value = $values[ $positions[$i] ];
			if ( $value < 0 ) {
				print sprintf( " %1.4e", $value );
			}
			else {
				print sprintf( " %1.5e", $value );
			}
		}
		print "\n";
	}
	close IN;
}
