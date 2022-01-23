#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <GPAC standard output file> <opt: time 0 in yr>\n"
  if ( $narg < 1 );

my $filename = shift;
die "cannot open $filename" if ( !open( IN, "<$filename" ) );

my $time0 = 0;
if ( $narg > 1 )
{
	$time0 = shift;
}

my $sec_to_yr = 1.0 / 31557600;
while (<IN>) {
	if (/^eq_event/) {
		my $line = $_;
		$line =~ s/^\s+//;
		chomp($line);
		my @arr = split( /\s+/, $line );

		#eq_event: t M E A d x y z
		die "$0 : cannot parse line: $_"
		  if ( $#arr != 8 );
		$arr[1] *= $sec_to_yr;
		$arr[1] -= $time0;
		print "$arr[6] $arr[7] $arr[8] $arr[1] $arr[2] $arr[3] $arr[4] $arr[5]\n";
	}
}
close IN;
