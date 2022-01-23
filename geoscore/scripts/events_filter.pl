#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <GPAC standard output file> <time 0 - yr> <dt - yr>\n"
  if ( $narg != 3 );

my $filename = shift;
die "cannot open $filename" if ( !open( IN, "<$filename" ) );
my $t0 = shift;
my $t1 = shift;
$t1 += $t0;

my $sec_to_yr = 1.0 / 31557600;
my $on = 0;
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
		if($arr[1] > $t1)
		{
			exit;
		}
		if($arr[1] >= $t0)
		{
			print $_;
		}
	}
}
close IN;
