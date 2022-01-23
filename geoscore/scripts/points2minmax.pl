#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <filenames ...>\n" if ( $narg < 1 );

my @min = ( 1e100,  1e100,  1e100 );
my @max = ( -1e100, -1e100, -1e100 );

for ( my $i = 0 ; $i < $narg ; $i++ ) {
	my $filename = shift;
	die "$0 : cannot open $filename!\n" if(!open( IN, "<$filename" ));
	while (<IN>) {
		my $line = $_;
		$line=~s/^\s+//;
		$line=~s/\s+$//;
		my @arr = split(/\s+/, $line );
		for ( my $j = 0 ; $j < 3 ; $j++ ) {
			if ( $max[$j] < $arr[$j] ) {
				$max[$j] = $arr[$j];
			}
			if ( $min[$j] > $arr[$j] ) {
				$min[$j] = $arr[$j];
			}
		}
	}
	close IN;
}

for ( my $k = 0 ; $k < 3 ; $k++ ) {
	print "$k: $min[$k] $max[$k]\n";
}
