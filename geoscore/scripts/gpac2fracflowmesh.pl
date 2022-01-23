#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <GPAC input filename>\n" if($narg != 1);

my $fname = shift;
die "cannot open $fname\n" if ( !open( IN, "<$fname" ) );

#get the necessary file names
my $meshfile = "";
my $propfile = "";
while (<IN>) {
	if (/Mesh\s+file\s*\=\s*\"(\S+)\"/) {
		$meshfile = $1;
	}
	elsif (/ElementPropertyPath\s*\=\s*\"(\S+)\"/) {
		$propfile = $1;
	}
}
close IN;

#get the properties in a hash
my %properties = ();
die "cannot open $propfile\n" if ( !open( IN, "<$propfile" ) );
while (<IN>) {
	my $line = $_;
	$line =~ s/^\s+//;
	chomp($line);
	my @arr = split( /\,\s*/, $line );
	die "cannot read property: $line\n" if ( $#arr != 2 );
	$properties{ $arr[0] } = [ $arr[1], $arr[2] ];
}
close IN;

#get the nodes and elements then write to FRACFLOW format
my $inode = 1;
my %nodes = ();
die "cannot open $meshfile\n" if ( !open( IN, "<$meshfile" ) );
my $node_on = 0;
my $el_on   = 0;
while (<IN>) {
	if (/^\s*\*NODE/) {
		$node_on = 1;
		$el_on   = 0;
	}
	elsif (/^\s*\*ELEMENT/) {
		$el_on   = 1;
		$node_on = 0;
	}
	elsif (/^\s*\*/) {
		$el_on = 0;
		$node_on = 0;
	}
	elsif ( $node_on == 1 ) {
		my $line = $_;
		chomp($line);
		$line =~ s/^\s+//;
		my @arr = split( /\,\s*/, $line );
		die "node entry not recognized: $line\n" if ( $#arr < 2 );
		$nodes{ $arr[0] } = [ $inode, $arr[1], $arr[2] ];
		print "nod $arr[1] $arr[2]\n";
		#note: z is ignored
		$inode++;
	}
	elsif ( $el_on == 1 ) {
		my $line = $_;
		chomp($line);
		$line =~ s/^\s+//;
		my @arr = split( /\,\s*/, $line );
		die "element entry not recognized: $line\n" if ( $#arr != 3 );
		my $n0 = ${$nodes{$arr[1]}}[0];
		my $n1 = ${$nodes{$arr[2]}}[0];
		my $n2 = ${$nodes{$arr[3]}}[0];
		my $p0 = ${$properties{$arr[0]}}[0];
		my $p1 = ${$properties{$arr[0]}}[1];
		print "tri $n0 $n1 $n2 $p0 $p1\n";
	}
}
print "END\n";
close IN;
