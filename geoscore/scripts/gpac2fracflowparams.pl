#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <GPAC input filename>\n" if($narg != 1);

my $fname = shift;
die "cannot open $fname\n" if ( !open( IN, "<$fname" ) );

#get the necessary items
my $params_on = 0;
while (<IN>) {
    if(/^\s*\</)
    {
	if($params_on == 1)
	{
	    print "END\n";
	    exit 0;
	}
    }
    if(/^\s*\<Small/)
    {
	$params_on = 1;
    }
    elsif($params_on == 1)
    {
	if(/Path/)
	{
	}
	elsif(/\>/)
	{
	}
	else
	{
	    my $line = $_;
	    chomp($line);
	    my @arr = split(/\=/,$line);
	    my $name = $arr[0];
	    $name=~s/\s+//g;
	    my $value = $arr[1];
	    $value=~s/\"//g;
	    print "$name $value\n";
	}
    }
}
close IN;
exit 1;
