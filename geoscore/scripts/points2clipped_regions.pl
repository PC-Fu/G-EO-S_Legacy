#!/bin/perl
use strict;

my $narg = $#ARGV + 1;

die "usage: $0 <points file> <polyhedra file>\n" if($narg != 2);

my $fnamep = shift;
die "Cannot open ARG1: $fnamep" if(!open(IN,"<$fnamep"));

my $fname = shift;
die "Cannot open ARG2: $fname" if(!open(POLY,"<$fname"));
my @polys = ();
my $i = 0;
while(<POLY>)
{
    $i++;
    my $line = $_;
    $line=~s/^\s+//;
    $line=~s/\s+$//;
    my @arr = split(/\s+/,$line);
    die "Cannot read $fname: line $i - $line\n" if($#arr != 35);
    push(@polys,\@arr);
}
close POLY;

$i = 0;
while(<IN>)
{
    $i++;
    my $line = $_;
    $line=~s/^\s+//;
    $line=~s/\s+$//;
    my @arr = split(/\s+/,$line);
    die "Cannot read $fnamep: line $i - $line\n" if($#arr < 2);
    my $process = &ProcessPoly(\@arr, \@polys);
    print "$line $process\n";
}
close IN;

sub ProcessPoly
{
    my $xptr = shift;
    my $pptr = shift;
    my $i = 0;
    foreach my $poly (@$pptr)
    {
	my $ok = &ProcessSinglePoly($xptr, $poly);
	if($ok == 1)
	{
	    return $i;
	}
	++$i;
    }
    return -1;
}

sub ProcessSinglePoly
{
    my $xptr = shift;
    my @pt=@$xptr;
    my $pptr = shift;
    my @p = @$pptr;
    #print STDERR "p[0] $p[0]\n";
    for(my $i = 0; $i < 6; $i++)
    {
	my $iix = 6 * $i; 
	my $iin = 6 * $i + 3;
	my @x = ($p[$iix],$p[$iix+1],$p[$iix+2]);
	my @n = ($p[$iin],$p[$iin+1],$p[$iin+2]); 
	my $dot = 0;
	for(my $j = 0; $j < 3; $j++)
	{
	    my $dx = $pt[$j] - $x[$j];
	    $dot += $n[$j] * $dx;
	    #print STDERR "FACE ($i): DIM $j: $pt[$j] - $x[$j] = $dx --> * $n[$j]\n";
	}
	if($dot > 0)
	{
	    return 0;
	}
    } 
    return 1;
}
