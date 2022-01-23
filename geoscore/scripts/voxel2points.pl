#!/bin/perl
use strict; 

my $narg = $#ARGV + 1;
die "usage: $0 <voxel filename> <x filename> <y filename> <z filename> <opt: t filename>\n" if($narg <4 || $narg > 5);

my $filename = shift;
die "cannot open $filename" if(!open(IN,"<$filename"));


###########################
#fill coordinate arrays
###########################
my @xs = ();
my @ys = ();
my @zs = ();
my @ts = ();
for(my $i = 0; $i < ($narg-1); $i++)
{
    my $fn = shift;
    die "cannot open $fn" if(!open(VIN,"<$fn"));
    while(<VIN>)
    {
	my $xx = $_;
	$xx=~s/^\s+//;
	$xx=~s/\s+$//;
	my @arr = split(/\s+/,$xx);
	for(my $ii=0; $ii <= $#arr; $ii++)
	{
	    $xx = $arr[$ii]; 
	    $xx *= 1.0;
	    if($i==0)
	    {
		push(@xs, $xx);
	    }
	    elsif($i==1)
	    {
		push(@ys, $xx);
	    }
	    elsif($i==2)
	    {
		push(@zs, $xx);
	    }
	    else
	    {
		push(@ts, $xx);
	    }
	}
    }
    close VIN;
}
if($#ts==-1)
{
    push(@ts, 0);
}

###########################
#check counts
###########################
{
    my $line = <IN>;
    chomp($line);
    $line=~s/^\s+//;
    my @arr = split(/\s+/,$line);
    print "READING VOXEL HEADING: ";
    for(my $i = 0; $i <= $#arr; $i++)
    {
	print "  $i -> $arr[$i]\n";
	die "too many x entries expected in $filename: $arr[$i] not ".($#xs+1)."\n" if($i==0 && ($#xs+1) != $arr[$i]);
	die "too many y entries expected in $filename: $arr[$i] not ".($#ys+1)."\n" if($i==1 && ($#ys+1) != $arr[$i]);
	die "too many z entries expected in $filename: $arr[$i] not ".($#zs+1)."\n" if($i==2 && ($#zs+1) != $arr[$i]);
	die "too many t entries expected in $filename: $arr[$i] not ".($#ts+1)."\n" if($i==3 && ($#ts+1) != $arr[$i]);
	die "too many entries expected in header of $filename: ".($#arr + 1)."\n" if($i>3);
    }
}


###########################
#create points
###########################
my @pparr = ();
for(my $it = 0; $it <= $#ts; $it++)
{
    open(OUT,sprintf(">tmp_%04d",$it));
    open(OUT3,sprintf(">tmp_%04d.3D",$it));
    print OUT3 "x y z var\n";
    for(my $k = 0; $k <= $#zs; $k++)
    {
	my $z = $zs[$k];
	for(my $j = 0; $j <= $#ys; $j++)
	{
	    my $y = $ys[$j];
	    for(my $i = 0; $i <= $#xs; $i++)
	    {
		my $x = $xs[$i];
		if($#pparr<0)
		{
		    my $line = <IN>;
		    chomp($line);
		    $line=~s/^\s+//;
		    my @arr = split(/\s+/,$line);		    
		    push(@pparr, @arr);
		}
		my $pp = shift(@pparr);
		print OUT "$x $y $z $pp\n";		
		print OUT3 "$x $y $z $pp\n";		
	    }
	}
    }
    close OUT;
    close OUT3;
}
