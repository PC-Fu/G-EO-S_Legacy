#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <abaqus file name> <path point file> <fault top> <fault bottom> <opt: u-axis (strike direction)> <opt: v-axis (dip direction)> (narg = $narg currently)\n"
  if ( $narg < 4 || $narg > 6);

#--------------------------------
#READ IN ARGUMENTS
#--------------------------------

#ABAQUS IN FILE
{
	my $fname = shift;
	die "cannot open $fname\n" if ( !open( IN, "<$fname" ) );
}

#POINTS IN FILE AND CUMULATIVE ARC LENGTH
my @pts = ();
my @cumulativearclength = (0);
{
	my $fname = shift;
	die "cannot open $fname\n" if ( !open( INPTS, "<$fname" ) );
	my @lastpt = ();
	my $arclength = 0;
	while(<INPTS>)
	{
		my $line=$_;
		$line =~s/^\s+//;
		chomp($line);
		my @arr = split(/\s+/, $line);
		die "should be a space-delimited 2D point: $#arr --> \"$line\"" if($#arr != 1);
		push(@pts, \@arr);
		if($#lastpt > -1)
		{
			my $dx = $arr[0] - $lastpt[0];
			$dx *= $dx;
			my $dy = $arr[1] - $lastpt[1];
			$dy *= $dy;
			$dx += $dy;
			$dx = sqrt($dx);
			die "arc length increment must be non-negative (dx = $dx): $#cumulativearclength" if($dx <= 0);
			$arclength += $dx;
			push(@cumulativearclength, $arclength);
		}
		$lastpt[0] = $arr[0];
		$lastpt[1] = $arr[1];
	}
	close INPTS;
}

#GET FAULT TOP/BOTTOM BOUNDS
my $top = shift;
my $bottom = shift;
die "top < bottom! ($top < $bottom)" if($top < $bottom);

#SET AXES
my $u = 0;
my $v = 1;
my $w = 2;
{
	if($narg > 2)
	{
		$u = shift;
		die "u must be integer in [0,2]" if($u != 0 && $u != 1 && $u != 2);
	}
	if($narg > 2)
	{
		$v = shift;
		die "u must be integer in [0,2]" if($v != 0 && $v != 1 && $v != 2);
	}
	die "u must not be v: $u == $v" if($u == $v);
	if($u == 2 || $v == 2)
	{
		$w = 1;
	}
	if($u == 1 || $v == 1)
	{
		$w = 0;
	}
}

#--------------------------------
#READ ABAQUS (FIRST TIME)
#--------------------------------
my $xumax = -1e100;
my $xvmax = -1e100;
my $xumin = 1e100;
my $xvmin = 1e100;
{
	my $node_on = 0;
	my $icurr   = 0;
	while (<IN>) {
		if (/^\s*\*\s*(\S+)/) {
			if (/NODE/) {
				$node_on = 1;
			}
			else {
				$node_on = 0;
			}
		}
		else {
			if ($node_on) {
				my $line = $_;
				chomp($line);
				$line =~ s/^\s+//;
				my @arr = split( /\,\s*/, $line );
				die "node_on but this doesn't look like a node: $_" if ( $#arr != 3 );

				my @x = ( $arr[1], $arr[2], $arr[3] );
				if($x[$u] < $xumin)
				{
					$xumin = $x[$u];
				}
				if($x[$u] > $xumax)
				{
					$xumax = $x[$u];
				}
				if($x[$v] < $xvmin)
				{
					$xvmin = $x[$v];
				}
				if($x[$v] > $xvmax)
				{
					$xvmax = $x[$v];
				}
			}
		}
	}
}
#die "xu_max: $xumax
#xv_max: $xvmax
#xu_min: $xumin
#xv_min: $xvmin
#";

#--------------------------------
#READ ABAQUS (SECOND TIME)
#--------------------------------
{
	my $node_on = 0;
	my $icurr   = 0;
	seek(IN, 0, 0);
	while (<IN>) {
		if (/^\s*\*\s*(\S+)/) {
			if (/NODE/) {
				$node_on = 1;
			}
			else {
				$node_on = 0;
			}
			print $_;
		}
		else {
			if ($node_on) {
				my $line = $_;
				chomp($line);
				$line =~ s/^\s+//;
				my @arr = split( /\,\s*/, $line );
				die "node_on but this doesn't look like a node: $_" if ( $#arr != 3 );

				my @x = ( $arr[1], $arr[2], $arr[3] );
				
				#get the v-component
				{
					#normalize the v-component
					my $xv = ($x[$v] - $xvmin) / ($xvmax - $xvmin);
					#transform to fault
					$xv *= ($top - $bottom);
					$xv += $bottom;
					#die "$xv $x[$v]\n";
					$x[$v] = $xv;
				}				

				#normalize the u-component
				my $xu = ($x[$u] - $xumin) / ($xumax - $xumin);

				#find the associated fault strike direction coordinates
				{	
					#get the current point's associated arclength				
					my $arc = $xu * $cumulativearclength[$#cumulativearclength];
					
					#get the cumulative arc lengths on either side
					{
						my $i1 = -1;
						my $arc1 = 0;
						for(my $i = 0; $i <= $#cumulativearclength; $i++)
						{
							if($cumulativearclength[$i] >= $arc)
							{
								$i1 = $i;
								$arc1 = $cumulativearclength[$i];
								last;
							}
						}
						die "impossible: i1 < 0 or i1 > $#cumulativearclength: $i1 ($arc: $cumulativearclength[0] $cumulativearclength[$#cumulativearclength])\n" if($i1 < 0 || $i1 > $#cumulativearclength);
						
						my @uw1 = ($pts[$i1][0], $pts[$i1][1]);
						if($i1 > 0)
						{
							#the point doesn't lie on the initial terminus
							my $i0 = $i1 - 1;
							my $arc0 = $cumulativearclength[$i0];

							#normalize the arclength within the interval
							my $fct = $arc;
							$fct -= $arc0;
							$fct /= ($arc1 - $arc0);
							
							#get the bounding points
							my @uw0 = ($pts[$i0][0], $pts[$i0][1]);
							
							#fct is the fractional distance between pts 0 and 1, so set x[$u] and x[$w]
							$x[$u] = $uw0[0] + ($uw1[0] - $uw0[0]) * $fct;
							$x[$w] = $uw0[1] + ($uw1[1] - $uw0[1]) * $fct;					
						}
						else
						{
							#the point doesn't lie on the initial terminus
							$x[$u] = $uw1[0];
							$x[$w] = $uw1[1];					
						}
					}
				}
				
				print "$arr[0], $x[$u], $x[$w], $x[$v]\n";
			}
			else {
				print $_;
			}
		}
	}
}
close IN;
