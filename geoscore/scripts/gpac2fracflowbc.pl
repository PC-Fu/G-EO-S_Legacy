#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <GPAC input filename> <opt: additional nodeset offset>\n"
  if ( $narg < 1 || $narg > 2 );

my $fname  = shift;
my $offset = 0;
if ( $narg > 1 ) {
	$offset = shift;
}
die "cannot open $fname\n" if ( !open( IN, "<$fname" ) );

#get extents of each nodeset
my %nodes = ();
my %nsets = ();
{

	#get mesh filename
	my $meshfile = "";
	while (<IN>) {
		if (/Mesh\s+file\s*\=\s*\"(\S+)\"/) {
			$meshfile = $1;
		}
	}
	close IN;

	die "cannot open $meshfile\n" if ( !open( IN, "<$meshfile" ) );
	my @xmax = ( 0, 0 );
	my @xmin = ( 0, 0 );
	my $node_on = 0;
	my $ns_on   = 0;
	my $nsname  = "";

	while (<IN>) {
		if (/^\s*\*/) {
			if ( $ns_on == 1 ) {
				$nsets{$nsname} = [
					$xmin[0] - $offset,
					$xmin[1] - $offset,
					$xmax[0] + $offset,
					$xmax[1] + $offset
				];
			}
		}
		if (/^\s*\*NODE/) {
			$node_on = 1;
			$ns_on   = 0;
		}
		elsif (/^\s*\*NSET/) {
			for ( my $i = 0 ; $i < 2 ; $i++ ) {
				$xmax[$i] = -1e100;
				$xmin[$i] = 1e100;
			}
			$ns_on   = 1;
			$node_on = 0;
			my $line = $_;
			chomp($line);
			my @arr = split( /\=/, $line );
			$nsname = $arr[1];
		}
		elsif (/^\s*\*/) {
			$ns_on   = 0;
			$node_on = 0;
		}
		elsif ( $node_on == 1 ) {
			my $line = $_;
			chomp($line);
			$line =~ s/^\s+//;
			my @arr = split( /\,\s*/, $line );
			die "node entry not recognized: $line\n" if ( $#arr < 2 );
			$nodes{ $arr[0] } = [ $arr[1], $arr[2] ];
		}
		elsif ( $ns_on == 1 ) {
			my $line = $_;
			chomp($line);
			$line =~ s/^\s+//;
			my @arr = split( /\,\s*/, $line );
			for ( my $i = 0 ; $i <= $#arr ; $i++ ) {
				if ( length( $arr[$i] ) > 0 ) {
					if ( $xmax[0] < ${ $nodes{ $arr[$i] } }[0] ) {
						$xmax[0] = ${ $nodes{ $arr[$i] } }[0];
					}
					if ( $xmax[1] < ${ $nodes{ $arr[$i] } }[1] ) {
						$xmax[1] = ${ $nodes{ $arr[$i] } }[1];
					}
					if ( $xmin[0] > ${ $nodes{ $arr[$i] } }[0] ) {
						$xmin[0] = ${ $nodes{ $arr[$i] } }[0];
					}
					if ( $xmin[1] > ${ $nodes{ $arr[$i] } }[1] ) {
						$xmin[1] = ${ $nodes{ $arr[$i] } }[1];
					}
				}
			}
		}
	}
	if ( $ns_on == 1 ) {
		$nsets{$nsname} = [
			$xmin[0] - $offset,
			$xmin[1] - $offset,
			$xmax[0] + $offset,
			$xmax[1] + $offset
		];
	}
}

#get the tables
my %tables = ();
{
	die "cannot open $fname\n" if ( !open( IN, "<$fname" ) );
	my $table_on   = 0;
	my @keysvalues = ();
	my $tname      = "";
	while (<IN>) {
		if (/\s*\</) {
			if ( $table_on == 1 ) {
				my @arr = ();
				push(@arr, @keysvalues);
				$tables{$tname} = \@arr;
			}
			$table_on = 0;
		}
		if (/\<Table1D/) {
			$table_on = 1;
		}
		if ( $table_on == 1 ) {
			my $i = -1;
			if (/coord/) {
				$i = 0;
			}
			elsif (/value/) {
				$i = 1;
			}
			elsif (/name/) {
				$i = 2;
			}
			if ( $i >= 0 ) {
				my $line = $_;
				$line =~ s/\s+//g;
				chomp($line);
				my @arr1 = split( /\"/, $line );
				if ( $i < 2 ) {
					my @arr = split( /\,/, $arr1[1] );
					$keysvalues[$i] = \@arr;
				}
				else {
					$tname = $arr1[1];
				}
			}
		}
	}
	close IN;
}

#get the boundary conditions and write
my %bc_pumps = ();
print
"#bc type(1 for velocity, and 2 for stress), dimensionI (1 for x and 2 for y), boundaries (x1,x2,y1,y2), value of bc\n";
{
	die "cannot open $fname\n" if ( !open( IN, "<$fname" ) );
	my $bc_on = 0;

	my $bc_type      = 0;
	my $bc_component = 0;
	my $bc_scale     = 1;

	my $bc_field = "";
	my $bc_nset  = 0;
	my $bc_table = "";

	while (<IN>) {
		if (/\s*\</) {
			if ( $bc_on == 1 ) {

				#get the domain extents
				my $x1 = ${ $nsets{$bc_nset} }[0];
				my $y1 = ${ $nsets{$bc_nset} }[1];
				my $x2 = ${ $nsets{$bc_nset} }[2];
				my $y2 = ${ $nsets{$bc_nset} }[3];

				#write the boundary condition
				if ( $bc_type < 3 ) {
					my $bc_value =
					  $bc_scale * ${ ${ $tables{$bc_table} }[1] }[0];
					print
					  "bc $bc_type $bc_component $x1 $x2 $y1 $y2 $bc_value\n";
				}
				else {
					if ( $bc_field =~ m/pressure/ )
					{
						print "bc $bc_type $bc_component $x1 $x2 $y1 $y2 0\n";
					}
					if ( exists $bc_pumps{$bc_component} ) {
						for (
							my $i = 0 ;
							$i <= $#{ $tables{$bc_table}[0] } ;
							$i++
						  )
						{
							my $bc_value =
							  $bc_scale * ${ $tables{$bc_table}[1] }[$i];
							if ( $bc_field =~ m/pressure/ ) {
								${ ${ $bc_pumps{$bc_component} }[$i] }[1] =
								  $bc_value;
							}
							else {
								${ ${ $bc_pumps{$bc_component} }[$i] }[2] =
								  $bc_value;
							}
						}
					}
					else {
						my @arr = ();
						my $ilast = $#{ $tables{$bc_table}[0] };
						#if($bc_table =~m/1/) { die "bc_table is $bc_table : $ilast\n"; }
						for (
							my $i = 0 ;
							$i <= $ilast ;
							$i++
						  )
						{
							my $bc_value =
							  $bc_scale * ${ $tables{$bc_table}[1] }[$i];
							if ( $bc_field =~ m/pressure/ ) {
								push(
									@arr,
									[
										${ ${ $tables{$bc_table} }[0] }[$i],
										$bc_value, 0
									]
								);
							}
							else {
								push(
									@arr,
									[
										${ ${ $tables{$bc_table} }[0] }[$i],
										0, $bc_value
									]
								);
							}
						}
						$bc_pumps{$bc_component} = \@arr;
					}
				}
			}
			$bc_on    = 0;
			$bc_scale = 1;
		}
		if (/\<BoundaryCondition\s+/) {
			$bc_on = 1;
		}
		if ( $bc_on == 1 ) {
			my $line = $_;
			my @arr1 = split( /\"/, $line );
			if (/fieldname/) {
				my $field = $arr1[1];
				$field =~ s/\s+//g;
				chomp($field);
				$bc_field = lc($field);
				if ( $bc_field eq "velocity" ) {
					$bc_type = 1;
				}
				elsif ( $bc_field eq "pressure" ) {
					$bc_type = 2;
				}
				else {
					$bc_type = 3;
				}
			}
			elsif (/setname/) {
				$bc_nset = $arr1[1];
			}
			elsif (/direction/) {
				my @arr = split( /\s+/, $arr1[1] );
				if ( $arr[1] != 0 ) {
					$bc_scale *= $arr[1];
					$bc_component = 2;
				}
				elsif ( $arr[0] != 0 ) {
					$bc_scale *= $arr[0];
					$bc_component = 1;
				}
			}
			elsif (/scale/) {
				$bc_scale *= $arr1[1];
			}
			elsif (/timetable/) {
				$bc_table = $arr1[1];
			}
			elsif (/component/) {
				$bc_component = $arr1[1];
			}
		}
	}
	close IN;
}
foreach my $key ( sort( keys %bc_pumps ) ) {
	print "begin_pump $key\n";
	my $ilast = $#{ $bc_pumps{$key} };
	for ( my $i = 0 ; $i <= $ilast ; $i++ ) {
		my $time         = ${ ${ $bc_pumps{$key} }[$i] }[0];
		my $pumppressure = ${ ${ $bc_pumps{$key} }[$i] }[1];
		my $pumpflow     = ${ ${ $bc_pumps{$key} }[$i] }[2];
		print "        add $time $pumppressure $pumpflow\n";
	}
	print "end_pump $key\n";
}
print "END\n";
