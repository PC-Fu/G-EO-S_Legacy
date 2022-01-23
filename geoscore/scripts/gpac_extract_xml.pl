#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die
"usage: $0 <parent node name> <filename of bash script used to generate constitutive models>\n"
  if ( $narg != 2 );

my $pname = shift;
my $ff    = shift;

my $suffix  = "Type";
my $suffix2 = "";
if ( $pname =~ m/Material/ ) {
	$suffix  = "MaterialType";
	$suffix2 = "Material";
}

#############################################
#FILL DATA STRUCTURES

my @leaves    = ();
my %bases     = ();
my %variables = ();

die "cannot open $ff\n" if ( !open( IN, "<$ff" ) );
while (<IN>) {
	my $line = $_;
	$line =~ s/^\s+//;
	chomp($line);
	if (/gpac_constitutive_model/) {
		my @arr = split( /\s+/, $line );
		my $arg0 = -1;
		for ( my $i = 0 ; $i <= $#arr ; $i++ ) {
			if ( $arr[$i] =~ m/gpac_constitutive_model/ ) {
				$arg0 = $i;
				last;
			}
		}
		die "could not find the first argument\n" if ( $arg0 < 0 );

		#GET CLASS AND BASE NAMES
		my $className = $arr[ $arg0 + 1 ];
		my $baseName  = $arr[ $arg0 + 2 ];

		#SET CLASS TYPE - SHOULD REFLECT SAME IN gpac_constitutive_model.pl
		my $classType = "leaf";
		if ( $baseName =~ m/ConstitutiveBase/ ) {
			$classType = "base";
		}
		elsif ( $className =~ m/Intermediate/ ) {
			$classType = "intermediate";
		}
		else {
			push( @leaves, $className );
		}

		#GET ALL REAL VARIABLES
		my @arr0 = split( /_r0\ /, $line );
		my @vars = ();
		for ( my $i = 1 ; $i <= $#arr0 ; $i++ ) {
			my @arrtmp = split( /\s+/, $arr0[$i] );
			my $name = $arrtmp[1];
			push( @vars, $name );
		}

		#CACHE IN HASH TABLES
		$variables{$className} = \@vars;
		$bases{$className}     = $baseName;
	}
}
close IN;

#Put in stub for constitutive base
{
	my @vars = ();
	$variables{"ConstitutiveBase"} = \@vars;
}

#############################################
#WRITE OUT LEAF CLASS XML

for ( my $i = 0 ; $i <= $#leaves ; $i++ ) {
	my $className = $leaves[$i];
	my $typename  = $className . $suffix;
	print " <xsd:complexType name=\"$typename\">
";
	&PrintVariables( $className, $variables{$className} );
	if ( exists $bases{$className} ) {
		my $cn = $bases{$className};
		&PrintVariables( $cn, $variables{$cn} );
		while ( exists $bases{$cn} ) {
			$cn = $bases{$cn};
			&PrintVariables( $cn, $variables{$cn} );
		}
	}
	print " </xsd:complexType>\n\n";
}

#############################################
#WRITE OUT PARENT CLASS XML

print " <xsd:complexType name=\"$pname\">
  <xsd:choice>
";
for ( my $i = 0 ; $i <= $#leaves ; $i++ ) {
	my $n = $leaves[$i] . $suffix2;
	my $t = $leaves[$i] . $suffix;
	print "   <xsd:element name=\"$n\" type=\"$t\"/>\n";
}
print "  </xsd:choice>
";

if ( $pname =~ m/Contact/ ) {
	print "  <!--tolerance-->
  <xsd:attribute name=\"penetrationTol\" type=\"xsd:double\"/>
  <xsd:attribute name=\"spatialTol\" type=\"xsd:double\"/>
  <xsd:attribute name=\"areaTol\" type=\"xsd:double\"/>
  <xsd:attribute name=\"cosMinTol\" type=\"xsd:double\"/>
  <xsd:attribute name=\"searchRadiusFactor\" type=\"xsd:double\"/>
  <xsd:attribute name=\"searchRadiusVelocityFactor\" type=\"xsd:double\"/>
  <xsd:attribute name=\"feParentSolnTol\" type=\"xsd:double\"/>
  <xsd:attribute name=\"maximumSeparation\" type=\"unittype\"/>
  <!--options-->
  <xsd:attribute name=\"active\" type=\"xsd:integer\"/>
  <xsd:attribute name=\"smoothed\" type=\"xsd:boolean\"/>
  <xsd:attribute name=\"allowSelfContact\" type=\"xsd:boolean\"/>\n";
}
else {
	print
	  "  <xsd:attribute name=\"name\" type=\"xsd:string\" use=\"required\"/>\n";
}
print " </xsd:complexType>\n\n";

{
	my $n = $pname;
	$n =~ s/Type/sType/;
	my $nn = $pname;
	$nn =~ s/Type//;
	print " <xsd:complexType name=\"$n\">
  <xsd:sequence>
   <xsd:element name=\"$nn\" type=\"$pname\" maxOccurs=\"unbounded\"/>
  </xsd:sequence>
 </xsd:complexType>
 ";
}

sub PrintVariables {
	my $className = shift;
	my $varptr    = shift;
	#print "=====>$className\n";
	my @vars      = @$varptr;
	print "   <!-- $className -->\n";
	for ( my $i = 0 ; $i <= $#vars ; $i++ ) {
		print "   <xsd:attribute name=\"$vars[$i]\" type=\"unittype\"/>\n";
	}
}
