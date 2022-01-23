#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die
"usage: $0 <class name> <filename of bash script used to generate constitutive models>\n"
  if ( $narg != 2 );

my $pname = shift;
my $ff    = shift;

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


#############################################
#WRITE OUT HEADER FILE
{
	my $ns = "Interfaces";
	my $suffix = "Cnt";
	my $baseName = "InterfaceBase";
	if($pname=~m/Material/)
	{
		$ns = "Materials";
		$suffix = "Mat";
		$baseName = "MaterialBase";
	}
	elsif($pname=~m/CohesiveZone/)
	{
		$ns = "CohesiveZones";
		$suffix = "CZ";
		$baseName = "CohesiveZoneBase";		
	}
	my $hf = $pname . ".h";
	my $type = $pname . "Type";
	my $vv = uc($pname);
	$vv .= "_H_";
	die "cannot open header file: $hf\n" if(!open(OUT,">$hf"));
	my $heading = &Heading($pname);

	print OUT "$heading

#ifndef $vv
#define $vv

";
for(my $i = 0; $i <= $#leaves; $i++)
{
	print OUT "#include \"$leaves[$i].h\"\n";
}
print OUT "\nnamespace $ns
{
  enum $type
  {
";
	for(my $i = 0; $i <= $#leaves; $i++)
	{
		print OUT "    $leaves[$i]$suffix,\n";
	}
	print OUT"    NONE,
  };

  void ReadXML( TICPP::HierarchicalDataNode& node,
                std::unique_ptr<$baseName>& cz );


  void Allocate( std::unique_ptr<$baseName>& cz,
                 const $type type );

  $type StringToModelType( const std::string& name );
  std::string ModelTypeToString( const $type type );


}

#endif /* $vv */	
";
	close OUT;
}

#############################################
#WRITE OUT CLASS FILE
{
	my $ns = "Interfaces";
	my $suffix = "Cnt";
	my $baseName = "InterfaceBase";
	my $str_suffix = "";
	if($pname=~m/Material/)
	{
		$ns = "Materials";
		$suffix = "Mat";
		$baseName = "MaterialBase";
		$str_suffix = "Material";
	}
	elsif($pname=~m/CohesiveZone/)
	{
		$ns = "CohesiveZones";
		$suffix = "CZ";
		$baseName = "CohesiveZoneBase";		
		$str_suffix = "";
	}
	my $hf = $pname . ".cpp";
	my $type = $pname . "Type";
	my $vv = uc($pname);
	$vv .= "_H_";
	die "cannot open header file: $hf\n" if(!open(OUT,">$hf"));
	my $heading = &Heading($pname);

	print OUT "$heading

#include \"$pname.h\"
#include \"IO/ticpp/HierarchicalDataNode.h\"
#include \"DataStructures/EncapsulatedObjects/EncapsulatedObjectBase.h\"

namespace $ns
{

  void ReadXML( TICPP::HierarchicalDataNode& node,
                std::unique_ptr<$baseName>& mat )
  {

    TICPP::HierarchicalDataNode* childNode = node.Next(true);
    if(!childNode)
      throw GPException(\"No material defined for the element region!\");

    $type type = StringToModelType( childNode->Heading() );

    Allocate( mat, type );

    if( mat )
      mat->ReadXML( *childNode );

  }



  void Allocate( std::unique_ptr<$baseName>& mat,
                 const $type type )
  {
    switch( type )
    {
";
	for(my $i = 0; $i <= $#leaves; $i++)
	{
		print OUT "      case $leaves[$i]$suffix:
        mat = std::unique_ptr<$baseName>( new $leaves[$i] );
        break;
";
	}
	print OUT "      default:
        break;
    }

  }


  $type StringToModelType( const std::string& name )
  {
    $type rval;

";
	for(my $i = 0; $i <= $#leaves; $i++)
	{
		my $elseif = "if";
		if($i > 0)
		{
			$elseif = "else if";
		}
		print OUT "    $elseif"."( !(name.compare(\"$leaves[$i]$str_suffix\")) )
    {
      rval = $leaves[$i]$suffix;
    }
";
	}
	print OUT "    else
    {
      rval = NONE;
//      throw GPException(\"$pname"."::StringToModelType()\");
    }

    return rval;
  }

  std::string ModelTypeToString( const $type type )
  {
    std::string rval;

    switch( type )
    {
";
	for(my $i = 0; $i <= $#leaves; $i++)
	{
		print OUT "      case $leaves[$i]$suffix:
        rval = \"$leaves[$i]$str_suffix\";
        break;
";
	}
	print OUT "      default:
        throw GPException(\"$pname"."::ModelTypeToString()\");
    }

    return rval;
  }
}
";
	close OUT;
}

sub Heading {
	my $filename = shift;

	my $date = "";
	{
		die "cannot open date function\n" if ( !open( DATE, "date|" ) );
		while (<DATE>) {
			$date = $_;
		}
		close DATE;
		chomp($date);
		$date =~ s/^\s+//;
	}

	my $str = "/*
 * $filename
 *
 *  Created on: $date
 *      Author: johnson346, settgast
 */
 ";
	return $str;
}
