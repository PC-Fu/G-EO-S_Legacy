#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die
"usage: $1 <class_name> <base_name> <base_type> <variables (2-tuples): type varname ...>\n"
  if ( $narg < 3 );

#get class and base types
my $className = shift;
my $baseName  = shift;
my $baseType  = shift;

my $classType = "leaf";
if ( $baseName =~ m/ConstitutiveBase/ ) {
	$classType = "base";
}
elsif ( $className =~ m/Intermediate/ ) {
	$classType = "intermediate";
}

###################################
# GET VARIABLE NAMES
###################################

my @s_intVars = ();
my @s_r0Vars  = ();
my @s_r1Vars  = ();
my @s_r2Vars  = ();
my @s_r2sVars = ();

my @p_intVars = ();
my @p_r0Vars  = ();
my @p_r1Vars  = ();
my @p_r2Vars  = ();
my @p_r2sVars = ();

my @s_intAliases = ();
my @s_r0Aliases  = ();
my @s_r1Aliases  = ();
my @s_r2Aliases  = ();
my @s_r2sAliases = ();

my @p_intAliases = ();
my @p_r0Aliases  = ();
my @p_r1Aliases  = ();
my @p_r2Aliases  = ();
my @p_r2sAliases = ();

#INCLUDE FILES
my @p_public    = ();
my @p_protected = ();
my @p_private   = ();

my @s_public    = ();
my @s_protected = ();
my @s_private   = ();

my @b_public    = ();
my @b_protected = ();
my @b_private   = ();

my %p_tanVars = ();
my %s_tanVars = ();

#PARSE ARGUMENTS
for ( my $i = 0 ; $i < ( $narg - 3 ) ; $i += 3 ) {
	my $type     = shift;
	my $varName  = shift;
	my $varAlias = shift;

	my $lctype = lc($type);
	if ( $lctype =~ m/^s_int/ ) {
		push( @s_intVars,    $varName );
		push( @s_intAliases, $varName );
	}
	elsif ( $lctype =~ m/^s_r0/ ) {
		push( @s_r0Vars,    $varName );
		push( @s_r0Aliases, $varAlias );
	}
	elsif ( $lctype =~ m/^s_r1/ ) {
		push( @s_r1Vars,    $varName );
		push( @s_r1Aliases, $varAlias );
	}
	elsif ( $lctype =~ m/^s_tan/ ) {
		push( @s_r1Vars,    $varName );
		push( @s_r1Aliases, $varAlias );
		$s_tanVars{$varName} = 1;
	}
	elsif ( $lctype =~ m/^s_r2s/ ) {
		push( @s_r2sVars,    $varName );
		push( @s_r2sAliases, $varAlias );
	}
	elsif ( $lctype =~ m/^s_r2/ ) {
		push( @s_r2Vars,    $varName );
		push( @s_r2Aliases, $varAlias );
	}
	elsif ( $lctype =~ m/^p_int/ ) {
		push( @p_intVars,    $varName );
		push( @p_intAliases, $varAlias );
	}
	elsif ( $lctype =~ m/^p_r0/ ) {
		push( @p_r0Vars,    $varName );
		push( @p_r0Aliases, $varAlias );
	}
	elsif ( $lctype =~ m/^p_r1/ ) {
		push( @p_r1Vars,    $varName );
		push( @p_r1Aliases, $varAlias );
	}
	elsif ( $lctype =~ m/^p_tan/ ) {
		push( @p_r1Vars,    $varName );
		push( @p_r1Aliases, $varAlias );
		$p_tanVars{$varName} = 1;
	}
	elsif ( $lctype =~ m/^p_r2s/ ) {
		push( @p_r2sVars,    $varName );
		push( @p_r2sAliases, $varAlias );
	}

	elsif ( $lctype =~ m/^p_public/ ) {
		push( @p_public, $varName );
	}
	elsif ( $lctype =~ m/^p_protected/ ) {
		push( @p_protected, $varName );
	}
	elsif ( $lctype =~ m/^p_private/ ) {
		push( @p_private, $varName );
	}

	elsif ( $lctype =~ m/^s_public/ ) {
		push( @s_public, $varName );
	}
	elsif ( $lctype =~ m/^s_protected/ ) {
		push( @s_protected, $varName );
	}
	elsif ( $lctype =~ m/^s_private/ ) {
		push( @s_private, $varName );
	}

	elsif ( $lctype =~ m/^b_public/ ) {
		push( @b_public, $varName );
	}
	elsif ( $lctype =~ m/^b_protected/ ) {
		push( @b_protected, $varName );
	}
	elsif ( $lctype =~ m/^b_private/ ) {
		push( @b_private, $varName );
	}

	else {
		die "cannot determine type: $type\n";
	}
}

my @p_function_files = ( \@p_public, \@p_protected, \@p_private );
my @s_function_files = ( \@s_public, \@s_protected, \@s_private );
my @b_function_files = ( \@b_public, \@b_protected, \@b_private );
my @function_files =
  ( \@p_function_files, \@s_function_files, \@b_function_files );

&PrintHeaderFile(
	$classType,  $className, $baseName,  $baseType, 
	\@s_intVars, \@s_r0Vars, \@s_r1Vars, \@s_r2Vars, \@s_r2sVars,
	\@p_intVars, \@p_r0Vars, \@p_r1Vars,  \@p_r2Vars, \@p_r2sVars,
	\@function_files, \%s_tanVars, \%p_tanVars,
	\@s_intAliases, \@s_r0Aliases, \@s_r1Aliases,
	\@s_r2Aliases, \@s_r2sAliases,
	\@p_intAliases, \@p_r0Aliases, \@p_r1Aliases,
	\@p_r2Aliases, \@p_r2sAliases
);

&PrintClassFile($classType, $className, $baseName, $baseType, \@function_files);

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

###################################
# PRINT HEADER FILE
###################################

sub PrintHeaderFile {
	my $classType = shift;
	my $className = shift;
	my $baseName  = shift;
	my $baseType  = shift;

	my $ucClassNameH = uc($className) . "_H_";

	my $s_intVars = shift;
	my $s_r0Vars  = shift;
	my $s_r1Vars  = shift;
	my $s_r2Vars  = shift;
	my $s_r2sVars = shift;

	my $p_intVars = shift;
	my $p_r0Vars  = shift;
	my $p_r1Vars  = shift;
	my $p_r2Vars  = shift;
	my $p_r2sVars = shift;

	my $ffiles = shift;

	my $s_tanVars = shift;
	my $p_tanVars = shift;
	
	my $s_intAliases = shift;
	my $s_r0Aliases  = shift;
	my $s_r1Aliases  = shift;
	my $s_r2Aliases  = shift;
	my $s_r2sAliases = shift;

	my $p_intAliases = shift;
	my $p_r0Aliases  = shift;
	my $p_r1Aliases  = shift;
	my $p_r2Aliases  = shift;
	my $p_r2sAliases = shift;
	
	my $heading = &Heading( $className . ".h" );

	my $include = "#include \"Utilities/GeometryUtilities.h\"";
	my $baseType = "";
	if ( $classType =~ m/base/ ) {
		$include .= "\n#include \"Constitutive/ConstitutiveBase.h\"";
		my @arr_tmp = split(/Base/, $className);
		$baseType = $arr_tmp[0];
	}
	else {
		$include .= "\n#include \"$baseName" . ".h\"";
	}
	if( $className=~m/Dieterich/)
	{
		$include .= "
#include \"PhysicsSolvers/Seismicity/BoundaryElementDataManagerT.h\"
#include \"DataStructures/Tables/Table.h\"
#if GPAC_MPI
#include <mpi.h>
#endif
";
	}
	elsif( $className=~m/GeodynMaterial/)
	{
		$include .= "
#include \"libmmlib.h\"
";
	}
	
	open( OUT, ">$className" . ".h" );
	my $disclaimer = &Disclaimer();
	print OUT "$disclaimer
#ifndef $ucClassNameH
#define $ucClassNameH

$include

$heading

";

	#print parameter class declaration
	{
		my $str = &DataClassDeclaration(
			$classType,                  $className . "ParameterData",
			$baseName . "ParameterData", $p_intVars,
			$p_r0Vars,                   $p_r1Vars,
			$p_r2Vars,                   $p_r2sVars,
			$$ffiles[0],                 $p_tanVars,
			$p_intAliases,               $p_r0Aliases,
			$p_r1Aliases,                $p_r2Aliases,
			$p_r2sAliases);
		print OUT $str;
	}

	#print state class declaration
	{
		my $str = &DataClassDeclaration(
			$classType,              $className . "StateData",
			$baseName . "StateData", $s_intVars,
			$s_r0Vars,               $s_r1Vars,
			$s_r2Vars,               $s_r2sVars,
			$$ffiles[1],             $s_tanVars,
			$s_intAliases,           $s_r0Aliases,
			$s_r1Aliases,            $s_r2Aliases,
            $s_r2sAliases);
		print OUT $str;
	}

	#print algorithmic class declaration
	if ( $classType =~ m/leaf/ ) {
		my $str = &ClassDeclarationLeaf( $className, $baseName, $$ffiles[2] );
		print OUT $str;
	}
	else {
		my $str =
		  &ClassDeclaration( $classType, $className, $baseName, $$ffiles[2], $baseType );
		print OUT $str;
	}

	print OUT "#endif /* $ucClassNameH */\n";
	close OUT;
}

#___________________________________________

sub DataClassDeclaration {
	my $classType = shift;
	my $className = shift;
	my $baseName  = shift;
	my $intVars   = shift;
	my $r0Vars    = shift;
	my $r1Vars    = shift;
	my $r2Vars    = shift;
	my $r2sVars   = shift;

	my $filecats = shift;

	my $public_str    = &GetFunctionDeclarations( $$filecats[0] );
	my $protected_str = &GetFunctionDeclarations( $$filecats[1] );
	my $private_str   = &GetFunctionDeclarations( $$filecats[2] );
	if ( length($protected_str) > 0 ) {
		my $tstr = "protected:\n$protected_str";
		$protected_str = $tstr;
	}
	if ( length($private_str) > 0 ) {
		my $tstr = "private:\n$private_str";
		$private_str = $tstr;
	}

	my $tangent_variables = shift;
	
	my $intAliases   = shift;
	my $r0Aliases    = shift;
	my $r1Aliases    = shift;
	my $r2Aliases    = shift;
	my $r2sAliases   = shift;

	my $inherit_str = "";
	if ( $classType =~ m/base/ ) {

	}
	else {
		$inherit_str = ": public $baseName";
	}

	my $str = "
//**********************************************************************************************************************
//**********************************************************************************************************************


class $className $inherit_str
{

public:\n";

	if ( $classType =~ m/base/ ) {

	}
	else {
		$str .= "\n  typedef $baseName base;\n";
	}

	$str .= &Constructors(
		$classType, $className, $intVars, $r0Vars,
		$r1Vars,    $r2Vars,    $r2sVars
	);
	$str .= "\n$public_str";

	$str .=
	  &VariableCountsAndNames( $classType, $intVars, $r0Vars, $r1Vars, $r2Vars,
		$r2sVars, $intAliases, $r0Aliases, $r1Aliases, $r2Aliases,
		$r2sAliases );
	$str .= "\n";

	$str .= &SerializeDeserialize(
		$classType, $className, $intVars, $r0Vars,
		$r1Vars,    $r2Vars,    $r2sVars
	);

	$str .= &Operators(
		$classType, $className, $baseName, $intVars,
		$r0Vars,    $r1Vars,    $r2Vars,   $r2sVars
	);

	$str .= &MapFunctions(
		$classType, $className, $baseName,
		$intVars,   $r0Vars,    $r1Vars,
		$r2Vars,    $r2sVars,   $tangent_variables
	);
	
	if($className=~/StateData/)
	{
	}
	else
	{
		$str .= "  virtual void ReadXML( const TICPP::HierarchicalDataNode& node )
  {\n";
		if ( $classType =~ m/base/ ) {
		}
		else {
			$str .= "    $baseName" . "::ReadXML( node );\n";
                }
      		for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
				$str .= "    $$r0Vars[$i] = node.GetAttributeOrDefault(\"$$r0Aliases[$i]\", 0.0);\n";
      		}
		if ( $classType =~ m/leaf/ ) {
			$str .= "    PostReadXML( node );\n";
		}
		$str .= "
  }
";
		if ( $classType =~ m/base/ ) {
			$str .= "  virtual void PostReadXML( const TICPP::HierarchicalDataNode& ) {}\n"
		}
	}
	$str .= "\n$protected_str\n$private_str\n";

	$str .= "};\n";
}

#___________________________________________

sub Constructors {
	my $classType = shift;
	my $className = shift;
	my $intVars   = shift;
	my $r0Vars    = shift;
	my $r1Vars    = shift;
	my $r2Vars    = shift;
	my $r2sVars   = shift;

	my $classOnly = $className;
	$classOnly =~ s/StateData//;
	$classOnly =~ s/ParameterData//;

	my ( $varDec, $varInit, $varInitCopy ) =
	  &VariableDeclarations( $classType, $intVars, $r0Vars, $r1Vars, $r2Vars,
		$r2sVars );
	my $str = "$varDec

  $className()$varInit  {}

  $className( const $className" . "& source)$varInitCopy  {}

  ~$className() {}\n";

	if ( $classType =~ m/base/ ) {
	}
	else {
		$str .= "  friend class ConstitutiveBase;
  friend class $classOnly;\n";
	}
	return $str;
}

sub VariableDeclarations {
	my $classType = shift;
	my $intVars   = shift;
	my $r0Vars    = shift;
	my $r1Vars    = shift;
	my $r2Vars    = shift;
	my $r2sVars   = shift;

	my $varDec      = "";
	my $varInit     = "";
	my $varInitCopy = "";

	my $ok = 0;
	if ( $classType =~ m/base/ ) {
		if (   @$intVars > 0
			|| @$r0Vars > 0
			|| @$r1Vars > 0
			|| @$r2Vars > 0
			|| @$r2sVars > 0 )
		{
			$varInit     .= ":";
			$varInitCopy .= ":";
		}
	}
	else {
		$varInit     .= ":\n    base()";
		$varInitCopy .= ":\n    base( source )";
		if (   @$intVars > 0
			|| @$r0Vars > 0
			|| @$r1Vars > 0
			|| @$r2Vars > 0
			|| @$r2sVars > 0 )
		{
			$varInit     .= ",";
			$varInitCopy .= ",";
		}
	}
	$varInit     .= "\n";
	$varInitCopy .= "\n";

	#INTEGERS
	my $type = "int";
	for ( my $i = 0 ; $i < @$intVars ; $i++ ) {
		$varDec      .= "  $type $$intVars[$i];\n";
		$varInit     .= "    $$intVars[$i]" . "(0),\n";
		$varInitCopy .= "    $$intVars[$i](source.$$intVars[$i]),\n";
	}

	#REALS
	$type = "realT";
	for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
		$varDec      .= "  $type $$r0Vars[$i];\n";
		$varInit     .= "    $$r0Vars[$i]" . "(0),\n";
		$varInitCopy .= "    $$r0Vars[$i](source.$$r0Vars[$i]),\n";
	}

	#R1's
	$type = "R1Tensor";
	for ( my $i = 0 ; $i < @$r1Vars ; $i++ ) {
		$varDec      .= "  $type $$r1Vars[$i];\n";
		$varInit     .= "    $$r1Vars[$i]" . "(0),\n";
		$varInitCopy .= "    $$r1Vars[$i](source.$$r1Vars[$i]),\n";
	}

	#R2's
	$type = "R2Tensor";
	for ( my $i = 0 ; $i < @$r2Vars ; $i++ ) {
		$varDec      .= "  $type $$r2Vars[$i];\n";
		#$varInit     .= "    $$r2Vars[$i]" . "(0),\n";
		$varInitCopy .= "    $$r2Vars[$i](source.$$r2Vars[$i]),\n";
	}

	#R2S's
	$type = "R2SymTensor";
	for ( my $i = 0 ; $i < @$r2sVars ; $i++ ) {
		$varDec      .= "  $type $$r2sVars[$i];\n";
		$varInit     .= "    $$r2sVars[$i]" . "(0),\n";
		$varInitCopy .= "    $$r2sVars[$i](source.$$r2sVars[$i]),\n";
	}

	if (   @$intVars > 0
		|| @$r0Vars > 0
		|| @$r1Vars > 0
		|| @$r2Vars > 0
		|| @$r2sVars > 0 )
	{
		my $tmp = substr( $varInit, 0, -2 );
		$varInit     = "$tmp\n";
		$tmp         = substr( $varInitCopy, 0, -2 );
		$varInitCopy = "$tmp\n";
		return ( $varDec, $varInit, $varInitCopy );
	}
	else {
		return ( $varDec, $varInit, $varInitCopy );
	}
}

#___________________________________________
sub VariableCountsAndNames {
	my $classType = shift;
	
	my $intVars   = shift;
	my $r0Vars    = shift;
	my $r1Vars    = shift;
	my $r2Vars    = shift;
	my $r2sVars   = shift;
	
	my $intAliases   = shift;
	my $r0Aliases    = shift;
	my $r1Aliases    = shift;
	my $r2Aliases    = shift;
	my $r2sAliases   = shift;

	my $op       = "=";
	my $str = "  static void GetVariableCounts( localIndex&";
	if ( $classType =~ m/base/ ) 
	{
		if(@$intVars > 0)
		{
			$str .= "intVarCounts";
		}
		$str .= ",\n                                 localIndex& ";
		if ( @$r0Vars > 0 ) {
			$str .= "realVarCounts";
		}
		$str .= ",\n                                 localIndex& ";
		if ( @$r1Vars > 0 ) {
			$str .= "R1TensorVarCounts";
		}
		$str .= ",\n                                 localIndex& ";
		if ( @$r2Vars > 0 ) {
			$str .= "R2TensorVarCounts";
		}
		$str .= ",\n                                 localIndex& ";
		if ( @$r2sVars > 0 ) {
			$str .= "R2SymTensorVarCounts";
		}
		$str .= " )\n";
		$str .= "  {
";
	}
	else 
	{
		$op       = "+=";
		$str .= " intVarCounts,
                                 localIndex& realVarCounts,
                                 localIndex& R1TensorVarCounts,
                                 localIndex& R2TensorVarCounts,
                                 localIndex& R2SymTensorVarCounts )
  {
    base::GetVariableCounts( intVarCounts,
                             realVarCounts,
                             R1TensorVarCounts,
                             R2TensorVarCounts,
                             R2SymTensorVarCounts );
";
	}
	if ( @$intVars > 0 ) {
		$str .= "    intVarCounts $op " . (@$intVars) . ";\n";
	}
	if ( @$r0Vars > 0 ) {
		$str .= "    realVarCounts $op " . (@$r0Vars) . ";\n";
	}
	if ( @$r1Vars > 0 ) {
		$str .= "    R1TensorVarCounts $op " . (@$r1Vars) . ";\n";
	}
	if ( @$r2Vars > 0 ) {
		$str .= "    R2TensorVarCounts $op " . (@$r2Vars) . ";\n";
	}
	if ( @$r2sVars > 0 ) {
		$str .= "    R2SymTensorVarCounts $op " . (@$r2sVars) . ";\n";
	}
	$str .= "\n  }

  static void GetVariableNames( sArray1d&";
  
	if ( $classType =~ m/base/ ) 
	{
		if(@$intVars > 0)
		{
			$str .= "intNames";
		}
		$str .= ",\n                                 sArray1d& ";
		if ( @$r0Vars > 0 ) {
			$str .= "realNames";
		}
		$str .= ",\n                                 sArray1d& ";
		if ( @$r1Vars > 0 ) {
			$str .= "R1TensorNames";
		}
		$str .= ",\n                                 sArray1d& ";
		if ( @$r2Vars > 0 ) {
			$str .= "R2TensorNames";
		}
		$str .= ",\n                                 sArray1d& ";
		if ( @$r2sVars > 0 ) {
			$str .= "R2SymTensorNames";
		}
		$str .= " )\n";
		$str .= "  {
";
	}
	else 
	{
		$str .= " intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
";
	}
	for ( my $i = 0 ; $i < @$intVars ; $i++ ) {
		$str .= "    intNames.push_back(\"$$intAliases[$i]\");\n";
	}
	for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
		$str .= "    realNames.push_back(\"$$r0Aliases[$i]\");\n";
	}
	for ( my $i = 0 ; $i < @$r1Vars ; $i++ ) {
		$str .= "    R1TensorNames.push_back(\"$$r1Aliases[$i]\");\n";
	}
	for ( my $i = 0 ; $i < @$r2Vars ; $i++ ) {
		$str .= "    R2TensorNames.push_back(\"$$r2Aliases[$i]\");\n";
	}
	for ( my $i = 0 ; $i < @$r2sVars ; $i++ ) {
		$str .= "    R2SymTensorNames.push_back(\"$$r2sAliases[$i]\");\n";
	}
	$str .= "  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>&";
  
	if ( $classType =~ m/base/ ) 
	{
		if(@$intVars > 0)
		{
			$str .= "intOffsets";
		}
		$str .= ",\n                                 std::map<std::string, size_t>& ";
		if ( @$r0Vars > 0 ) {
			$str .= "realOffsets";
		}
		$str .= ",\n                                 std::map<std::string, size_t>& ";
		if ( @$r1Vars > 0 ) {
			$str .= "R1TensorOffsets";
		}
		$str .= ",\n                                 std::map<std::string, size_t>& ";
		if ( @$r2Vars > 0 ) {
			$str .= "R2TensorOffsets";
		}
		$str .= ",\n                                 std::map<std::string, size_t>& ";
		if ( @$r2sVars > 0 ) {
			$str .= "R2SymTensorOffsets";
		}
		$str .= " ) const\n";
		$str .= "  {
";
	}
	else 
	{
		$str .= " intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
";
	}
	for ( my $i = 0 ; $i < @$intVars ; $i++ ) {
		$str .= "    intOffsets[\"$$intAliases[$i]\"] = (char*)(&"."$$intVars[$i]) - (char*)this;\n";
	}
	for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
		$str .= "    realOffsets[\"$$r0Aliases[$i]\"] = (char*)(&"."$$r0Vars[$i]) - (char*)this;\n";
	}
	for ( my $i = 0 ; $i < @$r1Vars ; $i++ ) {
		$str .= "    R1TensorOffsets[\"$$r1Aliases[$i]\"] = (char*)(&"."$$r1Vars[$i]) - (char*)this;\n";
	}
	for ( my $i = 0 ; $i < @$r2Vars ; $i++ ) {
		$str .= "    R2TensorOffsets[\"$$r2Aliases[$i]\"] = (char*)(&"."$$r2Vars[$i]) - (char*)this;\n";
	}
	for ( my $i = 0 ; $i < @$r2sVars ; $i++ ) {
		$str .= "    R2SymTensorOffsets[\"$$r2sAliases[$i]\"] = (char*)(&"."$$r2sVars[$i]) - (char*)this;\n";
	}
	$str .= "  }

  virtual void GetVariableValues( std::map<std::string, int>&";
  
	if ( $classType =~ m/base/ ) 
	{
		if(@$intVars > 0)
		{
			$str .= "intValues";
		}
		$str .= ",\n                                 std::map<std::string, realT>& ";
		if ( @$r0Vars > 0 ) {
			$str .= "realValues";
		}
		$str .= ",\n                                 std::map<std::string, R1Tensor>& ";
		if ( @$r1Vars > 0 ) {
			$str .= "R1TensorValues";
		}
		$str .= ",\n                                 std::map<std::string, R2Tensor>& ";
		if ( @$r2Vars > 0 ) {
			$str .= "R2TensorValues";
		}
		$str .= ",\n                                 std::map<std::string, R2SymTensor>& ";
		if ( @$r2sVars > 0 ) {
			$str .= "R2SymTensorValues";
		}
		$str .= " )\n";
		$str .= "  {
";
	}
	else 
	{
		$str .= " intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
";
	}
	for ( my $i = 0 ; $i < @$intVars ; $i++ ) {
		$str .= "    intValues[\"$$intAliases[$i]\"] = $$intVars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
		$str .= "    realValues[\"$$r0Aliases[$i]\"] = $$r0Vars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r1Vars ; $i++ ) {
		$str .= "    R1TensorValues[\"$$r1Aliases[$i]\"] = $$r1Vars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r2Vars ; $i++ ) {
		$str .= "    R2TensorValues[\"$$r2Aliases[$i]\"] = $$r2Vars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r2sVars ; $i++ ) {
		$str .= "    R2SymTensorValues[\"$$r2sAliases[$i]\"] = $$r2sVars[$i];\n";
	}
	$str .= "  }\n";
	
	return $str;
}

#___________________________________________

sub SerializeDeserialize {
	my $classType = shift;
	my $className = shift;
	my $intVars = shift;
	my $r0Vars  = shift;
	my $r1Vars  = shift;
	my $r2Vars  = shift;
	my $r2sVars = shift;
	
	my $signature = "const localIndex index,";
	my $val0      = "0";
	my $incr      = "++";
	my $ivar      = "index";
	if ( $className =~ m/StateData/ ) {
		$signature .= "\n                  const unsigned int stride,
                  const localIndex elemNum,";
		$val0 = "index";
		$incr = " += stride";
		$ivar = "elemNum";
	}
	my $signature2 = "";
	if( $classType =~ m/base/) 
	{
		$signature2 = ",\n                  localIndex& ";
        if ( @$intVars > 0 ) {
			$signature2 .= "intVarCounts";
		}
		$signature2 .= ",\n                  localIndex& ";
		if ( @$r0Vars > 0 ) {
			$signature2 .= "realVarCounts";
		}
		$signature2 .= ",\n                  localIndex& ";
		if ( @$r1Vars > 0 ) {
			$signature2 .= "R1TensorVarCounts";
		}
		$signature2 .= ",\n                  localIndex& ";
		if ( @$r2Vars > 0 ) {
			$signature2 .= "R2TensorVarCounts";
		}
		$signature2 .= ",\n                  localIndex& ";
		if ( @$r2sVars > 0 ) {
			$signature2 .= "R2SymTensorVarCounts";
		}
	}
	else
	{	
		$signature2 = ",
                  localIndex& intVarCounts,
                  localIndex& realVarCounts,
                  localIndex& R1TensorVarCounts,
                  localIndex& R2TensorVarCounts,
                  localIndex& R2SymTensorVarCounts";
	}

	my $str = "  void Serialize($signature
                  Array1dT<iArray1d*>& ";
	if ( $classType =~ m/base/ ) {
		if ( @$intVars > 0 ) {
			$str .= "intVars";
		}
		$str .= ",\n                  Array1dT<rArray1d*>& ";
		if ( @$r0Vars > 0 ) {
			$str .= "realVars";
		}
		$str .= ",\n                  Array1dT<Array1dT<R1Tensor>*>& ";
		if ( @$r1Vars > 0 ) {
			$str .= "R1Vars";
		}
		$str .= ",\n                  Array1dT<Array1dT<R2Tensor>*>& ";
		if ( @$r2Vars > 0 ) {
			$str .= "R2Vars";
		}
		$str .= ",\n                  Array1dT<Array1dT<R2SymTensor>*>& ";
		if ( @$r2sVars > 0 ) {
			$str .= "R2SymVars";
		}		
		$str .= "$signature2  ) const
  {
";	
#		if ( @$intVars > 0 ) {
#			$str .= "    intVarCounts = $val0;\n";
#		}
#		if ( @$r0Vars > 0 ) {
#			$str .= "    realVarCounts = $val0;\n";
#		}
#		if ( @$r1Vars > 0 ) {
#			$str .= "    R1TensorVarCounts = $val0;\n";
#		}
#		if ( @$r2Vars > 0 ) {
#			$str .= "    R2TensorVarCounts = $val0;\n";
#		}
#		if ( @$r2sVars > 0 ) {
#			$str .= "    R2SymTensorVarCounts = $val0;\n";
#		}
	}
	else {
		$str .= "intVars,
                  Array1dT<rArray1d*>& realVars,
                  Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                  Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                  Array1dT<Array1dT<R2SymTensor>*>& R2SymVars$signature2  ) const
  {
";
		my $arg_tmp = "index";
		if ( $className =~ m/StateData/ ) {
			$arg_tmp .= ", stride, elemNum";
		}
		$str .=
"    base::Serialize($arg_tmp, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                    intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );\n";
	}

	for ( my $i = 0 ; $i < @$intVars ; $i++ ) {
		$str .=
		  "    (*(intVars[intVarCounts]))[$ivar] = $$intVars[$i]; intVarCounts"
		  . "$incr;\n";
	}
	for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
		$str .=
"    (*(realVars[realVarCounts]))[$ivar] = $$r0Vars[$i]; realVarCounts"
		  . "$incr;\n";
	}
	for ( my $i = 0 ; $i < @$r1Vars ; $i++ ) {
		$str .=
"    (*(R1Vars[R1TensorVarCounts]))[$ivar] = $$r1Vars[$i]; R1TensorVarCounts"
		  . "$incr;\n";
	}
	for ( my $i = 0 ; $i < @$r2Vars ; $i++ ) {
		$str .=
"    (*(R2Vars[R2TensorVarCounts]))[$ivar] = $$r2Vars[$i]; R2TensorVarCounts"
		  . "$incr;\n";
	}
	for ( my $i = 0 ; $i < @$r2sVars ; $i++ ) {
		$str .=
"    (*(R2SymVars[R2SymTensorVarCounts]))[$ivar] = $$r2sVars[$i]; R2SymTensorVarCounts"
		  . "$incr;\n";
	}
	$str .= "  }


  void  Deserialize( $signature
                     const Array1dT<iArray1d*>& ";
	if ( $classType =~ m/base/ ) {
		if ( @$intVars > 0 ) {
			$str .= "intVars";
		}
		$str .= ",\n                  const Array1dT<rArray1d*>& ";
		if ( @$r0Vars > 0 ) {
			$str .= "realVars";
		}
		$str .= ",\n                  const Array1dT<Array1dT<R1Tensor>*>& ";
		if ( @$r1Vars > 0 ) {
			$str .= "R1Vars";
		}
		$str .= ",\n                  const Array1dT<Array1dT<R2Tensor>*>& ";
		if ( @$r2Vars > 0 ) {
			$str .= "R2Vars";
		}
		$str .= ",\n                  const Array1dT<Array1dT<R2SymTensor>*>& ";
		if ( @$r2sVars > 0 ) {
			$str .= "R2SymVars";
		}		
		$str .= "$signature2  )
  {
";	
#		if ( @$intVars > 0 ) {
#			$str .= "    intVarCounts = $val0;\n";
#		}
#		if ( @$r0Vars > 0 ) {
#			$str .= "    realVarCounts = $val0;\n";
#		}
#		if ( @$r1Vars > 0 ) {
#			$str .= "    R1TensorVarCounts = $val0;\n";
#		}
#		if ( @$r2Vars > 0 ) {
#			$str .= "    R2TensorVarCounts = $val0;\n";
#		}
#		if ( @$r2sVars > 0 ) {
#			$str .= "    R2SymTensorVarCounts = $val0;\n";
#		}
	}
	else {
		$str .= "intVars,
                     const Array1dT<rArray1d*>& realVars,
                     const Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                     const Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                     const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars$signature2 )
  {
";		my $arg_tmp = "index";
		if ( $className =~ m/StateData/ ) {
			$arg_tmp .= ", stride, elemNum";
		}
		$str .=
"    base::Deserialize($arg_tmp, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                    intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );\n";
	}
	for ( my $i = 0 ; $i < @$intVars ; $i++ ) {
		$str .=
		  "    $$intVars[$i] = (*(intVars[intVarCounts]))[$ivar]; intVarCounts"
		  . "$incr;\n";
	}
	for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
		$str .=
"    $$r0Vars[$i] = (*(realVars[realVarCounts]))[$ivar]; realVarCounts"
		  . "$incr;\n";
	}
	for ( my $i = 0 ; $i < @$r1Vars ; $i++ ) {
		$str .=
"    $$r1Vars[$i] = (*(R1Vars[R1TensorVarCounts]))[$ivar]; R1TensorVarCounts"
		  . "$incr;\n";
	}
	for ( my $i = 0 ; $i < @$r2Vars ; $i++ ) {
		$str .=
"    $$r2Vars[$i] = (*(R2Vars[R2TensorVarCounts]))[$ivar]; R2TensorVarCounts"
		  . "$incr;\n";
	}
	for ( my $i = 0 ; $i < @$r2sVars ; $i++ ) {
		$str .=
"    $$r2sVars[$i] = (*(R2SymVars[R2SymTensorVarCounts]))[$ivar]; R2SymTensorVarCounts"
		  . "$incr;\n";
	}
	$str .= "  }
";
	return $str;
}

#___________________________________________

sub Operators {
	my $classType = shift;
	my $className = shift;
	my $baseName  = shift;
	my $prodsig   = "";
	my $plussig   = "";
	my $esig      = "";
	if ( $classType =~ m/base/ ) {
	}
	else {
		$prodsig = "
    base::operator*=(factor);";
		$plussig = "
    base::operator+=(datum);";
		$esig = "
    base::operator=(datum);";
	}

	my $intVars = shift;
	my $r0Vars  = shift;
	my $r1Vars  = shift;
	my $r2Vars  = shift;
	my $r2sVars = shift;

	my $str = "  inline $className" . "&
  operator*=(const realT factor)
  {$prodsig
";
	for ( my $i = 0 ; $i < @$intVars ; $i++ ) {
		$str .= "    $$intVars[$i] *= factor;\n";
	}
	for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
		$str .= "    $$r0Vars[$i] *= factor;\n";
	}
	for ( my $i = 0 ; $i < @$r1Vars ; $i++ ) {
		$str .= "    $$r1Vars[$i] *= factor;\n";
	}
	for ( my $i = 0 ; $i < @$r2Vars ; $i++ ) {
		$str .= "    $$r2Vars[$i] *= factor;\n";
	}
	for ( my $i = 0 ; $i < @$r2sVars ; $i++ ) {
		$str .= "    $$r2sVars[$i] *= factor;\n";
	}
	$str .= "    return *this;
  }

  inline $className" . "&
  operator=(const $className" . "& datum)
  {$esig
";
	for ( my $i = 0 ; $i < @$intVars ; $i++ ) {
		$str .= "    $$intVars[$i] = datum.$$intVars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
		$str .= "    $$r0Vars[$i] = datum.$$r0Vars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r1Vars ; $i++ ) {
		$str .= "    $$r1Vars[$i] = datum.$$r1Vars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r2Vars ; $i++ ) {
		$str .= "    $$r2Vars[$i] = datum.$$r2Vars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r2sVars ; $i++ ) {
		$str .= "    $$r2sVars[$i] = datum.$$r2sVars[$i];\n";
	}
	$str .= "    return *this;
  }

  inline $className" . "&
  operator+=(const $className" . "& datum)
  {$plussig
";
	for ( my $i = 0 ; $i < @$intVars ; $i++ ) {
		$str .= "    $$intVars[$i] += datum.$$intVars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
		$str .= "    $$r0Vars[$i] += datum.$$r0Vars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r1Vars ; $i++ ) {
		$str .= "    $$r1Vars[$i] += datum.$$r1Vars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r2Vars ; $i++ ) {
		$str .= "    $$r2Vars[$i] += datum.$$r2Vars[$i];\n";
	}
	for ( my $i = 0 ; $i < @$r2sVars ; $i++ ) {
		$str .= "    $$r2sVars[$i] += datum.$$r2sVars[$i];\n";
	}
	$str .= "    return *this;
  }
";
	return $str;
}

#___________________________________________

sub MapFunctions {
	my $classType = shift;
	my $className = shift;
	my $baseName  = shift;
	my $to_base   = "";
	my $from_base = "";
	if ( $classType =~ m/base/ ) {
	}
	else {
		$to_base = "
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);";
		$from_base = "
    base::MapFromRegion(p0, p1, fct0, fct1);";
	}

	my $intVars = shift;
	my $r0Vars  = shift;
	my $r1Vars  = shift;
	my $r2Vars  = shift;
	my $r2sVars = shift;

	my $tVarsPtr = shift;
	my %tVars    = %$tVarsPtr;

	my $str = "  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, $className" . "& p0, $className" . "& p1)
  {$to_base
";
	for ( my $i = 0 ; $i < @$intVars ; $i++ ) {
		my $var = $$intVars[$i];
		$str .= sprintf(
"    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, %s, p0.%s, p1.%s);\n",
			$var, $var, $var );
	}
	for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
		my $var = $$r0Vars[$i];
		$str .= sprintf(
"    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, %s, p0.%s, p1.%s);\n",
			$var, $var, $var );
	}
	for ( my $i = 0 ; $i < @$r1Vars ; $i++ ) {
		my $var = $$r1Vars[$i];
		if ( exists $tVars{$var} ) {
			$str .= sprintf(
"    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, %s, p0.%s, p1.%s, true);\n",
				$var, $var, $var );
		}
		else {
			$str .= sprintf(
"    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, %s, p0.%s, p1.%s);\n",
				$var, $var, $var );
		}
	}
	for ( my $i = 0 ; $i < @$r2Vars ; $i++ ) {
		my $var = $$r2Vars[$i];
		$str .= sprintf(
"    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, %s, p0.%s, p1.%s);\n",
			$var, $var, $var );
	}
	for ( my $i = 0 ; $i < @$r2sVars ; $i++ ) {
		my $var = $$r2sVars[$i];
		$str .= sprintf(
"    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, %s, p0.%s, p1.%s);\n",
			$var, $var, $var );
	}
	$str .= "\n  }\n\n";

	$str .=
	    "  void MapFromRegion(const $className"
	  . "& p0, const $className"
	  . "& p1, const realT fct0,
                     const realT fct1)
  {$from_base
";
	for ( my $i = 0 ; $i < @$intVars ; $i++ ) {
		my $var = $$intVars[$i];
		$str .= sprintf(
"    GeometryUtilities::MapFromRegion(p0.%s, p1.%s, fct0, fct1, %s);\n",
			$var, $var, $var );
	}
	for ( my $i = 0 ; $i < @$r0Vars ; $i++ ) {
		my $var = $$r0Vars[$i];
		$str .= sprintf(
"    GeometryUtilities::MapFromRegion(p0.%s, p1.%s, fct0, fct1, %s);\n",
			$var, $var, $var );
	}
	for ( my $i = 0 ; $i < @$r1Vars ; $i++ ) {
		my $var = $$r1Vars[$i];
		if ( exists $tVars{$var} ) {
			$str .= sprintf(
"    GeometryUtilities::MapFromRegion(p0.%s, p1.%s, fct0, fct1, %s, true);\n",
				$var, $var, $var );
		}
		else {
			$str .= sprintf(
"    GeometryUtilities::MapFromRegion(p0.%s, p1.%s, fct0, fct1, %s);\n",
				$var, $var, $var );
		}
	}
	for ( my $i = 0 ; $i < @$r2Vars ; $i++ ) {
		my $var = $$r2Vars[$i];
		$str .= sprintf(
"    GeometryUtilities::MapFromRegion(p0.%s, p1.%s, fct0, fct1, %s);\n",
			$var, $var, $var );
	}
	for ( my $i = 0 ; $i < @$r2sVars ; $i++ ) {
		my $var = $$r2sVars[$i];
		$str .= sprintf(
"    GeometryUtilities::MapFromRegion(p0.%s, p1.%s, fct0, fct1, %s);\n",
			$var, $var, $var );
	}
	$str .= "\n  }\n";

	return $str;
}

#___________________________________________

sub ClassDeclaration {
	my $classType = shift;
	my $className = shift;
	my $baseName  = shift;

	my $sdName = $className . "StateData";
	my $pdName = $className . "ParameterData";

	my $filecats = shift;

	my $public_str    = &GetFunctionDeclarations( $$filecats[0] );
	my $protected_str = &GetFunctionDeclarations( $$filecats[1] );
	my $private_str   = &GetFunctionDeclarations( $$filecats[2] );
	if ( length($protected_str) > 0 ) {
		my $tstr = "protected:\n$protected_str";
		$protected_str = $tstr;
	}
	
	my $baseType  = shift;
	my $baseNameFunctionDeclaration = "\n";
	if(length($baseType) > 0)
	{
		$baseNameFunctionDeclaration = 	"\n  inline std::string BaseName() { return \"$baseType\"; }\n";	
	}

	my $sizes = "";
	my $initMatState = "";
	if ( $classType =~m/base/ ) {
		$sizes = "
  const int m_paramSize;
  const int m_stateSize;\n";
  		$initMatState = "\n  virtual void InitializeStates( const localIndex index ){}\n";
	}

	my $str = "

//**********************************************************************************************************************
//**********************************************************************************************************************


class $className: public $baseName
{
public:$sizes
  
  typedef $className" . "ParameterData ParameterClass;
  typedef $className" . "StateData     StateClass;
  $baseNameFunctionDeclaration
  $className" . "( const int paramSize, const int stateSize );

  virtual ~$className" . "();
  
  virtual void ReadXML( const TICPP::HierarchicalDataNode& node ) = 0;

  virtual void resize( const localIndex num ) = 0;
  
  virtual void resize( const localIndex num0,
                       const localIndex num1 ) = 0;
  
  virtual void insert( const localIndex num ) = 0;

  virtual void erase( const localIndex num ) = 0;
  $initMatState
  virtual const $sdName* StateData( const localIndex index0,
                                                  const localIndex index1 ) const = 0;
  virtual       $sdName* StateData( const localIndex index0,
                                                  const localIndex index1 )  = 0;

  virtual const $pdName* ParameterData( const localIndex index ) const = 0;
  virtual       $pdName* ParameterData( const localIndex index ) = 0;
  
  inline void IncrementPtr( const $sdName" . "* ptr ) const
  {
    ptr = reinterpret_cast<const $sdName"
	  . "*>( reinterpret_cast<const char*>(ptr) + m_stateSize );
  }

  inline void IncrementPtr( const $pdName" . "* ptr ) const
  {
    ptr = reinterpret_cast<const $pdName"
	  . "*>( reinterpret_cast<const char*>(ptr) + m_paramSize );
  }
  
  inline void MapToRegion(const realT fctNormal, const realT fct0, const realT fct1, 
                          const localIndex from0, const localIndex from1,
                          StateClass& s0, StateClass& s1)
  {
    StateData(from0, from1)->MapToRegion(fctNormal, fct0, fct1, s0, s1);
    //ParameterData(from)->MapToRegion(fctNormal, fct0, fct1, p0, p1);
  }

  inline void MapFromRegion(const realT fct0, const realT fct1, 
                          const StateClass& s0, const StateClass& s1, 
                          const localIndex to0, const localIndex to1)
  {
    StateData(to0, to1)->MapFromRegion(s0, s1, fct0, fct1);
    //ParameterData(to)->MapFromRegion(p0, p1, fct0, fct1);
  }
  
  inline void MapToRegion(const realT fctNormal, const realT fct0, const realT fct1, 
                          const localIndex from,
                          ParameterClass& p0, ParameterClass& p1)
  {
    //StateData(from0, from1)->MapToRegion(fctNormal, fct0, fct1, s0, s1);
    ParameterData(from)->MapToRegion(fctNormal, fct0, fct1, p0, p1);
  }

  inline void MapFromRegion(const realT fct0, const realT fct1, 
                          const ParameterClass& p0, const ParameterClass& p1, 
                          const localIndex to)
  {
    //StateData(to0, to1)->MapFromRegion(s0, s1, fct0, fct1);
    ParameterData(to)->MapFromRegion(p0, p1, fct0, fct1);
  }
  
  virtual void ZeroStates() = 0;
    
  virtual localIndex NumStateIndex0() const = 0;
  virtual localIndex NumStateIndex1() const = 0;

  virtual localIndex NumParameterIndex0() const = 0;
  virtual localIndex NumParameterIndex1() const = 0;

$public_str
$protected_str  
private:
  $className" . "();
  $className" . "( const $className" . "& );
  $className" . "& operator=( const $className" . "& );
  
$private_str  
};
";
	return $str;
}

sub ClassDeclarationLeaf {
	my $className = shift;
	my $baseName  = shift;

	my $filecats = shift;

	my $public_str    = &GetFunctionDeclarations( $$filecats[0] );
	my $protected_str = &GetFunctionDeclarations( $$filecats[1] );
	my $private_str   = &GetFunctionDeclarations( $$filecats[2] );
	if ( length($protected_str) > 0 ) {
		my $tstr = "protected:\n$protected_str";
		$protected_str = $tstr;
	}

	my $str = "
//**********************************************************************************************************************
//**********************************************************************************************************************
class $className" . ": public $baseName
{
public:

  typedef $className" . "ParameterData ParameterClass;
  typedef $className" . "StateData     StateClass;

  typedef Array1dT<ParameterClass> ParameterArrayType;
  typedef Array2dT<StateClass>     StateArrayType;


  StateArrayType m_stateData;
  ParameterArrayType m_parameterData;

  localIndex NumStateIndex0() const { return m_stateData.Dimension(0); }
  localIndex NumStateIndex1() const { return m_stateData.Dimension(1); }

  localIndex NumParameterIndex0() const { return m_parameterData.size(); }
  localIndex NumParameterIndex1() const { return 1; }

  static std::string Name() { return \"$className\"; }

  $className" . "();
  virtual ~$className" . "();

  inline void MapToRegion(const realT fctNormal, const realT fct0, const realT fct1, 
                          const localIndex from0, const localIndex from1,
                          StateClass& s0, StateClass& s1)
  {
    StateData(from0, from1)->MapToRegion(fctNormal, fct0, fct1, s0, s1);
    //ParameterData(from)->MapToRegion(fctNormal, fct0, fct1, p0, p1);
  }

  inline void MapFromRegion(const realT fct0, const realT fct1, 
                          const StateClass& s0, const StateClass& s1, 
                          const localIndex to0, const localIndex to1)
  {
    StateData(to0, to1)->MapFromRegion(s0, s1, fct0, fct1);
    //ParameterData(to)->MapFromRegion(p0, p1, fct0, fct1);
  }
 
  inline void MapToRegion(const realT fctNormal, const realT fct0, const realT fct1, 
                          const localIndex from,
                          ParameterClass& p0, ParameterClass& p1)
  {
    //StateData(from0, from1)->MapToRegion(fctNormal, fct0, fct1, s0, s1);
    ParameterData(from)->MapToRegion(fctNormal, fct0, fct1, p0, p1);
  }

  inline void MapFromRegion(const realT fct0, const realT fct1, 
                          const ParameterClass& p0, const ParameterClass& p1, 
                          const localIndex to)
  {
    //StateData(to0, to1)->MapFromRegion(s0, s1, fct0, fct1);
    ParameterData(to)->MapFromRegion(p0, p1, fct0, fct1);
  }
  
  
  virtual void ZeroStates()
  {
    for(localIndex j = 0; j < m_stateData.Dimension(1); j++)
    {
      for(localIndex i = 0; i < m_stateData.Dimension(0); i++)
      {
        m_stateData(i,j) *= 0.0;
      }
    }
  }
  
  virtual void SetVariableParameters(const bool varParams, const localIndex newSize = 0)
  { SetVariableParametersFromDerived<$className>(varParams, newSize); }

  virtual void ReadXML( const TICPP::HierarchicalDataNode& node )
  { ReadXMLFromDerived<$className>( node ); }

  virtual void resize( const localIndex num )
  { ResizeFromDerived<$className>( num ); }

  virtual void resize( const localIndex num0, const localIndex num1 )
  {
    m_stateData.resize2(num0, num1);
    ResizeFromDerived<$className>( num0 );
  }
  
  virtual void insert( const localIndex num )
  { InsertFromDerived<$className>( num ); }

  virtual void erase( const localIndex num )
  { EraseFromDerived<$className>( num ); }
 
  void GetVariableNames( sArray1d& intVars, sArray1d& realVars, sArray1d& R1TensorVars, sArray1d& R2TensorVars, sArray1d& R2SymTensorVars ) const
  { GetVariableNamesFromDerived<$className>(intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars ); }

  size_t GetStateOffset( const std::string& name, const int type ) const
  { return GetStateOffsetFromDerived<$className>(name, type); }
  
  size_t GetParameterOffset( const std::string& name, const int type ) const
  { return GetParameterOffsetFromDerived<$className>(name, type ); }
  
  bool GetStateValues( const std::string& name, rArray1d& values ) const
  { return GetStateValuesFromDerived<$className>(name, values); }

  bool GetParameterValues( const std::string& name, rArray1d& values ) const
  { return GetParameterValuesFromDerived<$className>(name, values); }

  bool SetStateValues( const std::string& name, const rArray1d& values )
  { return SetStateValuesFromDerived<$className>(name, values); }

  bool SetParameterValues( const std::string& name, const rArray1d& values )
  { return SetParameterValuesFromDerived<$className>(name, values); }

  virtual void Serialize( Array1dT<iArray1d*>& intVars, Array1dT<rArray1d*>& realVars, Array1dT<Array1dT<R1Tensor>*>& R1Vars, Array1dT<Array1dT<R2Tensor>*>& R2Vars, Array1dT<Array1dT<R2SymTensor>*>& R2SymVars ) const
  { SerializeFromDerived<$className>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual void Deserialize( const Array1dT<iArray1d*>& intVars, const Array1dT<rArray1d*>& realVars, const Array1dT<Array1dT<R1Tensor>*>& R1Vars, const Array1dT<Array1dT<R2Tensor>*>& R2Vars, const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars  )
  { DeserializeFromDerived<$className>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<$className>( localIndices, buffer, doBufferPacking ); }
  unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<$className>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<$className>( localIndices, buffer, doBufferPacking ); }
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<$className>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer )
  { return UnpackFromDerived<$className>( localIndices, buffer ); }

  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer )
  { return UnpackFromDerived<$className>( localIndices, buffer ); }

  const StateClass* StateData( const localIndex index0, const localIndex index1 ) const
  { return &( m_stateData(index0,index1) );  }
  StateClass* StateData( const localIndex index0, const localIndex index1 )
  { return &( m_stateData(index0,index1) );  }

  const ParameterClass* ParameterData( const localIndex index ) const
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  ParameterClass* ParameterData( const localIndex index )
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  
$public_str
$protected_str
private:
  $className" . "(const $className" . "&);
  $className" . "& operator=(const $className" . "&);
  
$private_str
};
";
	return $str;
}

###################################
# PRINT CLASS FILE
###################################

sub PrintClassFile {
	my $classType = shift;
	my $className = shift;
	my $baseName  = shift;
	my $baseType = shift;
	my $ffiles = shift;

	my $heading = &Heading( $className . ".cpp" );

	open( OUT, ">$className" . ".cpp" );

	my $disclaimer = &Disclaimer();
	
	my $factory = ucfirst($baseType);
	$factory = "Constitutive/$factory/$factory"."Factory.h";
	my $registration = uc($baseType);
	$registration = "
/// Register class in the class factory
REGISTER_$registration( $className )
";

	print OUT "$disclaimer

$heading

#include \"Utilities/GeometryUtilities.h\"
#include \"Utilities/FindRoots.h\"
#include \"Utilities/MaterialUtilities.h\"
#include \"IO/ticpp/HierarchicalDataNode.h\"

#include \"$className" . ".h\"
";
	if ( $classType =~ m/base/ ) {

		print OUT "
$className" . "::$className" . "( const int paramSize, const int stateSize ):
$baseName" . "(),
m_paramSize( paramSize ),
m_stateSize( stateSize )
{

}
";
	}
	elsif ( $classType =~ m/intermediate/ ) {
		print OUT "#include <typeinfo>
#include <assert.h>

$className" . "::$className" . "( const int paramSize, const int stateSize ):
$baseName" . "( paramSize, stateSize )
{
  // TODO Auto-generated constructor stub

}
";
	}
	else {
		print OUT "#include \"IO/ticpp/HierarchicalDataNode.h\"
#include \"$factory\"
#include <typeinfo>
#include <assert.h>

$className" . "::$className" . "( ):
$baseName" . "( sizeof(ParameterClass), sizeof(StateClass) )
{
  // TODO Auto-generated constructor stub
}
";
	}

	print OUT "
$className" . "::~$className" . "()
{

}

";

	#Print the files directly into the class file
	for ( my $i = 0 ; $i <= $#$ffiles ; $i++ ) {
		my $files = $$ffiles[$i];
		for ( my $j = 0 ; $j <= $#$files ; $j++ ) {
			my $str = &GetFunctionImplementations( $$files[$j] );
			print OUT $str;
		}
	}
	
	if ( $classType =~ m/leaf/ )
	{
		print OUT $registration;
	}
	close OUT;
}

sub GetFunctionImplementations {
	my $filenames = shift;
	my $str       = "";
	for ( my $i = 0 ; $i <= $#$filenames ; $i++ ) {
		my $filename = $$filenames[$i];

		#die "GetFunctionImplementations $filename\n";
		die "cannot open $filename implementation include file\n"
		  if ( !open( INC, "<$filename" ) );
		while (<INC>) {
			if (/FUNCTION_BEGIN_PARSE/) {
			}
			else {
				my $line = $_;
				$line =~ s/virtual_//;
				$line =~ s/static_realT/realT/;
				$line =~ s/static_void/void/;
				$line =~ s/static_bool/bool/;
				$str .= $line;
			}
		}
		close INC;
	}
	return $str;
}

sub GetFunctionDeclarations {
	my $filenames = shift;
	my $str       = "";
	for ( my $i = 0 ; $i <= $#$filenames ; $i++ ) {
		my $filename = $$filenames[$i];

		#die "GetFunctionDeclarations $filename\n";
		die "cannot open $filename declaration include file\n"
		  if ( !open( INC, "<$filename" ) );
		my $read_on = 0;
		my $tstr    = "";
		while (<INC>) {
			my $line = $_;
			if (/FUNCTION_BEGIN_PARSE/) {
				$read_on = 1;
			}
			elsif ( $read_on == 1 ) {
				if ( $line =~ m/\)/ ) {

					#find the ending parenthesis
					my @arr = split( /\)/, $line );
					$tstr .= "  $arr[0]";
					if ( $arr[1] =~ m/\s*const/ ) {
						$tstr .= ") const;\n\n";
					}
					else {
						$tstr .= ");\n\n";
					}
					$read_on = 0;

					#remove the class scoping
					@arr = split( /\:\:/, $tstr );
					die "I'm confused ... we found no :: ... $tstr\n"
					  if ( $#arr < 1 );
					$tstr = $arr[0];
					$tstr =~ s/\S+$//;
					$tstr .= $arr[1];				
					for(my $j = 2; $j <=$#arr; $j++)
					{
						$tstr .= "::$arr[$j]";
					}
					#add the temporary string to the return value and reset
					$str .= $tstr;
					$tstr = "";
				}
				else {
					$line =~ s/virtual_/virtual /;
					$line =~ s/static_realT/static realT/;
					$line =~ s/static_void/static void/;
					$line =~ s/static_bool/static bool/;
					$tstr .= "  $line";
				}
			}
		}
		close INC;
	}
	return $str;
}

sub Disclaimer {
	my $str =
"//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2013, Lawrence Livermore National Security, LLC.
//
//  Produced at the Lawrence Livermore National Laboratory
//
//  Written by:
//  Randolph Settgast (settgast1\@llnl.gov)     Stuart Walsh(walsh24\@llnl.gov)
//  Scott Johnson (johnson346\@llnl.gov)        Pengcheng Fu (fu4\@llnl.gov)
//  Joshua White (white230\@llnl.gov)           
//
//  LLNL-CODE-618232
//  GPAC, Version 2.0
//
//  All rights reserved.
//
//  This file is part of GPAC. For details, please contact Scott Johnson or Randolph Settgast. Please also read \"Additional BSD Notice\" below.
//
//  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.
//  * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer (as noted below) in the 
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Additional BSD Notice
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
";
	return $str;
}
