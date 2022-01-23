#!/bin/perl
use strict;
use Math::Trig;

my $narg = $#ARGV + 1;
die "usage: $0 <point file where z along the bore hole axis> <poissons ratio> <sigma_v> <sigma_H> <sigma_h> <wellbore pressure> <borehole radius> <angle between well bore and sigma_H direction - degrees>\n" if($narg != 8);

my $fname = shift;
die "Cannot open $fname\n" if(!open(IN,"<$fname"));

my $nu = shift;
my $sigv = shift;
my $sigH = shift;
my $sigh = shift;
my $pp = shift;
my $radius = shift;
my $alpha = shift;
$alpha *= pi / 180;

my @params = ($nu,$sigv,$sigH,$sigh,$pp,$radius,$alpha); 

#print "x y z vx vy vz vv\n";
while(<IN>)
{
    my $line = $_;
    $line=~s/^\s+//;
    $line=~s/\s+$//;
    my @arr = split(/\s+/,$line);
    die "Cannot read line: \"$line\" - it does not have at least 3 entries\n" if($#arr < 2);
    my ($r,$theta,$z) = &CartesianToCylindrical(@arr);

    my @sigcyl = &CylindricalStress($r,$theta,$z,@params);
    #my @sigcyl = ($rr,$tt,$zz,$rt,$tz,$rz);
    if(0)
    {
	my $str = "SIGCYL: $sigcyl[0]";
	for(my $i = 1; $i <= $#sigcyl; $i++)
	{
	    $str .= " $sigcyl[$i]";
	}
	$str .= "\n";
	die $str;
    }

    my @sigcart = &CylindricalStressToCartesian($theta, @sigcyl);
    if(0)
    {
	my $str = "SIGCART: $sigcart[0]";
	for(my $i = 1; $i <= $#sigcart; $i++)
	{
	    $str .= " $sigcart[$i]";
	}
	$str .= "\n";
	die $str;
    }
    #print "$arr[0] $arr[1] $arr[2] $sigcart[0] $sigcart[1] $sigcart[2] $sigcart[3] $sigcart[4] $sigcart[5]\n";

    my @lambda = &EigenVals(@sigcart);
    my @eig = &Eigen(@lambda,@sigcart);

    #print "$arr[0] $arr[1] $arr[2] $eig[0] $eig[1] $eig[2] $lambda[0]\n";
    {
	my @xx = ($arr[0], $arr[1], $arr[2]);
	my $dot = 0;
	for(my $i = 0; $i <= $#xx; $i++)
	{
	    $dot += $xx[$i] * $eig[$i];
	}
	my $fct = 1;
	if($dot < 0)
	{
	    $fct = -1;
	}
	for(my $i = 0; $i <= $#xx; $i++)
	{
	    $eig[$i] *= $fct;
	}
    }

    print "$arr[0] $arr[1] $arr[2] $sigcart[0] $sigcart[1] $sigcart[2] $sigcart[3] $sigcart[4] $sigcart[5] $eig[0] $eig[1] $eig[2] $lambda[0]\n";
}
close IN;

sub CartesianToCylindrical
{
    my $x = shift;
    my $y = shift;
    my $z = shift;

    my $r = sqrt($x * $x + $y * $y);
    die "x y z: $x $y $z -> r = 0\n" if($r == 0);
    my $theta = atan2($y,$x);
    return ($r,$theta,$z);
}

sub CylindricalStress
{
    my $r = shift;
    my $theta = shift;
    my $z = shift;
    my $nu = shift;
    my $sigv = shift;
    my $sigH = shift;
    my $sigh = shift;
    my $pp = shift;
    my $radius = shift;
    my $alpha = shift;

    #die "$r $theta $z $nu $sigv $sigH $sigh $pp $radius $alpha\n";

    my $R2divr2 = $radius * $radius / ($r * $r);
    my $R4divr4 = $R2divr2 * $R2divr2;
    my $pR2divr2 = $pp * $R2divr2;

    my $s2a = sin($alpha);
    $s2a *= $s2a;
    my $c2a = cos($alpha);
    $c2a *= $c2a;

    my $sigHs2a = $sigH * $s2a;
    my $sighc2a = $sigh * $c2a;

    my $c2t = cos(2 * $theta);
    my $s2t = sin(2 * $theta);

    my $srr = 0.5 * ($sigv + $sigHs2a + $sighc2a) * (1 - $R2divr2);
    $srr += 0.5 * ($sigv - $sigHs2a - $sighc2a) * (1 - 4 * $R2divr2 + 3 * $R4divr4) * $c2t;
    $srr += $pR2divr2;

    my $stt = 0.5 * ($sigv + $sigHs2a + $sighc2a) * (1 + $R2divr2);
    $stt -= 0.5 * ($sigv - $sigHs2a - $sighc2a) * (1 + 3 * $R4divr4) * $c2t;
    $stt -= $pR2divr2;

    my $szz = $sigH * $c2a;
    $szz += $sigh * $s2a;
    $szz -= 2 * $nu * ($sigv - $sigHs2a - $sighc2a) * $R2divr2 * $c2t;

    my $srt = 0.5 * ($sigHs2a + $sighc2a - $sigv);
    $srt *= 1 + 2 * $R2divr2 - 3*$R4divr4;
    $srt *= $s2t;

    my $stz = 0.5 * ($sigh - $sigH);
    $stz *= 1 + $R2divr2;
    $stz *= sin(2*$alpha) * cos($theta);

    my $srz = 0.5 * ($sigh - $sigH);
    $srz *= 1 - $R2divr2;
    $srz *= sin(2*$alpha) * sin($theta);

    return ($srr,$stt,$szz,$srt,$srz,$stz);
}

sub CylindricalStressToCartesian
{
    my $theta = shift;
    $theta *= -1;

    my $srr = shift;
    my $stt = shift;
    my $szz = shift;
    my $srt = shift;
    my $srz = shift;
    my $stz = shift;

    #transform through angle theta
    my $ct = cos($theta);
    my $st = sin($theta);
    my $ct2  = $ct * $ct;
    my $st2  = $st * $st;
    my $ctst = $ct * $st;

    my $sxx = $stt * $ct2 + $srr * $st2 + 2 * $srt * $ctst;
    my $syy = $stt * $st2 + $srr * $ct2 - 2 * $srt * $ctst;
    my $sxy = ($srr - $stt) * $ctst + $srt * ($ct2 - $st2);

    #transform orthogonal directions
    my $sxz = $ct * $stz + $st * $srz;
    my $syz = $ct * $srz - $st * $stz; 
    return ($sxx, $syy, $szz, $sxy, $sxz, $syz); 
}

sub EigenVals
{
    my $A11 = shift;
    my $scale = $A11;
    my $A22 = shift;
    if($A22 > $scale)
    {
	$scale = $A22;
    }
    my $A33 = shift;
    if($A33 > $scale)
    {
	$scale = $A33;
    }
    my $A12 = shift;
    if($A12 > $scale)
    {
	$scale = $A12;
    }
    my $A13 = shift;
    if($A13 > $scale)
    {
	$scale = $A13;
    }
    my $A23 = shift;
    if($A23 > $scale)
    {
	$scale = $A23;
    }

    my $Inv3 = 1 / 3.0;
    my $sigma = $Inv3 * ($A11 + $A22 + $A33);
    my @eigenvals = (0,0,0);
    if($scale == 0.0)
    {
	return @eigenvals;
    }

    my $iscale = 1 / $scale;
    $A11 -= $sigma;
    $A22 -= $sigma;
    $A33 -= $sigma;

    $A11 *= $iscale;
    $A22 *= $iscale;
    $A33 *= $iscale;
    $A12 *= $iscale;
    $A13 *= $iscale;
    $A23 *= $iscale;

    my $A12A12 = $A12 * $A12;
    my $A13A13 = $A13 * $A13;
    my $A23A23 = $A23 * $A23;
    my $A12A13A23 = $A12 * $A13 * $A23;

    my $II = 0.5 * ($A11 * $A11 + $A22 * $A22 + $A33 * $A33) + ($A12A12 + $A13A13 + $A23A23);
    my $III = $A11 * ($A22 * $A33 - $A23A23) - ($A12A12 * $A33 - $A12A13A23) + ($A12A13A23 - $A13A13 * $A22);

    my @theta = (0,0,0);

    if ($II > 0.0)
    {
	my $temp = 3.0 / $II;
	my $sqrt_3divII = sqrt( $temp );
	my $arg = 0.5 * $III * $temp * $sqrt_3divII;

	if ($arg > 1.0)
	{
	    $arg = 1.0;
	}
	if ($arg < -1.0)
	{
	    $arg = -1.0;
	}

	$temp = 2 * pi * $Inv3;
	$theta[0] = acos( $arg ) * $Inv3;
	$theta[1] = $theta[0] + $temp;
	$theta[2] = $theta[0] - $temp;

	$temp = 2.0 / $sqrt_3divII * $scale;
	for (my $i = 0; $i < 3; $i++)
	{
	    $eigenvals[$i] = $temp * cos( $theta[$i] ) + $sigma;
	}
    }
    else
    {
	for (my $i = 0; $i < 3; $i++)
	{
	    $eigenvals[$i] = $sigma;
	}
    }
    for (my $i = 0; $i < 3; $i++)
    {
	for (my $j = 0; $j < 3; $j++)
	{
	    if ($eigenvals[$j] < $eigenvals[$j + 1])
	    {
		my $temp = $eigenvals[$j];
		$eigenvals[$j] = $eigenvals[$j + 1];
		$eigenvals[$j + 1] = $temp;
	    }
	}
    }
    return @eigenvals;
}

sub EigenVecs
{
    my @lambda = (0,0,0);
    for(my $i = 0; $i < 3; $i++)
    {
	$lambda[$i] = shift; 
    }
    my @M = (0,0,0,0,0,0);
    for(my $i = 0; $i < 6; $i++)
    {
	$M[$i] = shift; 
    }
    
    my @v1 = (0,0,0);
    my @v2 = (0,0,0);
    my @v3 = (0,0,0);

    my $tol = 1e-14;
    #compare eigenvalues for multiplicity                                                                                                                
    #due to the fact that eigenvalues must be ordered                                                                                                    
    #either e1=e2, e2=e3 or e1=e2=e3                                                                                                                     
    my $multiplicity_flag = 0;
    #if e1=e2 then multiplicity flag = 1;        
    my $tmp = $lambda[0]-$lambda[1];
    if($tmp < 0)
    {
	$tmp = -$tmp;
    }
    if ($tmp < $tol)
    {
	$multiplicity_flag = 1;
    }
    #if e2=e3 then multiplicity flag += 2;
    $tmp = $lambda[1] - $lambda[2];
    if($tmp < 0)
    {
	$tmp = -$tmp;
    }
    if ($tmp < $tol)
    {
	$multiplicity_flag += 2;
    }
    #if all eigenvalues are the same, then eigenvectors are any 3 orthoganal                                                                             
    #vectors.                                                                                                                                            
    if ($multiplicity_flag == 3)
    {
	$v1[0] = 1;
	$v2[1] = 1;
	$v3[2] = 1;
    }
    else
    {
	my @vt1 = (0,0,0);
	my @vt2 = (0,0,0);
	my @vt3 = (0,0,0);
	
	my @v_mag = (0,0,0);
	my $maxmag = 0;
	
	$M[0] -= $lambda[1];
	$M[1] -= $lambda[1];
	$M[2] -= $lambda[1];
	
	$maxmag = 0.0;
	if ($multiplicity_flag > 0)
	{
	    my $imax = -1;
	    for (my $i = 0; $i < 6; $i++)
	    {
		$tmp = $M[$i];
		if($tmp < 0)
		{
		    $tmp = -$tmp;
		}
		if ($maxmag < $tmp)
		{
		    $maxmag = $tmp;
		    $imax = $i;
		}
	    }
	    if ($imax == 0 || $imax == 1)
	    {
		$vt1[0] = -$M[1];
		$vt1[1] = $M[0];
		$vt1[2] = 0;

		$vt2[0] = -$M[3] * $M[0];
		$vt2[1] = -$M[3] * $M[1];
		$vt2[2] = $M[0] * $M[0] + $M[1] * $M[1];
		
		$tmp = 0;
		my $tmp2 = 0;
		for(my $i = 0; $i < 3; $i++)
		{    
		    $tmp += $vt1[$i] * $vt1[$i];
		    $tmp2 += $vt2[$i] * $vt2[$i];
		}
		$tmp = 1.0 / sqrt($tmp);
		$tmp2 = 1.0 / sqrt($tmp2);
		for(my $i = 0; $i < 3; $i++)
		{    
		    $vt1[$i] *= $tmp;
		    $vt2[$i] *= $tmp;
		}
	    }
	    elsif ($imax == 3)
	    {
		$vt1[0] = $M[3];
		$vt1[1] = 0;
		$vt1[2] = -$M[0];
		
		$vt2[0] = -$M[1] * $M[0];
		$vt2[1] = $M[0] * $M[0] + $M[3] * $M[3];
		$vt2[2] = -$M[1] * $M[3];

		$tmp = 0;
		my $tmp2 = 0;
		for(my $i = 0; $i < 3; $i++)
		{    
		    $tmp += $vt1[$i] * $vt1[$i];
		    $tmp2 += $vt2[$i] * $vt2[$i];
		}
		$tmp = 1.0 / sqrt($tmp);
		$tmp2 = 1.0 / sqrt($tmp2);
		for(my $i = 0; $i < 3; $i++)
		{    
		    $vt1[$i] *= $tmp;
		    $vt2[$i] *= $tmp;
		}		
	    }
	    elsif ($imax == 2 || $imax == 4)
	    {
		$vt1[0] = 0;
		$vt1[1] = $M[4];
		$vt1[2] = $M[2];
		
		$vt2[0] = $M[2] * $M[2] + $M[4] * $M[4];
		$vt2[1] = -$M[1] * $M[2];
		$vt2[2] = -$M[1] * $M[4];
		
		$tmp = 0;
		my $tmp2 = 0;
		for(my $i = 0; $i < 3; $i++)
		{    
		    $tmp += $vt1[$i] * $vt1[$i];
		    $tmp2 += $vt2[$i] * $vt2[$i];
		}
		$tmp = 1.0 / sqrt($tmp);
		$tmp2 = 1.0 / sqrt($tmp2);
		for(my $i = 0; $i < 3; $i++)
		{    
		    $vt1[$i] *= $tmp;
		    $vt2[$i] *= $tmp;
		}		
	    }
	    elsif ($imax == 5)
	    {
		$vt1[0] = 0;
		$vt1[1] = -$M[5];
		$vt1[2] = $M[4];
		
		$vt2[0] = $M[4] * $M[4] + $M[5] * $M[5];
		$vt2[1] = -$M[3] * $M[4];
		$vt2[2] = -$M[3] * $M[5];
		
		$tmp = 0;
		my $tmp2 = 0;
		for(my $i = 0; $i < 3; $i++)
		{    
		    $tmp += $vt1[$i] * $vt1[$i];
		    $tmp2 += $vt2[$i] * $vt2[$i];
		}
		$tmp = 1.0 / sqrt($tmp);
		$tmp2 = 1.0 / sqrt($tmp2);
		for(my $i = 0; $i < 3; $i++)
		{    
		    $vt1[$i] *= $tmp;
		    $vt2[$i] *= $tmp;
		}		
	    }
	    
	    if ($multiplicity_flag == 1)
	    {
		for(my $i = 0; $i < 3; $i++)
		{
		    $v1[$i] = $vt1[$i];
		    $v2[$i] = $vt2[$i];
		}
		$v2[0] = $v1[1] * $v2[2] - $v1[2] * $v2[1];
		$v2[1] = $v1[2] * $v2[0] - $v1[0] * $v2[2];
		$v2[2] = $v1[0] * $v2[1] - $v1[1] * $v2[0];
	    }
	    elsif ($multiplicity_flag == 2)
	    {
		for(my $i = 0; $i < 3; $i++)
		{
		    $v2[$i] = $vt1[$i];
		    $v3[$i] = $vt2[$i];
		}
		$v1[0] = $v2[1] * $v3[2] - $v2[2] * $v3[1];
		$v1[1] = $v2[2] * $v3[0] - $v2[0] * $v3[2];
		$v1[2] = $v2[0] * $v3[1] - $v2[1] * $v3[0];                                                                                                    
	    }
	}
	else
	{
	    for (my $i = 0; $i < 3; $i++)
	    {
		$maxmag = 0.0;
		
		$vt1[0] = $M[2] * $M[5] - $M[4] * $M[4];
		$vt1[1] = $M[3] * $M[4] - $M[1] * $M[5];
		$vt1[2] = $M[1] * $M[4] - $M[3] * $M[2];
		
		$vt2[0] = $vt1[1];
		$vt2[1] = $M[0] * $M[5] - $M[3] * $M[3];
		$vt2[2] = $M[3] * $M[1] - $M[4] * $M[0];
		
		$vt3[0] = $vt1[2];
		$vt3[1] = $vt2[2];
		$vt3[2] = $M[0] * $M[2] - $M[1] * $M[1];


		my $j = 0;

		#1
		$tmp = 0;
		for(my $k = 0; $k < 3; $k++)
		{
		    $tmp += $vt1[$k] * $vt1[$k]; 
		}
		$tmp = sqrt($tmp);
		$v_mag[$j] = $tmp;

		if ($maxmag < $v_mag[$j])
		{
		    $maxmag = $v_mag[$j];
		    $v1[0] = $vt1[0] / $v_mag[$j];
		    $v1[1] = $vt1[1] / $v_mag[$j];
		    $v1[2] = $vt1[2] / $v_mag[$j];
		}

		#2
		$j+=1;
		$tmp = 0;
		for(my $k = 0; $k < 3; $k++)
		{
		    $tmp += $vt2[$k] * $vt2[$k]; 
		}
		$tmp = sqrt($tmp);
		$v_mag[$j] = $tmp;

		if ($maxmag < $v_mag[$j])
		{
		    $maxmag = $v_mag[$j];
		    $v2[0] = $vt2[0] / $v_mag[$j];
		    $v2[1] = $vt2[1] / $v_mag[$j];
		    $v2[2] = $vt2[2] / $v_mag[$j];
		}

		#3
		$j+=1;
		$tmp = 0;
		for(my $k = 0; $k < 3; $k++)
		{
		    $tmp += $vt3[$k] * $vt3[$k]; 
		}
		$tmp = sqrt($tmp);
		$v_mag[$j] = $tmp;

		if ($maxmag < $v_mag[$j])
		{
		    $maxmag = $v_mag[$j];
		    $v3[0] = $vt3[0] / $v_mag[$j];
		    $v3[1] = $vt3[1] / $v_mag[$j];
		    $v3[2] = $vt3[2] / $v_mag[$j];
		}
	    }
	    
	    $tmp = $v1[0] * $v2[0] + $v1[1] * $v2[1] + $v1[2] * $v2[2];   
	    if($tmp < 0)
	    {
		$tmp = -$tmp;
	    }
 	    if ( $tmp > $tol)
	    {
		$v2[0] = $v1[1] * $v3[2] - $v1[2] * $v3[1];
		$v2[1] = $v1[2] * $v3[0] - $v1[0] * $v3[2];
		$v2[2] = $v1[0] * $v3[1] - $v1[1] * $v3[0];
	    }	    
	    else
	    {
		$tmp = $v1[0] * $v3[0] + $v1[1] * $v3[1] + $v1[2] * $v3[2];
		if($tmp < 0)
		{
		    $tmp = -$tmp;
		}
		if ($tmp > $tol)
		{
		    $v1[0] = $v2[1] * $v3[2] - $v2[2] * $v3[1];
		    $v1[1] = $v2[2] * $v3[0] - $v2[0] * $v3[2];
		    $v1[2] = $v2[0] * $v3[1] - $v2[1] * $v3[0];
		}
		else
		{
		    $tmp = $v3[0] * $v2[0] + $v3[1] * $v2[1] + $v3[2] * $v2[2];   
		    if($tmp < 0)
		    {
			$tmp = -$tmp;
		    }
		    if ($tmp > $tol)
		    {
			$v3[0] = $v1[1] * $v2[2] - $v1[2] * $v2[1];
			$v3[1] = $v1[2] * $v2[0] - $v1[0] * $v2[2];
			$v3[2] = $v1[0] * $v2[1] - $v1[1] * $v2[0];
		    }
		}
	    }
	}
    }
    return (@v1, @v2, @v3);
}



sub Eigen
{
    my $sxx = shift;
    my $syy = shift;
    my $szz = shift;
    my $sxy = shift;
    my $sxz = shift;
    my $syz = shift;

    my $p1 = $sxy * $sxy + $sxz * $sxz + $syz * $syz;
    my @eig = (0,0,0,0);
    #die "P1: $p1\n";
    if ($p1 < 1e-5)
    { 
	if($sxx < $syy && $sxx < $szz)
	{
	    $eig[0] = $sxx;
	    $eig[1] = 1;
	}
	elsif($syy < $sxx && $syy < $szz)
	{
	    $eig[0] = $syy;
	    $eig[2] = 1;
	}
	else
	{
	    $eig[0] = $szz;
	    $eig[3] = 1;
	}
    }
    else
    {
	my $q = ($sxx + $syy * $szz)/3.0;
	my $p2 = ($sxx - $q)**2 + ($syy - $q)**2 + ($szz - $q)**2 + 2*$p1;
	#(A(1,1) - q)^2 + (A(2,2) - q)^2 + (A(3,3) - q)^2 + 2 * p1
	my $p = sqrt($p2 / 6.0);

	my $invp = 1.0 / $p;
	my $bxx = $invp * ($sxx - $q);
	my $byy = $invp * ($syy - $q);
	my $bzz = $invp * ($szz - $q);
	my $bxy = $invp * $sxy;
	my $bxz = $invp * $sxz;
	my $byz = $invp * $syz;

	#die "B: 1/p = $invp -- $bxx $byy $bzz $bxy $bxz $byz\n";

	my $detB = $bxx * $byy * $bzz;
	$detB += $bxy * $byz * $bxz;
	$detB += $bxz * $bxy * $byz;
	$detB -= $bxz * $byy * $bxz;
	$detB -= $bxy * $bxy * $bzz;
	$detB -= $bxx * $byz * $byz;

	#die "DETB: $detB\n";

	my $r = $detB * 0.5;
 
	# In exact arithmetic for a symmetric matrix  -1 <= r <= 1
	# but computation error can leave it slightly outside this range.
	my $phi = 0;
	if ($r <= -1) 
	{
	    $phi = pi / 3.0;
	}
	elsif ($r >= 1)
	{
	    $phi = 0;
	}
	else
	{
	    $phi = acos($r) / 3.0;
	}

	#the eigenvalues satisfy eig3 <= eig2 <= eig1
	$eig[0] = $q + 2 * $p * cos($phi);
	#$eig[2] = $q + 2 * $p * cos($phi + (2.0*pi/3));
	#$eig[1] = 3 * $q - $eig[0] - $eig[2];
	#since trace(A) = eig1 + eig2 + eig3

	#get the eigenvectors
	#B00*v0 + B01*v1 + B02*v2 = lam * v0
	#B01*v0 + B11*v1 + B12*v2 = lam * v1
	#B02*v0 + B12*v1 + B22*v2 = lam * v2

	#symmetric order: xx, xy, xz, yy, yz, zz
	my $alpha = $sxy / ($eig[0] - $sxx);
	my $beta = $sxz / ($eig[0] - $sxx);
	my $gamma = ($sxy * $beta + $syz) / ($eig[0] - $syy - $sxy * $alpha);
	$eig[1] = $alpha * $gamma + $beta;
	$eig[2] = $gamma;
	$eig[3] = 1;
	my $dd = 0;
	for(my $i = 1; $i <= 3; $i++)
	{
	    $dd += $eig[$i] * $eig[$i];
	}
	my $dd = 1.0 / sqrt($dd);
	for(my $i = 1; $i <= 3; $i++)
	{
	    $eig[$i] *= $dd;
	}
    }
    return @eig;
}
