#!/usr/bin/perl

my $narg = $#ARGV + 1;

if($narg < 1)
{
    print "usage: $0 <filenames ...>\n";
}    

my @min = (1e100,1e100,1e100);
my @max = (-1e100,-1e100,-1e100);

for(my $i=0;$i<$narg;$i++)
{
    my $filename = shift;
    open(IN,"<$filename");
    while(<IN>)
    {
	my $line = $_;
	chomp($line);
	my @arr = split(/ /,$line);
	for(my $j=0;$j<3;$j++)
	{
	    if($max[$j] < $arr[$j])
	    {
		$max[$j] = $arr[$j];
	    }
	    if($min[$j] > $arr[$j])
	    {
		$min[$j] = $arr[$j];
	    }
	}
    }
    close IN;
}

for(my $k=0;$k<3;$k++)
{
    print "$k: $min[$k] $max[$k]\n";
}
