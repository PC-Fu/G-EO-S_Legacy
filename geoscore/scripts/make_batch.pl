#!/usr/bin/perl
use strict;
use Getopt::Long;
my $help = "";
my $gpac = "~/gpac_trunk/src/GPAC.x";
my $input = "test.xml";
my $jobName = "gpac_run";
my $procs = -1;
my $time = -1;

GetOptions("exec:s" => \$gpac,
	   "i:s" => \$input,
	   "jobname:s" => \$jobName,
	   "h:s" => \$help,
	   "t:i" => \$time);

die "usage: $0 -exec /path/to/gpac -i /path/to/input.xml -jobname simulation_name -h help -t time_in_seconds_to_run\n" if($help ne "");

#check that the executable exists
if(-e $gpac)
{
}
else
{
    die "You must specify the path to your executable with -exec\n";
}

#GET USER
my $user = "";
{
    open(IN,"whoami|");
    $user = <IN>;
    chomp($user);
    $user =~s/^\s+//;
}

#GET PARALLEL
if(-e $input)
{
    die "Cannot open $input" if(!open(IN,"<$input"));
    $procs = 1;
    while(<IN>)
    {
	my $line = $_;
	if(/xpar\s*\=\s*\"(\S+)\"/)
	{
	    #print "xpar=$1\n";
	    $procs *= $1;
	}
	if(/ypar\s*\=\s*\"(\S+)\"/)
	{
	    #print "ypar=$1\n";
	    $procs *= $1;
	}
	if(/zpar\s*\=\s*\"(\S+)\"/)
	{
	    #print "zpar=$1\n";
	    $procs *= $1;
	}
    }
    close IN;
}
else
{
    die "The input file \"$input\" does not exist";
}
die "Too few processes requested: $procs" if $procs <= 0;

#GET LOCATION
my $cwd = "";
{
    open(IN,"pwd|");
    $cwd = <IN>;
    close IN;
    $cwd =~s/^\s+//;
    chomp($cwd);
}

#GET MACHINE
my $cluster = "";
{
    open(IN,"hostname|");
    while(<IN>)
    {
	my @arr = split(/\d+/,$_);
	$cluster = $arr[0];
    }
    close IN;
}

#GET RESOURCES
my $nodesMax = 1;
my $hoursMax = 10;
my $procsPerNode = 8;
{
    if($cluster eq "hera")
    {
	$procsPerNode = 16;
    }
    elsif($cluster eq "aztec")
    {
	$procsPerNode = 12;
    }
    if($cluster eq "ansel")
    {
	$nodesMax = 64;
	$hoursMax = 12;
    }
    else
    {
	open(IN,"news job.lim.$cluster|");
	while(<IN>)
	{
	    if(/pbatch\s+\S+\s+weekdays\s+(\S+)\s+(\d+)/)
	    {
		$nodesMax = $1;
		$hoursMax = $2;
	    }
	}
	close IN;
    }
}

#GET BANK
my $bank = "";
{
    open(IN,"mdiag -u $user|");
    while(<IN>)
    {
	if(/\=(\S+)\s*$/)
	{
	    my @arr = split(/\,/,$1);
	    my $nbanks = $#arr + 1;
	    if($nbanks <= 1)
	    {
		if($nbanks == 0)
		{
		    die "sorry, you don't have any banks!";
		}
		else
		{
		    $bank = $arr[0];
		}
	    }
	    else
	    {
		my $rr = rand();
		my $index = $#arr;
		if($rr > 0.5)
		{
		    $index++;
		}
		$bank = $arr[$index];
	    }
	    #die "$bank";
	}
    }
    close IN;
}

#SET VARIABLES
my $sMax = $hoursMax * 60 * 60;
if($cluster=~m/aztec/)
{
    $sMax = 360000;
}
if($time < 0 || $time > $sMax)
{
    $time = $sMax;
}

my $nodes = int($procs / $procsPerNode);
if($procs > $nodes * $procsPerNode)
{
    $nodes++;
}
if($nodes > $nodesMax)
{
    die "You asked for too many nodes: $nodes > $nodesMax\n";
}

my $pname = "nodes";
if($cluster=~m/aztec/)
{
    $pname = "ttc";
    $nodes = $procs;
    $bank = "ees";
}

print "#!/bin/csh
#MSUB -N $jobName
#MSUB -l walltime=$time
#MSUB -l $pname=$nodes
#MSUB -l partition=$cluster
#MSUB -q pbatch
#MSUB -A $bank
#MSUB -r n
#MSUB -m be

date > jstart
cd $cwd
";
if($procs == 1)
{
    print "  $gpac -i $input > prunout\n";
}
else
{
    print "  srun -n $procs -N $nodes $gpac -i $input > prunout\n";
}
print "date > jend
";
