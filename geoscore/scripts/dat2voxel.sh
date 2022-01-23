#!/bin/bash

#ARGS (1) input file
if [ $# -ne 1 ]; then
   echo "usage: ARGS (1) input file"
   exit 1
fi 

tail -n `cat $1 | wc | awk '{print ($1 - 1)}'` $1 > ppcurr

#create tmp perl script
echo 'use strict;' > tmp.pl
echo 'my $fname = shift;' > tmp.pl
echo 'open(IN,"<$fname");' >> tmp.pl
echo 'my $xlast = -1e100;' >> tmp.pl
echo 'while(<IN>)' >> tmp.pl
echo '{' >> tmp.pl
echo '  my $x = 1.0 * $_;' >> tmp.pl
echo '  if($x < $xlast)' >> tmp.pl
echo '  {' >> tmp.pl
echo '    last;' >> tmp.pl
echo '  }' >> tmp.pl
echo '  else' >> tmp.pl
echo '  {' >> tmp.pl
echo '    print $_;' >> tmp.pl
echo '    $xlast = $x;' >> tmp.pl
echo '  }' >> tmp.pl
echo '}' >> tmp.pl
echo 'close IN;' >> tmp.pl

#get xs
cat ppcurr | awk '{print $1}' > xs
perl tmp.pl xs > xpp
rm xs

#get ys
cat ppcurr | awk '{print $2}' > ys
perl tmp.pl ys > yy
uniq yy > ypp
rm ys yy

#get zs
cat ppcurr | awk '{print $3}' > zs
uniq zs > zpp
rm zs