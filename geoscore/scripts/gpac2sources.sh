#!/bin/bash

#ARGS (1) file
if [ $# -ne 1 ]; then
   echo "usage: ARGS (1) GPAC standard output file"
   exit 1
fi 

#check to make sure the file exists
if [[ ! -f $1 ]]; then
    echo "Cannot find $1\n";
    exit 1
fi

grep ^source $1 | awk '{split($3,a,"=");t=a[2]; split($4,a,"=");x=a[2]; split($5,a,"=");y=a[2]; split($6,a,"=");z=a[2]; split($7,a,"=");m=a[2];mag=(log(m)/log(10)-9.1)/1.5;print x" "y" "z" "mag" "t;}' 
