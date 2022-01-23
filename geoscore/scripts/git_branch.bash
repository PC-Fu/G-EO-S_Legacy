#!/bin/bash

#ARGS (1) input file
if [ $# -ne 1 ]; then
   echo "usage: <branch name>\n"
   exit 1
fi 

git checkout -b $1 origin/$1
