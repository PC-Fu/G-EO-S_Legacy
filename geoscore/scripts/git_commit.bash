#!/bin/bash

#ARGS (1) message (2-opt) remote
REMOTE=origin
if [ $# -eq 2 ]; then
    REMOTE=$2
fi

if [ $# -eq 0 ] || [ $# -gt 2 ]; then
   echo "usage: <message> <opt: remote>\n"
   exit 1
fi 

BRANCH=`git status | head -n 1 | awk '{print $4}'`
git commit -a -m "$1"

#only go forward if no error
if [ $? -eq 0 ]; then
    git push $REMOTE $BRANCH
else
    echo "Failed to commit: code $?"
fi
