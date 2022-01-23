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
git push $REMOTE $BRANCH
git fetch $REMOTE
git checkout master
git pull $REMOTE master
git merge --squash origin/$BRANCH
if [ $? -eq 0 ]; then
    git commit -m "$1"
    if [ $? -eq 0 ]; then
	git push $REMOTE master
	git branch -d $BRANCH
    else
	echo "Commit failed: $?";	
    fi
else
    echo "Merge failed: $?";
fi
