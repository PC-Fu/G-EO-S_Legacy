#!/bin/bash

#C
for f in `find ./ -name "*.c"`
do
    perl replace_license.pl $f > tmp
    mv tmp $f
done

#H
for f in `find ./ -name "*.h"`
do
    perl replace_license.pl $f > tmp
    mv tmp $f
done

#CPP
for f in `find ./ -name "*.cpp"`
do
    perl replace_license.pl $f > tmp
    mv tmp $f
done
