#!/bin/bash


for num_nodes in  1000 5000 10000 15000 20000 30000 40000 50000 100000 250000 500000 1000000
do
    for n in {1..3}
    do
#	echo $num_nodes
	./nodal_update_benchmark.x $num_nodes 1000 1
    done
done


#for n in {1..10}
#do
#./nodal_update_benchmark.x 100000 1000
#done