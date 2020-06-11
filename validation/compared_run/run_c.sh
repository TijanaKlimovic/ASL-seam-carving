#!/bin/bash

# Usage: ./run_c.sh

cd ../../
for b in "t1_1" "t1_5" "t2_5"
do
	git checkout $b
	gcc *.c -O3 -march=native -lm -g -o main
	cd validation/compared_run/resized

	for i in *.png
	do
		../../../main $i ../$b/$i 20 20 0
	done

	cd ../../../
done