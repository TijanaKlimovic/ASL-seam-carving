#!/bin/bash

gcc *.c -lm -g -o main
cd validation/no_alpha_images

for i in *.png
do
	../../main $i ../out/$1$2/$i $1 $2
done