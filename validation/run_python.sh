#!/bin/bash

cd validation/no_alpha_images

for i in *.png
do
	python3 ../../seam_carving.py $1 $2 $i ../out/$1$3p/$i
done