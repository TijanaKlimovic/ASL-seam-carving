#!/bin/bash

gcc ../main.c ../convolution.c ../parse_img.c ../min_seam.c ../optimal_image.c -o main -lm

for i in input/*.png
do
	# Usage: %s <image_path> <output_file_name> <width_diff> <height_diff> <timing boolean>
	fname=${i##*/}
	./main $i "output/5_10_${fname}" 5 10 0
	./main $i "output/50_0_${fname}" 50 0 0
	./main $i "output/0_50_${fname}" 0 50 0

done

# go over written out files and compare them
for out in output/*.png
do
	HASH=$(sha256sum $out | awk '{print $1}')
	FILENAME=$(echo "${out##*/}")
	REFERENCE_FILE="reference_output/$FILENAME"
	echo "$HASH  $REFERENCE_FILE" | sha256sum --check
done
