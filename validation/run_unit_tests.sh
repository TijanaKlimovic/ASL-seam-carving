#!/bin/bash

# build c small unit tests
gcc main_test.c ../convolution.c ../parse_img.c ../min_seam.c ../optimal_image.c -o small_tests -lm
echo "Running small unit tests"
./small_tests

# # build automatized c tests
# gcc tests.c ../convolution.c ../parse_img.c ../min_seam.c ../optimal_image.c -o tests -lm

# for img in unit_tests/input_automated/*.png
# do
# 	echo "Processing $img"

# 	# run c tests
# 	./tests $img

# 	# run python tests
# 	python3 test.py $img

# 	# go over written out files and compare them
# 	for out in unit_tests/out_c/*.txt
# 	do
# 		HASH=$(sha256sum $out | awk '{print $1}')
# 		FILENAME=$(echo "${out##*/}")
# 		REFERENCE_FILE="unit_tests/out_reference/$FILENAME"
# 		echo "$HASH  $REFERENCE_FILE" | sha256sum --check
# 	done

# 	echo ""

# done
