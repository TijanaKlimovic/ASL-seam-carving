#!/bin/bash

# Usage: ./seam_test.sh <image> <lower_bound%> <upper_bound%> <step%>
# Resizes image from 'lower_bound' width (keeping aspect ratio) to 'upper_bound' width by steps of 'step'.
# At each resized image removes 'seams_to_remove' seams.

if [ "$#" -ne 5 ]
then
    echo "Usage: ./time_test.sh <img_file> <lower_bound_img_width> <upper_bound_img_width> <step> <seams_to_remove>"
    exit
fi

# compile code without instrumentation
gcc -O3 -fno-tree-vectorize ../*.c -o seam_carving -lm
# compile code with instrumentation
gcc -D count_instr -O3 -fno-tree-vectorize ../*.c -o seam_carving_ctr -lm

img=$1
lower_bound=$2
upper_bound=$3
step=$4
seams=$5

filename=$(echo "${img##*/}")

# clean up outputs from prev runs
if [ -d "resized" ]
then
    rm -f resized/*
fi

if [ -d "out" ]
then
    rm -f out/*
fi

rm -f cycles.txt
rm -f counts.txt

mkdir -p resized #it puts here the resized images
mkdir -p out # it puts here the seam carved resized images -> if we don't have timing flag switched on

echo "Running seam carving:"

for size in $(seq $lower_bound $step $upper_bound)
  do 
     resized="resized/resized_${size}_${seams_to_remove}.png"
     convert $img -resize $size $resized
     out="out/${size}_${filename}"
     # run it with timing
     echo "./seam_carving $resized $out $seams $seams 1"
     ./seam_carving $resized $out $seams $seams 1 >> cycles.txt
     # run it without timing but with instrumentation
     echo "./seam_carving_ctr $resized $out $seams $seams 0"
     ./seam_carving_ctr $resized $out $seams $seams 0 >> counts.txt
 done


# call python script to parse results
python3 parse_out.py
