#!/bin/bash

# Usage: ./run_resize.sh <lower_bound> <upper_bound> <step>

lower_bound="$1"
upper_bound="$2"
step="$3"

for size in $(seq "$lower_bound" "$step" "$upper_bound")
do 
    image="resized_python/resized_${size}.jpg"
    res_image="ref/resized_${size}.jpg"
    new_size=$(($size - 20))
    echo $size
    time python3 CAIS.py -i $image -r $new_size $new_size -o $res_image
done