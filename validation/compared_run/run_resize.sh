#!/bin/bash

# Usage: ./run_resize.sh <lower_bound> <upper_bound> <step>

lower_bound="$1"
upper_bound="$2"
step="$3"

for size in $(seq "$lower_bound" "$step" "$upper_bound")
do 
    resized="resized/resized_${size}.png"
    convert 1.png -resize "$size" "$resized"
done