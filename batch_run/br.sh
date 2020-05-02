#!/bin/bash

#PARAMS [EDITABLE]
#--------------------------------------------------------
w_red_p=0.01 #width reduction percentage
h_red_p=0.02 #height reduction percentage

#these 2 steps are performed in a single for (not nested)
w_step=100 #width step (for iteration); w = w - w_step
h_step=200 #height step (for iteration); h = h - h_step

gcc_flags="-Ofast -Wall -lm -g"

timing=1 #passed to ./seam_carving
#--------------------------------------------------------
#END OF PARAMS [DON'T EDIT THE REST OF THE CODE]


if [ "$#" -ne 3 ]
then
	echo "usage: ./br.sh <img_file> <src_dir> <output_file.txt>"
	exit
fi

img="$1"
fname=$(basename "$img")
printf "%-10s %10s %10s %10s %10s %20s\n" "filename" "width" "height" "w_diff" "h_diff" "time[cycles]" >> "$3"

gcc "$2"/*.c $gcc_flags -o seam_carving

width=$(identify "$img" | awk '{print $3}' | tr "x" "\n" | head -n 1)
height=$(identify "$img" | awk '{print $3}' | tr "x" "\n" | tail -n 1)

w="$width"
h="$height"
mkdir inputs
mkdir outputs

while (( w > 0 && h > 0 ))
do
	new_name="$w"x"$h""_$fname"
	convert "$img" -resize "$w"x"$h"\! "inputs/$new_name"
	w_red=$(echo "$w * $w_red_p" | bc | tr "." "\n" | head -n 1)
	h_red=$(echo "$h * $h_red_p" | bc | tr "." "\n" | head -n 1)
	echo "./seam_carving inputs/$new_name outputs/$new_name $w_red $h_red $timing"
	result=$(./seam_carving "inputs/$new_name" "outputs/$new_name" "$w_red" "$h_red" "$timing" | awk '{print $3}') #time
	printf "%-10s %10d %10d %10d %10d %20d\n" "$fname" "$w" "$h" "$w_red" "$h_red" "$result" >> "$3"
	(( w = w - w_step ))
	(( h = h - h_step ))
done
