#!/bin/bash

if [ "$#" -ne 2 ]
then
	echo "usage: ./br.sh <img_file> <src_dir>"
	exit
fi

img="$1"

here=$(pwd)

gcc "$2"/*.c -Ofast -Wall -lm -g -o seam_carving

width=$(identify "$img" | awk '{print $3}' | tr "x" "\n" | head -n 1)
height=$(identify "$img" | awk '{print $3}' | tr "x" "\n" | tail -n 1)

i="$width"
fname=$(basename "$img")
mkdir inputs
mkdir outputs

while [ "$i" -gt 0 ]
do
	new_name="$i"x"$height""_$fname"
	convert "$img" -resize "$i"x"$height"\! "inputs/$new_name"
	dim_diff=$(echo "$i * 0.15" | bc | tr "." "\n" | head -n 1)
	echo "./seam_carving inputs/$new_name outputs/$new_name $dim_diff 0 1"
	result=$(./seam_carving "inputs/$new_name" "outputs/$new_name" "$dim_diff" 0 1 | awk '{print $3}') #time
	(( i = i - 200 ))
done
