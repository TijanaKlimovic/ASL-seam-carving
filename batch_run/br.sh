#!/bin/bash

#PARAMS [EDITABLE]
#--------------------------------------------------------
w_red_p=0.01 #width reduction percentage
h_red_p=0.02 #height reduction percentage

#these 2 steps are performed in a single for (not nested)
w_step=100 #width step (for iteration); w = w - w_step
h_step=200 #height step (for iteration); h = h - h_step

gcc_flags="-Ofast -Wall -lm -g"

timing=0 #passed to ./seam_carving
#--------------------------------------------------------
#END OF PARAMS [DON'T EDIT THE REST OF THE CODE]


if [ "$#" -ne 3 ]
then
	echo "usage: ./br.sh <img_file/img_folder> <src_dir> <output_file.txt>"
	exit
fi

if [ -f "$1" ]
then
	# single image passed
	files=("$1")
	f_step=0
	width=$(identify "$1" | awk '{print $3}' | tr "x" "\n" | head -n 1)
	height=$(identify "$1" | awk '{print $3}' | tr "x" "\n" | tail -n 1)
elif [ -d "$1" ]
then
	# directory passed
	files=()
	width=99999999
	height="$width"

	for f in "$1"/*.png
	do
		files+=("$f")
		tmp_w=$(identify "$f" | awk '{print $3}' | tr "x" "\n" | head -n 1)
		tmp_h=$(identify "$f" | awk '{print $3}' | tr "x" "\n" | tail -n 1)

		if [ "$tmp_w" -lt "$width" ]
		then
			width=$tmp_w
		fi

		if [ "$tmp_h" -lt "$height" ]
		then
			height=$tmp_h
		fi

	done

	f_step=1
else
	echo "$1 doesn't exist"
	echo "usage: ./br.sh <img_file/img_folder> <src_dir> <output_file.txt>"
	exit
fi

file_count=${#files[@]}

printf "%-10s %10s %10s %10s %10s %20s\n" "filename" "width" "height" "w_diff" "h_diff" "time[cycles]" >> "$3"

gcc "$2"/*.c $gcc_flags -o seam_carving

w="$width"
h="$height"
f_idx=0
mkdir inputs
mkdir outputs

while (( w > 0 && h > 0 && f_idx < file_count ))
do
	img="${files[f_idx]}"
	fname=$(basename "$img")
	new_name="$w"x"$h""_$fname"
	convert "$img" -resize "$w"x"$h"\! "inputs/$new_name"
	w_red=$(echo "$w * $w_red_p" | bc | tr "." "\n" | head -n 1)
	h_red=$(echo "$h * $h_red_p" | bc | tr "." "\n" | head -n 1)
	echo "./seam_carving inputs/$new_name outputs/$new_name $w_red $h_red $timing"
	result=$(./seam_carving "inputs/$new_name" "outputs/$new_name" "$w_red" "$h_red" "$timing" | awk '{print $3}') #time
	printf "%-10s %10d %10d %10d %10d %20d\n" "$fname" "$w" "$h" "$w_red" "$h_red" "$result" >> "$3"

	if [ "$f_step" -eq 0 ]
	then
		(( w = w - w_step ))
		(( h = h - h_step ))
	fi

	(( f_idx = f_idx + f_step ))
done
