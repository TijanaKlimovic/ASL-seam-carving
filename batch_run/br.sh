#!/bin/bash

#TODO: Compute FLOP based on image

#PARAMS [EDITABLE]
#--------------------------------------------------------
w_red_p=0.00 #width reduction percentage
h_red_p=0.15 #height reduction percentage

#these 2 steps are performed in a single for (not nested)
w_step=0 #width step (for iteration); w = w - w_step
h_step=200 #height step (for iteration); h = h - h_step

gcc_flags="-Ofast -Wall -lm -g"

timing=1 #passed to ./seam_carving
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
		# tmp_w=$(identify "$f" | awk '{print $3}' | tr "x" "\n" | head -n 1)
		# tmp_h=$(identify "$f" | awk '{print $3}' | tr "x" "\n" | tail -n 1)

		# if [ "$tmp_w" -lt "$width" ]
		# then
		# 	width=$tmp_w
		# fi

		# if [ "$tmp_h" -lt "$height" ]
		# then
		# 	height=$tmp_h
		# fi

	done

	width=1000
	height=1000

	f_step=1
else
	echo "$1 doesn't exist"
	echo "usage: ./br.sh <img_file/img_folder> <src_dir> <output_file.txt>"
	exit
fi

file_count="${#files[@]}"

printf "%-10s %10s %10s %10s %10s %20s %20s %15s\n" "filename" "width" "height" "w_diff" "h_diff" "time[cycles]" "cost_measure" "performance" >> "$3"

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
	m_diff=$(echo "$w * $w_red_p" | bc | tr "." "\n" | head -n 1)
	n_diff=$(echo "$h * $h_red_p" | bc | tr "." "\n" | head -n 1)
	echo "./seam_carving inputs/$new_name outputs/$new_name $m_diff $n_diff $timing"
	cycles=$(./seam_carving "inputs/$new_name" "outputs/$new_name" "$m_diff" "$n_diff" "$timing" | awk '{print $3}') #time
	m="$w"
	n="$h"

	opt_seam=$(echo "$m_diff + $n_diff + 3 * $m_diff * $n_diff" | bc) #adds in optimal seam
	echo "opt_seam = $opt_seam"

	if [ "$m_diff" -lt 2 ]
	then
		dora=0
	else
		dora=$(echo "($m_diff * $m - ($m_diff - 1) * ($m_diff - 2) / 2) * (115 * $n - 1)" | bc)
	fi

	echo "dora = $dora"

	if [ "$n_diff" -lt 2 ]
	then
		ioana=0
	else
		ioana=$(echo "($n_diff * $n - ($n_diff - 1) * ($n_diff - 2) / 2) * (115 * $m - 1)" | bc)
	fi

	echo "ioana = $ioana"
	aydin=$(echo "($n_diff * $m_diff * (2 * $m - $m_diff + 1) / 2) * (115 * (-$n_diff - 1 + 2 * $n) / 2 - 1)" | bc) #aydin's part
	echo "aydin = $aydin"
	tijana=$(echo "($m_diff * $n_diff * (2 * $n - $n_diff + 1) / 2) * (115 * (-$m_diff - 1 + 2 * $m) / 2 - 1)" | bc) #tijana's part
	echo "tijana = $tijana"
	cost=$(echo "$opt_seam + $dora + $ioana + $aydin + $tijana" | bc)
	performance=$(echo "scale=4; $cost / $cycles" | bc) #flops / cycle

	printf "%-10s %10d %10d %10d %10d %20d %20d %15s\n" "$fname" "$w" "$h" "$m_diff" "$n_diff" "$cycles" "$cost" "$performance" >> "$3"

	if [ "$f_step" -eq 0 ]
	then
		(( w = w - w_step ))
		(( h = h - h_step ))
	fi

	(( f_idx = f_idx + f_step ))
done
