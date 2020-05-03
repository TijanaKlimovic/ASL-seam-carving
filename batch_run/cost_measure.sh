#!/bin/bash

#PARAMS [EDITABLE]
#--------------------------------------------------------

#--------------------------------------------------------
#END OF PARAMS [DON'T EDIT THE REST OF THE CODE]

tmp_file="cost_measure_tmp_file.tmp"


if [ "$#" -ne 2 ]
then
	echo "usage: ./br.sh <input_image_dir> <br_file>"
	exit
fi

if [ ! -d "$1" ]
then
	echo "ERROR: not a directory: $1"
	echo "usage: ./br.sh <input_image_dir> <br_file>"
	exit
fi

if [ ! -f "$2" ]
then
	echo "ERROR: not a file: $2"
	echo "usage: ./br.sh <input_image_dir> <br_file>"
	exit
fi

check_str=$(awk '/cost_measure/' "$2")

if [ -n "$check_str" ]
then
	echo "ERROR: $2 already contains the cost_measure column"
	exit
fi

# printf "%20s %15s\n" "cost_measure" "performance" > "$tmp_file"

list=($(awk '{print $1}' "$2"))
size="${#list[@]}"
idx=2 #idx = 1 is headers

while (( idx <= size ))
do
	name=$(awk "NR == $idx {print \$1}" "$2")
	width=$(awk "NR == $idx {print \$2}" "$2")
	height=$(awk "NR == $idx {print \$3}" "$2")
	m_diff=$(awk "NR == $idx {print \$4}" "$2")
	n_diff=$(awk "NR == $idx {print \$5}" "$2")
	cycles=$(awk "NR == $idx {print \$6}" "$2")

	input_name="$width"x"$height"_"$name"
	inp_f="$1/$input_name"

	#number of columns, aka width
	m=$(identify "$inp_f" | awk '{print $3}' | tr "x" "\n" | head -n 1)

	#number of rows, aka height
	n=$(identify "$inp_f" | awk '{print $3}' | tr "x" "\n" | tail -n 1)

	opt_seam=$(echo "$m_diff + $n_diff + 3 * $m_diff * $n_diff" | bc) #adds in optimal seam
	# echo "opt_seam = $opt_seam"

	if [ "$m_diff" -lt 2 ]
	then
		dora=0
	else
		dora=$(echo "($m_diff * $m - ($m_diff - 1) * ($m_diff - 2) / 2) * (115 * $n - 1)" | bc)
	fi

	# echo "dora = $dora"

	if [ "$n_diff" -lt 2 ]
	then
		ioana=0
	else
		ioana=$(echo "($n_diff * $n - ($n_diff - 1) * ($n_diff - 2) / 2) * (115 * $m - 1)" | bc)
	fi

	# echo "ioana = $ioana"
	aydin=$(echo "($n_diff * $m_diff * (2 * $m - $m_diff + 1) / 2) * (115 * (-$n_diff - 1 + 2 * $n) / 2 - 1)" | bc) #aydin's part
	# echo "aydin = $aydin"
	tijana=$(echo "($m_diff * $n_diff * (2 * $n - $n_diff + 1) / 2) * (115 * (-$m_diff - 1 + 2 * $m) / 2 - 1)" | bc) #tijana's part
	# echo "tijana = $tijana"
	cost=$(echo "$opt_seam + $dora + $ioana + $aydin + $tijana" | bc)
	performance=$(echo "scale=4; $cost / $cycles" | bc) #flops / cycle

	# printf "%20d %15s\n" "$cost" "$performance" >> "$tmp_file"
	printf "%-20s %20s\n" "$name" "$performance"

	(( idx = idx + 1 ))
done

# paste "$2" "$tmp_file"
# rm "$tmp_file"
