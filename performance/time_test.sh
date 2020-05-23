#!/bin/bash

# Usage: ./seam_test.sh <image> <lower_bound%> <upper_bound%> <step%>
# Resizes image from 'lower_bound' width (keeping aspect ratio) to 'upper_bound' width by steps of 'step'.
# At each resized image removes 'seams_to_remove' seams.

if [ "$#" -ne 8 ]
then
    echo "Usage: ./time_test.sh <img_folder> <lower_bound_img_width> <upper_bound_img_width> <step> <seams_to_remove> <run_id> <vec_bool> <src_folder>"
    exit
fi

dir="$1"
lower_bound="$2"
upper_bound="$3"
step="$4"
seams="$5"
run_id="$6"
vec_bool="$7"
src_folder="$8"

if [ ! -d "$1" ]
then
    echo "$1 doesn't exist"
    echo "Usage: ./time_test.sh <img_folder> <lower_bound_img_width> <upper_bound_img_width> <step> <seams_to_remove> <run_id> <vec_bool> <src_folder>"
    exit
fi

if [ ! -d "$src_folder" ]
then
    echo "$src_folder doesn't exist"
    echo "Usage: ./time_test.sh <img_folder> <lower_bound_img_width> <upper_bound_img_width> <step> <seams_to_remove> <run_id> <vec_bool> <src_folder>"
    exit
fi

if [ -d "$run_id" ]
then
    rm -rf "$run_id" #remove old workspace
fi

mkdir "$run_id" #create new workspace

#compile code with optimization
if [ "$vec_bool" -eq 0 ]
then
    echo "gcc -Ofast -fno-tree-vectorize $src_folder/*.c -o $run_id/seam_carving -lm"
    gcc -Ofast -fno-tree-vectorize "$src_folder"/*.c -o "$run_id"/seam_carving -lm
else
    echo "gcc -Ofast -march=native $src_folder/*.c -o $run_id/seam_carving -lm"
    gcc -Ofast -march=native "$src_folder"/*.c -o "$run_id"/seam_carving -lm
fi

#compile code with no optimization
echo "gcc $src_folder/*.c -o $run_id/seam_carving_noop -lm"
gcc "$src_folder"/*.c -o "$run_id"/seam_carving_noop -lm
# compile code with instrumentation
echo "gcc -D count_instr -Ofast $src_folder/*.c -o $run_id/seam_carving_ctr -lm"
gcc -D count_instr -Ofast "$src_folder"/*.c -o "$run_id"/seam_carving_ctr -lm

mkdir -p "$run_id/resized" #it puts here the resized images
mkdir -p "$run_id/out" # it puts here the seam carved resized images -> if we don't have timing flag switched on
res_file="$run_id/result.txt"

echo "Running seam carving:"

for img in "$dir"/*.png
do
    fname=$(basename "$img")
    echo "$fname" >> "$res_file"
    echo -e "-O0\t \t \t \t-O3" >> "$res_file"
    echo -e "Img Width\tCycles/Time\tFlops\tPerformance [Flops/Cycle]\tImg Width\tCycles/Time\tFlops\tPerformance [Flops/Cycle]" >> "$res_file"

    for size in $(seq "$lower_bound" "$step" "$upper_bound")
    do 
        resized="$run_id/resized/resized_${size}.png"
        convert "$img" -resize "$size" "$resized"
        out="$run_id/out/${size}_${fname}"
        # run it without timing but with instrumentation
        echo "$run_id/seam_carving_ctr $resized $out $seams $seams 0"
        # cnt=$("$run_id"/seam_carving_ctr "$resized" "$out" "$seams" "$seams" 0 | tr "=ADSULT" " " | tr "M" "+" | bc)
        cnt=$("$run_id"/seam_carving_ctr "$resized" "$out" "$seams" "$seams" 0 | tail -n 1 | awk -F "=" '{print $2}')
        # run it with timing -O3
        echo "$run_id/seam_carving $resized $out $seams $seams 1"
        o3time=$("$run_id"/seam_carving "$resized" "$out" "$seams" "$seams" 1 | tail -n 1 | awk '{print $3}')
        #run it with timing -O0
        echo "$run_id/seam_carving_noop $resized $out $seams $seams 1"
        o0time=$($run_id/seam_carving_noop "$resized" "$out" "$seams" "$seams" 1 | tail -n 1 | awk '{print $3}')
        #performances
        o3perf=$(echo "scale=4; $cnt / $o3time" | bc)
        o0perf=$(echo "scale=4; $cnt / $o0time" | bc)

        echo -e "${size}\t${o0time}\t${cnt}\t${o0perf}\t${size}\t${o3time}\t${cnt}\t${o3perf}" >> "$res_file"
    done

    echo "" >> "$res_file"
done

