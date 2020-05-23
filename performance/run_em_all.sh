#!/bin/bash

if [ -f "time_test.sh" ]
then
	echo "Bad File!"
	echo "Please follow these instructions: (For each run)"
	echo "1. Put your team046 repository into a clean state; i.e. you should not have any 'Changes not staged for commit'"
	echo "2. Create a folder (outside the whole team046 repository folder). This folder will contain the results."
	echo "3. Copy this file (run_em_all.sh) to that folder."
	echo "4. Run the script from that folder."
	echo "5. Don't forget to turn off your CPU's Turbo Boost feature."
	exit
fi

if [ "$#" -ne 5 ]
then
    echo "Usage: ./run_em_all.sh <repo_folder> <lower_bound_img_width> <upper_bound_img_width> <step> <seams_to_remove>"
    exit
fi

myhome=$(pwd)

repo="$1"
lower_bound="$2"
upper_bound="$3"
step="$4"
seams="$5"

if [ ! -d "$repo" ]
then
	echo "$repo doesn't exist"
    echo "Usage: ./run_em_all.sh <repo_folder> <lower_bound_img_width> <upper_bound_img_width> <step> <seams_to_remove>"
    exit
fi

cd "$repo"
git checkout master
git pull
cp "performance/time_test.sh" "$myhome/"
cp -r "performance/imgs" "$myhome/"
cd "$myhome"

branches=("t1_1" "t1_2" "t1_3" "t1_4" "t1_5" "t2_3" "t2_3_2" "t_2_4" "t2_5" "t2_5_sr")
# instrumnt=("1" "1" "0" "1" "1" "1" "1" "0" "1" "0")
vec_bool=0

for br in ${branches[@]}
do

	if [ "$br" = "t_2_4" ]
	then
		vec_bool=1
	fi

	cd "$repo"
	git checkout "$br"
	git pull
	cd "$myhome"

	for ratio in imgs/*
	do
		id=$(basename "$ratio")
		./time_test.sh "$ratio" "$lower_bound" "$upper_bound" "$step" "$seams" "${br}_${id}" "$vec_bool" "${repo}"
	done

done
