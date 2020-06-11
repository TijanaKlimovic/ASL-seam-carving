#!/bin/bash

# Usage: run_compare.sh <target_dir> <ref_dir>

target="$1"
ref="$2"

cd $target
for i in *.png
do
    compare -verbose -metric MAE $i ../$ref/$i comp/$i &>> comp/results.txt
done