#!/bin/bash

# Usage: ./run_convert.sh

for i in ref_png/*.png
do
    convert $i -alpha off no_alpha_$i
done