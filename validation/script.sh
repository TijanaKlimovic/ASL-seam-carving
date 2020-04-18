#!/bin/bash

gcc *.c -lm -g -o main
cd validation

for i in images/*.png
do
	convert $i -resize 250 small_$i
done

for i in small_images/*.png
do
    convert $i -alpha off no_alpha_$i
done

cd no_alpha_small_images

for i in *.png
do
	../../main $i ../out/10x10/$i 10 10
	../../main $i ../out/50x50/$i 50 50
	../../main $i ../out/40x20/$i 40 20
done