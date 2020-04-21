#!/bin/bash

cd validation

for i in images/*.png
do
    convert $i -alpha off no_alpha_$i
done