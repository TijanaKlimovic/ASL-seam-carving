#!/bin/bash

cd validation/out/$1$2p
touch ../$1$2diff/results

for i in *.png
do
    compare -verbose -metric MAE $i ../$1$2/$i ../$1$2diff/$i &>> ../$1$2diff/results
done