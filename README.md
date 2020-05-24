# asl_seam_carving

Optimizing the hell out of seam carving.

## Running
Compile: gcc -Wall \*.c -o seam_carving -lm

Usage: ./seam_carving \<input_file_path\> \<output_file_name\> \<with_diff\> \<height_diff\>

Expects: a PNG image with 3 channels (without alpha channel)

## Project outline
The project contains two approach paths: 

* track 1: RGB as a 3-dimentional matrix with the channel as the first dimention, scalar replacement, unrolling, restructured computations, blocking
* track 2: RGB as 2-dimentional matrix of pair elements (r,g,b), scalar replacement, unrolling, restructured computations and if-else blocks, SIMD, blocking, inlining

The outlne of the project:

* track 1
	1. branch t1_1 - base
	2. branch t1_2 - optimized computations (convolution computations)
	3. branch t1_3 - scalar replacement, unrolling
	4. branch t1_4 - blocking
	5. branch t1_5 - improved blocking at higher stride
* track 2
	1. branch t1_1 - base
	2. branch t1_2 - optimal omputations
	3. branch t2_3 - RGB representation with pair elements (r,g,b)
	4. branch t2_3_2 - optimized computations (image rotation to find min seam)
	5. branch t_2_4 - SIMD
	6. branch t2_5 - blocking
	7. branch t2_5_sr - scalar replacement
	8. branch t2_inlined - inlining

The project consists of 3 main parts:

* importance / pixel: energy map for image constructed with convolution
* min seam identification: dynamic programming for wanted direction (horizontal / vertical)
* optimal removal order for min seams: ynamic programming for wanted number of seams to remove