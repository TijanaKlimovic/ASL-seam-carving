//#error Please comment out the next two lines under linux, then comment this error
//#include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler
//#include <windows.h> // Include if under windows

#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "tsc_x86.h"

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
//#define FREQUENCY 3.6e9
#define CALIBRATE

//--------------- necessary for benchmarking---------------

#include <stdio.h>
#include "parse_img.h"
#include "min_seam.h"
#include "optimal_image.h"

//--------------- necessary for run---------------

#include <string.h>

//--------------- extern count variables --------------

#include "count.h"

#ifdef count_instr 
unsigned long long add_count = 0;	//count the total number of add instructions
unsigned long long mult_count = 0; //count the total number of mult instructions
#endif

//---------------------------------------------------------

void init_tsc() {
	; // no need to initialize anything for x86
}

myInt64 start_tsc(void) {
    tsc_counter start;
    CPUID();
    RDTSC(start);
    return COUNTER_VAL(start);
}

myInt64 stop_tsc(myInt64 start) {
	tsc_counter end;
	RDTSC(end);
	CPUID();
	return COUNTER_VAL(end) - start;
}

/* 
takes in the image path and runs the seam carving algorithm
conforms to the reference python interface
*/

int run_python_validation(const char* path,const char* output_file_name, const char* col_row,
	double percentage) {
	int width, height;
	int *output;
	int *res;

	if (!load_image(path, &width, &height, &output)) {
		return 1;
	}

	int width_diff = 0, height_diff = 0;

	if (strcmp(col_row, "c")) {
		height_diff = height - floor(height * (0.01*percentage));
	} else {
		width_diff = width - floor(width * (0.01*percentage));
	}

	res = optimal_image(width, height, width_diff, height_diff, output);
	save_image(output_file_name, width - width_diff, height - height_diff, res);
	free(res);
	return 0;
}

/* 
takes in the image path and resuns the seam carving algorithm
*/

// Data structure to hold information in a T cell
struct cell_T {
	// optimal cost to get to this dimension
	int optimal_cost;
	// corresponding matrix for this dimension,
	// obtained based on optimal cost
	// size is 3 * (width - row) * (height - column)
	int *i;
};

int run(int width, int height, int *image, const char *output_file_name, int width_diff, int height_diff ) {
	width_diff++;
	height_diff++;

	struct cell_T *T = (struct cell_T *)malloc(width_diff * height_diff * sizeof(struct cell_T));
	T[0].optimal_cost = 0;
	T[0].i = image;

	// calc first row -> vertical seam only
	int j;
	for (j = 1; j < width_diff; ++j) {
		int T_index_left = j - 1;
		int *image_left = T[T_index_left].i;
		int image_width = width - j + 1;

		int *backtrack = (int *)malloc(height * sizeof(int));
		int optimal_cost = min_seam(height, image_width, image_left, 1, backtrack);

		T[j].optimal_cost = T[T_index_left].optimal_cost + optimal_cost;
		T[j].i = (int *)malloc(3 * height * (image_width - 1) * sizeof(int));

		int k, l;
		for (k = 0; k < height; ++k) { // construct each row at a time
			int crr_col = 0;
			for (l = 0; l < image_width; ++l)
				if (l != backtrack[k]) { // check for elem to remove
					T[j].i[k*(image_width-1)+crr_col] = T[T_index_left].i[k*image_width + l];
					T[j].i[height*(image_width-1)+k*(image_width-1)+crr_col] = T[T_index_left].i[height*image_width+k*image_width + l];
					T[j].i[2*height*(image_width-1)+k*(image_width-1)+crr_col] = T[T_index_left].i[2*height*image_width+k*image_width + l];
					crr_col++;
				}
		}
		free(backtrack);
	}
		
	// going through the rest of the dp matrix
	int i;
	for (i = 1; i < height_diff; ++i) {
		// first column -> horizontal seam only
		int T_index = i * width_diff;
		int T_index_up = (i - 1) * width_diff;
		int *image_up = T[T_index_up].i;
		int image_height = height - i + 1;

		int *backtrack = (int *)malloc(width * sizeof(int));
		int optimal_cost = min_seam(image_height, width, image_up, 0, backtrack);

		T[T_index].optimal_cost = T[T_index_up].optimal_cost + optimal_cost;
		T[T_index].i = (int *)malloc(3 * (image_height - 1) * width * sizeof(int));

		int k, l;
		for (k = 0; k < width; ++k) { // construct each column at a time
			int crr_row = 0;
			for (l = 0; l < image_height; ++l)
				if (l != backtrack[k]) { // check for elem to remove
					T[T_index].i[crr_row*width+k] = T[T_index_up].i[l*width + k];
					T[T_index].i[(image_height-1)*width+crr_row*width+k] = T[T_index_up].i[image_height*width+l*width + k];
					T[T_index].i[2*(image_height-1)*width+crr_row*width+k] = T[T_index_up].i[2*image_height*width+l*width + k];
					crr_row++;
				}
		}
		free(backtrack);

		// continue from the 2nd column
		for (j = 1; j < width_diff; ++j) {
			// height -> horizontal seam -> row
			// width -> vertical seam -> column
			int T_index = i * width_diff + j;
			int T_index_left = i * width_diff + j - 1;
			int *image_left = T[T_index_left].i;
			int image_left_width = width - j + 1;
			int image_left_height = height - i;

			int *backtrack_left = (int *)malloc(image_left_height * sizeof(int));
			int optimal_cost_left = min_seam(image_left_height, image_left_width, image_left, 1, backtrack_left);

			int T_index_up = (i - 1) * width_diff + j;
			int *image_up = T[T_index_up].i;
			int image_up_width = width - j;
			int image_up_height = height - i + 1;

			int *backtrack_up = (int *)malloc(image_up_width * sizeof(int));
			int optimal_cost_up = min_seam(image_up_height, image_up_width, image_up, 0, backtrack_up);

			// remove column
			if (T[T_index_left].optimal_cost + optimal_cost_left <=
				T[T_index_up].optimal_cost + optimal_cost_up) {
				T[T_index].optimal_cost = T[T_index_left].optimal_cost + optimal_cost_left;
				T[T_index].i = (int *)malloc(3 * image_left_height * (image_left_width - 1) * sizeof(int));
				
				int k, l;
				for (k = 0; k < image_left_height; ++k) { // construct each row at a time
					int crr_col = 0;
					for (l = 0; l < image_left_width; ++l)
						if (l != backtrack_left[k]) { // check for elem to remove
							T[T_index].i[k*(image_left_width-1)+crr_col] = T[T_index_left].i[k*image_left_width + l];
							T[T_index].i[image_left_height*(image_left_width-1)+k*(image_left_width-1)+crr_col] = T[T_index_left].i[image_left_height*(image_left_width)+k*image_left_width + l];
							T[T_index].i[2*image_left_height*(image_left_width-1)+k*(image_left_width-1)+crr_col] = T[T_index_left].i[2*image_left_height*(image_left_width)+k*image_left_width + l];
							crr_col++;
						}
				}

			// remove row
			} else {
				T[T_index].optimal_cost = T[T_index_up].optimal_cost + optimal_cost_up;
				T[T_index].i = (int *)malloc(3 * (image_up_height - 1) * image_up_width * sizeof(int));

				int k, l;
				for (k = 0; k < image_up_width; ++k) { // construct each column at a time
					int crr_row = 0;
					for (l = 0; l < image_up_height; ++l)
						if (l != backtrack_up[k]) { // check for elem to remove
							T[T_index].i[crr_row*image_up_width+k] = T[T_index_up].i[l*image_up_width + k];
							T[T_index].i[(image_up_height-1)*image_up_width+crr_row*image_up_width+k] = T[T_index_up].i[image_up_height*image_up_width+l*image_up_width + k];
							T[T_index].i[2*(image_up_height-1)*image_up_width+crr_row*image_up_width+k] = T[T_index_up].i[2*image_up_height*image_up_width+l*image_up_width + k];
							crr_row++;
						}
				}
			}
			free(backtrack_left);
			free(backtrack_up);
		
		}
	}
		
	// copy 
	int *res = malloc(3 * (width - width_diff + 1) * (height - height_diff + 1) * sizeof(int));
	memcpy(res, T[width_diff * height_diff - 1].i, 3 * (width - width_diff + 1) * (height - height_diff + 1) * sizeof(int));
	for (int i = 1; i < width_diff * height_diff; ++i) {
		free(T[i].i);
	}

	free(T);

	save_image(output_file_name, width - width_diff + 1, height - height_diff + 1, res);
	free(res);
	return 0;
}

/* 
benchmarking function using rdtsc instruction
*/

double rdtsc(int width, int height, int *output, const char *output_file_name, int width_diff, int height_diff) {
    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;


#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i) {
            run(width, height, output, output_file_name, width_diff, height_diff);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
		run(width, height, output, output_file_name, width_diff, height_diff);	
	}
    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}

/* main function , has two options, first is the image path second is -t/<anything except nothing> */ 

int main(int argc, char const *argv[]) {

	if (argc == 5) {
		// verification run - python interface
		printf("Usage should be: %s <image_path> <output_file_name> <c/r> <percentage>\n", argv[0]);

		double percentage = atoi(argv[4]);
		int to_retun = run_python_validation(argv[1], argv[2], argv[3], percentage);
		#ifdef count_instr 
		printf("\nADDS=%llu MULTS=%llu\n", add_count, mult_count);
		#endif
		return to_retun;
	}

	if (argc < 6) {
		printf("Usage: %s <image_path> <output_file_name> <width_diff> <height_diff> <timing boolean>\n", argv[0]);
		return 1;
	}
	// normal run

	int width_diff = atoi(argv[3]);
	int height_diff = atoi(argv[4]);
	int width, height;
	int *output;

	if (!load_image(argv[1], &width, &height, &output)) {
		return 1;
	}

	if(strcmp(argv[5],"1")){
		int out = run(width, height, output, argv[2], width_diff, height_diff);
		free(output);
		return out;
	}
	//time it
	double r = rdtsc(width, height, output, argv[2], width_diff, height_diff);

	free(output);
    printf("RDTSC instruction: %.0lf cycles measured\n", r);
    return 0;
}
