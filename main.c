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
#include "optimal_image.h"


//--------------- necessary for run---------------

#include <string.h>

//--------------- extern count variables --------------

#include "count.h"

//#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

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
	unsigned char *output;
	unsigned char *res;

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
takes in the image path and runs the seam carving algorithm
*/

int run(int width, int height, unsigned char *output, const char *output_file_name, int width_diff, int height_diff ) {
	unsigned char *res = optimal_image(width, height, width_diff, height_diff, output);
	save_image(output_file_name, width - width_diff, height - height_diff, res);
	free(res);
	//stbi_image_free(loaded); TODO
	return 0;
}

int run_timed(int width, int height, unsigned char *output, const char *output_file_name, int width_diff, int height_diff ) {
	unsigned char *res = optimal_image(width, height, width_diff, height_diff, output);
	free(res);
	return 0;
}

/* 
benchmarking function using rdtsc instruction
*/

double rdtsc(int width, int height, unsigned char *output, const char *output_file_name, int width_diff, int height_diff) {
    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;


#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i) {
            run_timed(width, height, output, output_file_name, width_diff, height_diff);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif
    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
		run_timed(width, height, output, output_file_name, width_diff, height_diff);	
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
	unsigned char *output;

	if (!load_image(argv[1], &width, &height, &output)) {
		return 1;
	}

	if(strcmp(argv[5],"1")){
		//print_matrix(output, width, height, 3);
		int out = run(width, height, output, argv[2], width_diff, height_diff);
		free(output);

		#ifdef count_instr
		printf("\nADDS=%llu MULTS=%llu\n", add_count, mult_count);
		#endif

		return out;
	}
	//time it
	double r = rdtsc(width, height, output, argv[2], width_diff, height_diff);

	free(output);

    printf("\nRDTSC instruction: %.0lf cycles measured\n", r);
	#ifdef count_instr
	printf("\nADDS=%llu MULTS=%llu\n", add_count, mult_count);
	#endif

    return 0;
}
