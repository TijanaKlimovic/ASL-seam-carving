//#error Please comment out the next two lines under linux, then comment this error
//#include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler
#include <windows.h> // Include if under windows

#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "tsc_x86.h"
#include "run.h"

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


extern int width, height;
extern unsigned char *original, *saved;
extern double *output;

//--------------- extern variables --------------


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
*/

int run(const char* path) {
	int width_diff = 40;
	int height_diff = 20;


	if (!load_image(path)) {
		return 1;
	}

	output = optimal_image(width, height, width_diff, height_diff, output);
	printf("finished optimal image\n");
	save_image("output.png", width - width_diff, height - height_diff, output, saved);
	free(output);
	free(saved);
	free(original);
	return 0;
}


/* 
benchmarking function using rdtsc instruction
*/

double rdtsc(const char* path) {
    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;


#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i) {
            run(path);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        run(path);
    }

    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}


/* main function , has two options, first is the image path second is -t/<anything except nothing> */ 

int main(int argc, char const *argv[]) {
	
if (argc < 3) {
		printf("Usage: %s <image_path>\n", argv[0]);
		return 1;
	}
	// normal run
	if(strcmp(argv[2],"-t")){
		return run(argv[1]);
	}
	//time it
	double r = rdtsc(argv[1]);
    printf("RDTSC instruction:\n %lf cycles measured", r);
    return 0;
}


