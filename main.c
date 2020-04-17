#include "run.h"
#include "benchmarking.h"
#include <string.h>

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