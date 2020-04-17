#include "run.h"


extern int width, height;
extern unsigned char *original, *saved;
extern double *output;

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