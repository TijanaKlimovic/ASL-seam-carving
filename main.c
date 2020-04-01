#include <stdio.h>
#include "parse_img.h"
#include "optimal_image.h"


extern int width, height;
extern unsigned char *original, *saved;
extern double *output;

int main(int argc, char const *argv[]) {
	int diff = 10;
	
	if (argc < 2) {
		printf("Usage: %s <image_path>\n", argv[0]);
	}

	if (!load_image(argv[1])) {
		return 1;
	}

	output = optimal_image(width, height, diff, diff, output);
	printf("finished optimal image\n");
	save_image("output.png", width - diff, height - diff, output, saved);
	free(output);
	free(saved);
	free(original);
	return 0;

}