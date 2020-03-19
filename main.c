#include <stdio.h>
#include "parse_img.h"
#include "optimal_image.h"


extern int width, height;
extern unsigned char *original, *saved;
extern double *output;

int main(int argc, char const *argv[]) {
	
	if (argc < 2) {
		printf("Usage: %s <image_path>\n", argv[0]);
	}

	if (!load_image(argv[1])) {
		return 1;
	}

	output = optimal_image(width, height, width - 50, height - 50, output);
	save_image("output.png", width - 50, height - 50);
	return 0;
}