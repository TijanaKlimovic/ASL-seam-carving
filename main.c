#include <stdio.h>
#include "parse_img.h"
#include "optimal_image.h"

int main(int argc, char const *argv[]) {
	int width, height;
	double *output;
	int width_diff = 10;
	int height_diff = 10;
	
	if (argc < 3) {
		printf("Usage: %s <image_path> <output_file_name>\n", argv[0]);
		return 1;
	}

	if (!load_image(argv[1], &width, &height, &output)) {
		return 1;
	}

	output = optimal_image(width, height, width_diff, height_diff, output);
	printf("Finished seam carving\n");
	save_image(argv[2], width - width_diff, height - height_diff, output);
	printf("Saved image as %s with size (%d, %d) \n", argv[2], width - width_diff, height - height_diff);
	free(output);
	return 0;

}