#include <stdlib.h>
#include <stdio.h>
#include "../parse_img.h"
#include "../optimal_image.h"
#include "../convolution.h"
#include <string.h>

int compare_matrix(double *m, double *expected, int h, int w) {
	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < w; ++j) {
			int index = i * w + j;
			if (m[index] != expected[index]) {
				return 1;
			}
		}
	}
	return 0;
}

void test_pad_image(double *img, double *expected, int h, int w) {
	double *out = padd0_image(h, w, img); // need to free out
	if (compare_matrix(out, expected, h, w) == 1) {
		printf("test_pad_image FAILED\n");
	}
	free(out);
}

void test_calc_RGB_energy(double *img, double *expected, int h, int w) {
	double *out = malloc(h * w * sizeof(double));
	double *padded = padd0_image(h, w, img);
	calc_RGB_energy(h+2, w+2, padded, out);
	if (compare_matrix(out, expected, h, w) == 1) {
		printf("test_calc_RGB_energy FAILED\n");
	}
	free(padded);
	free(out);
}


int main(int argc, char const *argv[]) {
	// TEST IMAGE 1
	double *img;
	int width, height;
	if (!load_image("small_tests/input/test_3.png", &width, &height, &img)) {
		printf("Cannot load image");
		return 1;
	}

	{
		double expected[] = {
							0,0,0,0,0,
							0,180,184,179,0,
							0,179,212,255,0,
							0,162,170,220,0,
							0,0,0,0,0,

							0,0,0,0,0,
							0,227,228,176,0,
							0,185,236,215,0,
							0,162,204,196,0,
							0,0,0,0,0,

							0,0,0,0,0,
							0,241,241,239,
							0,238,240,240,0,
							0,236,240,245,0,
							0,0,0,0,0
											};
		test_pad_image(img, expected, height, width);
	}

	{
		double expected[] = {
							3886, 2836,4102,
							2858,362,2738,
							3808,2998,4024

								};
		test_calc_RGB_energy(img, expected, height, width);
	}

    return 0;


}
