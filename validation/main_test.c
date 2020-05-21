#include <stdlib.h>
#include <stdio.h>
#include "../parse_img.h"
#include "../optimal_image.h"
#include "../min_seam.h"
#include "../convolution.h"
#include <string.h>

void print_matrix_short(short *matrix, int width, int height, int channels) {
	for (int k = 0; k < channels; ++k) {
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				printf("%12d ", matrix[k * width * height + i * width + j]);
			}
			printf("\n");
		}
		printf("\n\n");	
	}
}

int compare_matrix(unsigned char *m, unsigned char *expected, int h, int w, int channels) {
	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < w; ++j) {
			for (int k = 0; k < channels; k++) {
				int index = i * w * channels + j * channels + k;
				if (m[index] != expected[index]) {
					return 1;
				}
			}
		}
	}
	return 0;
}

int compare_matrix_int(int *m, int *expected, int h, int w, int channels) {
	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < w; ++j) {
			for (int k = 0; k < channels; k++) {
				int index = i * w * channels + j * channels + k;
				if (m[index] != expected[index]) {
					return 1;
				}
			}
		}
	}
	return 0;
}

int compare_matrix_short(short *m, short *expected, int h, int w, int channels) {
	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < w; ++j) {
			for (int k = 0; k < channels; k++) {
				int index = i * w * channels + j * channels + k;
				if (m[index] != expected[index]) {
					return 1;
				}
			}
		}
	}
	return 0;
}

void test_pad_image(unsigned char *img, short *expected, int h, int w) {
	short *out = padd0_image(h, w, img); // need to free out
	if (compare_matrix_short(out, expected, h+2, w+2, 3) == 1) {
		printf("test_pad_image FAILED\n");
		print_matrix_short(out, w+2, h+2, 3);
		print_matrix_short(expected, w+2, h+2, 3);
	} else {
		printf("test_pad_image PASSED\n");
	}
	
	free(out);
}

void test_calc_RGB_energy(unsigned char *img, int *expected, int h, int w) {
	int *out = malloc(h * w * sizeof(int));
	short *padded = padd0_image(h, w, img);
	//print_matrix(padded, w+2, h+2, 3);

	calc_RGB_energy(h+2, w+2, padded, out);
	if (compare_matrix_int(out, expected, h, w, 1) == 1) {
		printf("test_calc_RGB_energy FAILED\n");
	} else {
		printf("test_calc_RGB_energy PASSED\n");
	}
	//print_matrix(out, w, h, 1);
	free(padded);
	free(out);
}

void test_min_seam(unsigned char *img, int expected, int h, int w, int is_vertical) {
	int *backtrack;
	if (is_vertical) {
		backtrack = (int *)malloc(h * sizeof(int));
	} else {
		backtrack = (int *)malloc(w * sizeof(int));
	}
	
	int min_cost = min_seam(h, w, img, is_vertical, backtrack);
	if (min_cost != expected) {
		printf("test_min_seam FAILED\n");
	} else {
		printf("test_min_seam PASSED\n");
	}
	free(backtrack);
}

void test_optimal_image(unsigned char *img, unsigned char *expected, int h, int w, int h_diff, int w_diff) {
	unsigned char *out;
	out = optimal_image(w, h, w_diff, h_diff, img);
	if (compare_matrix(out, expected, h-h_diff, w-w_diff, 3) == 1) {
		printf("test_optimal_image FAILED\n");
		print_matrix(out, w, h, 3);
		print_matrix(expected, w, h, 3); 
	} else {
		printf("test_optimal_image PASSED\n");
	}
	free(out);
}

int main(int argc, char const *argv[]) {
	
	// ------------ TEST IMAGE 1 -------------
	{
	unsigned char *img;
	int width, height;
	if (!load_image("unit_tests/input_small/test_3.png", &width, &height, &img)) {
		printf("Cannot load image");
		return 1;
	}

	{
		short expected[] = {
							0,0,0, 0  ,0  ,0  , 0  ,0  ,0  , 0  ,  0,  0, 0,0,0,
							0,0,0, 180,227,241, 184,228,241, 179,176,239, 0,0,0,
							0,0,0, 179,185,238, 212,236,240, 255,215,240, 0,0,0,
							0,0,0, 162,162,236, 170,204,240, 220,196,245, 0,0,0,
							0,0,0, 0  ,0  ,0  , 0  ,0  ,0  , 0  ,  0,  0, 0,0,0
											};

		test_pad_image(img, expected, height, width);
	}

	{
		int expected[] = {
							3886, 2836,4102,
							2858,362,2738,
							3808,2998,4024
											};
		test_calc_RGB_energy(img, expected, height, width);
	}

	{	// test vertical seam
		int expected = 6196;
		test_min_seam(img, expected, height, width, 1);
	}

	{ // test horizontal seam
		int expected = 5958;
		test_min_seam(img, expected, height, width, 0);
	}

	{ // test 2nd horizontal seam
		unsigned char img[] = {
							180,227,241, 184,228,241, 179,176,239,
							162,162,236, 170,204,240, 220,196,245
											};
		int expected = 9954;
		test_min_seam(img, expected, height - 1, width, 0);
	}

	{ // test 2nd vertical seam
		unsigned char img[] = {
							180,227,241, 179,176,239,
							179,185,238, 255,215,240,
							162,162,236, 220,196,245
											};
		int expected = 10064;
		test_min_seam(img, expected, height, width - 1, 1);
	}

	{	// test edge case : remove 0 seam
		unsigned char expected[] = {
							180,227,241, 184,228,241, 179,176,239,
							179,185,238, 212,236,240, 255,215,240,
							162,162,236, 170,204,240, 220,196,245
											};
		test_optimal_image(img, expected, height, width, 0, 0);
	}

	{	// test removing a vertical seam
		// need to load image again as optimal_imgae free-s it up
		if (!load_image("unit_tests/input_small/test_3.png", &width, &height, &img)) {
			printf("Cannot load image");
			return 1;
		}
		unsigned char expected[] = {
							180,227,241, 179,176,239,
							179,185,238, 255,215,240,
							162,162,236, 220,196,245
											};
		test_optimal_image(img, expected, height, width, 0, 1);
	}

	{	// test removing a horizontal seam
		// need to load image again as optimal_imgae free-s it up
		if (!load_image("unit_tests/input_small/test_3.png", &width, &height, &img)) {
			printf("Cannot load image");
			return 1;
		}
		unsigned char expected[] = {
							180,227,241, 184,228,241, 179,176,239,
							162,162,236, 170,204,240, 220,196,245
											};
		test_optimal_image(img, expected, height, width, 1, 0);
	}

	{	// test removing a vertical and a horizontal seam
		// need to load image again as optimal_imgae free-s it up
		if (!load_image("unit_tests/input_small/test_3.png", &width, &height, &img)) {
			printf("Cannot load image");
			return 1;
		}
		unsigned char expected[] = {
							180,227,241, 179,176,239,
							162,162,236, 220,196,245
											};
		test_optimal_image(img, expected, height, width, 1, 1);
	}

	{	// test removing 2 horizontal seams
		// need to load image again as optimal_imgae free-s it up
		if (!load_image("unit_tests/input_small/test_3.png", &width, &height, &img)) {
			printf("Cannot load image");
			return 1;
		}
		unsigned char expected[] = {
							162,162,236, 170,204,240, 179,176,239, 
											};
		test_optimal_image(img, expected, height, width, 2, 0);
	}	

	{	// test removing 2 horizontal and vertical seams
		// need to load image again as optimal_imgae free-s it up
		if (!load_image("unit_tests/input_small/test_3.png", &width, &height, &img)) {
			printf("Cannot load image");
			return 1;
		}
		unsigned char expected[] = {
							162,162,236
											};
		test_optimal_image(img, expected, height, width, 2, 2);
	}
	}	

	// ------------ TEST IMAGE 2 -------------
	{
	unsigned char *img;
	int width, height;
	if (!load_image("unit_tests/input_small/test_2.png", &width, &height, &img)) {
		printf("Cannot load image");
		return 1;
	}

	{
		short expected[] = {
							0,0,0, 0  ,0 ,  0, 0  , 0, 0, 0  ,  0 ,0, 0,0,0,
							0,0,0, 191,53,105, 229,47,69, 251,105,64, 0,0,0,
							0,0,0, 162,58,137, 203,52,97, 238,48,60,  0,0,0,
							0,0,0, 0  ,0 ,  0, 0  , 0, 0, 0  ,  0 ,0, 0,0,0,

											};
		test_pad_image(img, expected, height, width);
	}

	{
		int expected[] = {
							2108, 1856, 2086,
							2092, 1898, 2234
											};
		test_calc_RGB_energy(img, expected, height, width);
	}

	{	// test vertical seam
		int expected = 3754;
		test_min_seam(img, expected, height, width, 1);
	}

	{ // test horizontal seam
		int expected = 6034;
		test_min_seam(img, expected, height, width, 0);
	}

	{	// test removing a vertical seam
		unsigned char expected[] = {
							191,53,105, 251,105,64,
							162,58,137, 238,48,60
											};
		test_optimal_image(img, expected, height, width, 0, 1);
	}

	{	// test removing a horizontal seam
		// need to load image again as optimal_imgae free-s it up
		if (!load_image("unit_tests/input_small/test_2.png", &width, &height, &img)) {
			printf("Cannot load image");
			return 1;
		}
		unsigned char expected[] = {
							191,53,105, 203,52,97, 238,48,60
											};
		test_optimal_image(img, expected, height, width, 1, 0);
	}

	{	// test removing a vertical and horizontal seam
		// need to load image again as optimal_imgae free-s it up
		if (!load_image("unit_tests/input_small/test_2.png", &width, &height, &img)) {
			printf("Cannot load image");
			return 1;
		}
		
		unsigned char expected[] = {
							191,53,105, 238,48,60
											};
		test_optimal_image(img, expected, height, width, 1, 1);
	}

	}
	
    return 0;
}
