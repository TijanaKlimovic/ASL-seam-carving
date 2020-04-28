#include <stdlib.h>
#include <stdio.h>
#include "../parse_img.h"
#include "../optimal_image.h"
#include "../min_seam.h"
#include "../convolution.h"
#include <string.h>

int compare_matrix(double *m, double *expected, int h, int w, int channels) {
	for (int k = 0; k < channels; ++k) {
		for (int i = 0; i < h; ++i) {
			for (int j = 0; j < w; ++j) {
				int index = k * h * w + i * w + j;
				if (m[index] != expected[index]) {
					return 1;
				}
			}
		}
	}
	return 0;
}

void test_pad_image(double *img, double *expected, int h, int w) {
	double *out = padd0_image(h, w, img); // need to free out
	if (compare_matrix(out, expected, h, w, 3) == 1) {
		printf("test_pad_image FAILED\n");
	}
	free(out);
}

void test_calc_RGB_energy(double *img, double *expected, int h, int w) {
	double *out = malloc(h * w * sizeof(double));
	double *padded = padd0_image(h, w, img);
	calc_RGB_energy(h+2, w+2, padded, out);
	if (compare_matrix(out, expected, h, w, 1) == 1) {
		printf("test_calc_RGB_energy FAILED\n");
	}
	free(padded);
	free(out);
}

void test_min_seam(double *img, double expected, int h, int w, int is_vertical) {
	int *backtrack;
	if (is_vertical) {
		backtrack = (int *)malloc(h * sizeof(int));
	} else {
		backtrack = (int *)malloc(w * sizeof(int));
	}
	
	double min_cost = min_seam(h, w, img, is_vertical, backtrack);
	if (min_cost != expected) {
		printf("test_min_seam FAILED\n");
	}
	free(backtrack);
}

void test_optimal_image(double *img, double *expected, int h, int w, int h_diff, int w_diff) {
	double *out;
	out = optimal_image(w, h, w_diff, h_diff, img);
	if (compare_matrix(out, expected, h-h_diff, w-w_diff, 3) == 1) {
		printf("test_optimal_image FAILED\n");
	}
	free(out);
}

int main(int argc, char const *argv[]) {
	// ------------ TEST IMAGE 1 -------------
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

	{	// test vertical seam
		double expected = 6196;
		test_min_seam(img, expected, height, width, 1);
	}

	{ // test horizontal seam
		double expected = 5958;
		test_min_seam(img, expected, height, width, 0);
	}

	{ // test 2nd horizontal seam
		double img[] = {
							180,184,179,
							162,170,220,

							227,228,176,
							162,204,196,

							241,241,239,
							236,240,245,
											};
		double expected = 9954;
		test_min_seam(img, expected, height - 1, width, 0);
	}

	{ // test 2nd vertical seam
		double img[] = {
							180,179,
							179,255,
							162,220,

							227,176,
							185,215,
							162,196,

							241,239,
							238,240,
							236,245,
											};
		double expected = 10064;
		test_min_seam(img, expected, height, width - 1, 1);
	}

	{	// test removing a vertical seam
		double expected[] = {
							180,179,
							179,255,
							162,220,

							227,176,
							185,215,
							162,196,

							241,239,
							238,240,
							236,245,
											};
		test_optimal_image(img, expected, height, width, 0, 1);
	}

	{	// test removing a horizontal seam
		// need to load image again as optimal_imgae free-s it up
		if (!load_image("small_tests/input/test_3.png", &width, &height, &img)) {
			printf("Cannot load image");
			return 1;
		}
		double expected[] = {
							180,184,179,
							162,170,220,

							227,228,176,
							162,204,196,

							241,241,239,
							236,240,245,
											};
		test_optimal_image(img, expected, height, width, 1, 0);
	}

	{	// test removing a vertical and a horizontal seam
		// need to load image again as optimal_imgae free-s it up
		if (!load_image("small_tests/input/test_3.png", &width, &height, &img)) {
			printf("Cannot load image");
			return 1;
		}
		double expected[] = {
							180,179,
							162,220,

							227,176,
							162,196,

							241,239,
							236,245,
											};
		test_optimal_image(img, expected, height, width, 1, 1);
	}

	{	// test removing 2 horizontal seams
		// need to load image again as optimal_imgae free-s it up
		if (!load_image("small_tests/input/test_3.png", &width, &height, &img)) {
			printf("Cannot load image");
			return 1;
		}
		double expected[] = {
							162,170,179,

							162,204,176,

							236,240,239,
											};
		test_optimal_image(img, expected, height, width, 2, 0);
	}	

    return 0;
}
