#include <stdlib.h>
#include <stdio.h>
#include "../parse_img.h"
#include "../optimal_image.h"
#include "../min_seam.h"
#include "../convolution.h"
#include <string.h>

void write_out_matrix(double *img, int h, int w, int channels, char *fname) {
	FILE *f = fopen(fname, "w");
	if (f == NULL){
	    printf("Error opening file %s!\n", fname);
	    exit(1);
	}
	for (int k = 0; k < channels; ++k) {
		for (int i = 0; i < h; ++i) {
			for (int j = 0; j < w; ++j) {
				fprintf(f, "%d ", (int)img[k * h * w + i * w + j]);
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void write_out(double val, char *fname) {
	FILE *f = fopen(fname, "w");
	if (f == NULL){
	    printf("Error opening file %s!\n", fname);
	    exit(1);
	}
	fprintf(f, "%d\n", (int) val);
	fclose(f);
}

void test_pad_image(double *img, int h, int w, char *fname) {
	double *padded = padd0_image(h, w, img);
	write_out_matrix(padded, h + 2, w + 2, 3, fname);
	free(padded);
}

void test_calc_rgb_energy(double *img, int h, int w, char *fname) {
	double *out = malloc(h * w * sizeof(double));
	double *padded = padd0_image(h, w, img);
	calc_RGB_energy(h+2, w+2, padded, out);
	write_out_matrix(out, h, w, 1, fname);
	free(padded);
	free(out);
}

void test_min_seam(double *img, int h, int w, int is_vertical, char *fname) {
	int *backtrack;
	if (is_vertical) {
		backtrack = (int *)malloc(h * sizeof(int));
	} else {
		backtrack = (int *)malloc(w * sizeof(int));
	}
	
	double min_cost = min_seam(h, w, img, is_vertical, backtrack);
	free(backtrack);
	write_out(min_cost, fname);
}

void test_optimal_image(double *img, int h, int w, int h_diff, int w_diff, char *fname) {
	double *out;
	out = optimal_image(w, h, w_diff, h_diff, img);
	write_out_matrix(out, h - h_diff, w - w_diff, 3, fname);
	free(out);
}

int main(int argc, char const *argv[]) {
	if (argc < 2) {
		printf("Usage: %s <image_path>\n", argv[0]);
		return 1;
	}

	const char *path = argv[1];

	double *img;
	int width, height;
	if (!load_image(path, &width, &height, &img)) {
		printf("Cannot load image");
		return 1;
	}

	{
	// test pad image
	test_pad_image(img, height, width, "unit_tests/out_c/padding.txt");
	}

	{
	// test calc rgb energy
	test_calc_rgb_energy(img, height, width, "unit_tests/out_c/energy_map.txt");
	}
	
	{
	// test min seam - vertical
	test_min_seam(img, height, width, 1, "unit_tests/out_c/min_seam_cost_v.txt");
	}

	{
	// test min seam - horizontal
	test_min_seam(img, height, width, 0, "unit_tests/out_c/min_seam_cost_h.txt");
	}

	{
	// test optimal image - remove 1 vertical seam
	test_optimal_image(img, height, width, 0, 1, "unit_tests/out_c/seam_carved_1.txt");
	}

	{
	// test optimal image - remove 1 horizontal seam
	double *img;
	int width, height;
	if (!load_image(path, &width, &height, &img)) {
		printf("Cannot load image");
		return 1;
	}
	test_optimal_image(img, height, width, 1, 0, "unit_tests/out_c/seam_carved_2.txt");
	}

	{
	// test optimal image - remove 2 vertical seam
	double *img;
	int width, height;
	if (!load_image(path, &width, &height, &img)) {
		printf("Cannot load image");
		return 1;
	}
	test_optimal_image(img, height, width, 0, 2, "unit_tests/out_c/seam_carved_3.txt");
	}

	{
	double *img;
	int width, height;
	if (!load_image(path, &width, &height, &img)) {
		printf("Cannot load image");
		return 1;
	}
	// test optimal image - remove 2 horizontal seam
	test_optimal_image(img, height, width, 2, 0, "unit_tests/out_c/seam_carved_4.txt");
	}

    return 0;
}
