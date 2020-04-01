// compile with : gcc -Wall parse_img.c  -o parse_img -lm
// usage: ./parse_img <filename>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "parse_img.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define C (3)

int width, height;
unsigned char *original, *saved;
double *output;

void print_matrix(double *matrix, int width, int height) {
	for (int k = 0; k < 3; ++k) {
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				printf("%lf ", matrix[k * width * height + i * width + j]);
			}
			printf("\n");
		}
		printf("\n\n");	
	}
}

void convert_double(unsigned char *src, double *dst, int width, int height) {
	int i, j, k;
	for (i = 0; i < height; ++i) {
		for (j = 0; j < width; ++j) {
			for (k = 0; k < C; ++k) {
				dst[k * width * height + i * width + j] = (double)src[i * width * C + j * C + k];
			}
		}
	}
}

void convert_uchar(double *src, unsigned char *dst, int width, int height) {
	int i, j, k;
	for (i = 0; i < height; ++i) {
		for (j = 0; j < width; ++j) {
			for (k = 0; k < C; ++k) {
				dst[C * width * i + C * j + k] = (unsigned char)src[k * width * height + i * width + j];
			}
		}
	}
}

int allocate_double_buffer(int width, int height, double **buffer) {
	*buffer = malloc(width * height * C * sizeof(double));
	if (*buffer == NULL) {
		printf("Failed to allocate buffer\n");
		return 0;
	}
	return 1;
}

int allocate_uchar_buffer(int width, int height, unsigned char **buffer) {
	*buffer = malloc(width * height * C * sizeof(unsigned char));
	if (*buffer == NULL) {
		printf("Failed to allocate buffer\n");
		return 0;
	}
	return 1;
}

int load_image(const char *filename) {
	int n;
	original = stbi_load(filename, &width, &height, &n, C);
	if (original == NULL) {
		printf("Failed to load image %s\n", filename);
		return 0;
	}
	assert(n == C);
	printf("Loaded image (%d x %d)\n", width, height);
	if (!allocate_double_buffer(width, height, &output)) 
		return 0;
	printf("Allocated buffer\n");
	convert_double(original, output, width, height);
	printf("Converted to double\n");
	return 1;
}

void save_image(char *filename, int new_width, int new_height, double *buffer, unsigned char *output) {
	allocate_uchar_buffer(new_width, new_height, &output);
	convert_uchar(buffer, output, new_width, new_height);
	stbi_write_png(filename, new_width, new_height, C, output, new_width * C);
}

// int main(int argc, char **argv) {
// 	if (argc < 2) {
// 		printf("Usage: %s <image_path>\n", argv[0]);
// 	}
// 	if (!load_image(argv[1])) {
// 		return 1;
// 	}
// 	print_matrix(output, width, height);
// 	save_image("saved.png", width, height);
// 	return 0;
// }