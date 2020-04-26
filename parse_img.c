#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "parse_img.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define C (3)

void print_matrix(double *matrix, int width, int height, int channels) {
	for (int k = 0; k < channels; ++k) {
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				printf("%12lf ", matrix[k * width * height + i * width + j]);
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

void convert_grayscale_uchar(double *src, unsigned char *dst, int width, int height) {
	int i, j;
	for (i = 0; i < height; ++i) {
		for (j = 0; j < width; ++j) {
			dst[width * i + j] = (unsigned char)src[i * width + j];
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

int load_image(const char *filename, int *width, int *height, double **output) {
	int n;
	unsigned char *loaded;
	loaded = stbi_load(filename, width, height, &n, C);
	if (loaded == NULL) {
		printf("Failed to load image %s\n", filename);
		return 0;
	}
	assert(n == C);
	printf("Loaded image (%d x %d)\n", *width, *height);
	if (!allocate_double_buffer(*width, *height, output)) 
		return 0;
	convert_double(loaded, *output, *width, *height);
	return 1;
}

void save_image(const char *filename, int new_width, int new_height, double *buffer) {
	unsigned char *output;
	allocate_uchar_buffer(new_width, new_height, &output);
	convert_uchar(buffer, output, new_width, new_height);
	stbi_write_png(filename, new_width, new_height, C, output, new_width * C);
	free(output);
}

// converts the image pixels into the range of [0-255]
unsigned char* normalize_image(double* image, int height, int width) {
	unsigned char *normalized = malloc(height * width * sizeof(unsigned char));
	int max = 1;
    for (int i = 0; i < height; i++)
    	for (int j = 0; j < width; j++) {
    		if (image[i*width + j] > max) {
    			max = image[i*width + j];
    		} 
    }

    for (int i = 0; i < height; ++i) {
    	for (int j = 0; j < width; ++j) {
    		normalized[i*width + j] = (double) image[i*width + j] / max * 255;
    	}
    }
    return normalized;
}

void save_as_grayscale_image(char *filename, int new_width, int new_height, double *image) {
	unsigned char *output = normalize_image(image, new_height, new_width);
	stbi_write_png(filename, new_width, new_height, 1, output, new_width);
}


