#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "parse_img.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define C (3)

void print_matrix(unsigned char *matrix, int width, int height, int channels) {
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < channels; ++k) {
				printf("%12d ", matrix[i*width*channels + j*channels + k]);
			}
		}
		printf("\n");	
	}
	printf("\n");
}

void print_matrix_int(int *matrix, int width, int height, int channels) {
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

void convert_from_uchar_to_int(unsigned char *src, int *dst, int width, int height) {
	int i, j, k;
	for (i = 0; i < height; ++i) {
		for (j = 0; j < width; ++j) {
			for (k = 0; k < C; ++k) {
				dst[k * width * height + i * width + j] = (int)src[i * width * C + j * C + k];
			}
		}
	}
}

void convert_from_int_to_uchar(int *src, unsigned char *dst, int width, int height) {
	int i, j, k;
	for (i = 0; i < height; ++i) {
		for (j = 0; j < width; ++j) {
			for (k = 0; k < C; ++k) {
				dst[C * width * i + C * j + k] = (unsigned char)src[k * width * height + i * width + j];
			}
		}
	}
}

void convert_from_rgb(unsigned char *src, unsigned char *dst, int width, int height) {
	int i, j, k;
	for (i = 0; i < height; ++i) {
		for (j = 0; j < width; ++j) {
			for (k = 0; k < C; ++k) {
				dst[k * width * height + i * width + j] = src[i * width * C + j * C + k];
			}
		}
	}
}

void convert_to_rgb(unsigned char *src, unsigned char *dst, int width, int height) {
	int i, j, k;
	for (i = 0; i < height; ++i) {
		for (j = 0; j < width; ++j) {
			for (k = 0; k < C; ++k) {
				dst[C * width * i + C * j + k] = src[k * width * height + i * width + j];
			}
		}
	}
}

int allocate_int_buffer(int width, int height, int **buffer) {
	*buffer = malloc(width * height * C * sizeof(int));
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

int load_image(const char *filename, int *width, int *height, unsigned char **output) {
	int n;
	//unsigned char *loaded = stbi_load(filename, width, height, &n, C);
	*output = stbi_load(filename, width, height, &n, C);
	if (output == NULL) {
		printf("Failed to load image %s\n", filename);
		return 0;
	}
	assert(n == C);
	//if (!allocate_uchar_buffer(*width, *height, output)) 
	//	return 0;
	//convert_from_rgb(loaded, *output, *width, *height);
	//stbi_image_free(loaded);
	//print_matrix(*output, *width, *height, 3);
	return 1;
}

void save_image(const char *filename, int new_width, int new_height, unsigned char *buffer) {
	//unsigned char *output;
	//allocate_uchar_buffer(new_width, new_height, &output);
	//convert_to_rgb(buffer, output, new_width, new_height);
	stbi_write_png(filename, new_width, new_height, C, buffer, new_width * C);
	//free(output);
}

// converts the image pixels into the range of [0-255]
unsigned char* normalize_image(int* image, int height, int width) {
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

void save_as_grayscale_image(char *filename, int new_width, int new_height, int *image) {
	unsigned char *output = normalize_image(image, new_height, new_width);
	stbi_write_png(filename, new_width, new_height, 1, output, new_width);
	free(output);
}
