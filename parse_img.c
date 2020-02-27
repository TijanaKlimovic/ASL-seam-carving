// compile with : gcc -Wall parse_img.c  -o parse_img -lm

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define C (3)

typedef unsigned char u8;

int width, height;
u8 *original, *r, *g, *b;

void copy_rgb(u8 *r, u8 *g, u8 *b, u8 *src) {
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < C; ++k) {
				// red
				if (k % C == 0)
					r[i * width + j] = src[i * width * C + j * C + k];
				// green
				else if (k % C == 1)
					g[i * width + j] = src[i * width * C + j * C + k];
				// blue
				else if (k % C == 2)
					b[i * width + j] = src[i * width * C + j * C + k];
			}
		}
	}
}

int allocate_buffer(int width, int height, u8 **buffer) {
	*buffer = malloc(width * height * sizeof(u8));
	if (*buffer == NULL) {
		printf("Failed to allocate buffer\n");
		return 0;
	}
	return 1;
}

int load_image(char *filename) {
	int n;
	original = stbi_load(filename, &width, &height, &n, C);
	if (original == NULL) {
		printf("Failed to load image %s\n", filename);
		return 0;
	}
	assert(n == C);
	printf("Loaded image (%d x %d)\n", width, height);
	if (!allocate_buffer(width, height, &r) 
	|| !allocate_buffer(width, height, &g)
	|| !allocate_buffer(width, height, &b)) 
		return 0;
	copy_rgb(r, g, b, original);
	return 1;
}

int main(int argc, char **argv) {
	if (argc < 2) {
		printf("Usage: %s <image_path>\n", argv[0]);
	}
	if (!load_image(argv[1])) {
		return 1;
	}
	return 0;
}