#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>	

void print_matrix(unsigned char *matrix, int width, int height, int channels);
void print_matrix_int(int *matrix, int width, int height, int channels);
int load_image(const char *filename, int *width, int *height, unsigned char **output);

void save_image(const char *filename, int width, int height, unsigned char *from);
void save_as_grayscale_image(char *filename, int new_width, int new_height, int *image);