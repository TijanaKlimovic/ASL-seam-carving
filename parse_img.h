#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>	

void print_matrix(int *matrix, int width, int height, int channels);
int load_image(const char *filename, int *width, int *height, int **output);

void save_image(const char *filename, int width, int height, int *from);
void save_as_grayscale_image(char *filename, int new_width, int new_height, int *image);