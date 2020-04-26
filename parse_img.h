#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void print_matrix(double *matrix, int width, int height, int channels);
int load_image(const char *filename, int *width, int *height, double **output);

void save_image(const char *filename, int width, int height, double *from);
void save_as_grayscale_image(char *filename, int new_width, int new_height, double *image);