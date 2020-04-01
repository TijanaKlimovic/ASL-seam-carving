#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include "min_seam.h"

// Data structure to hold information in a T cell
struct cell_T {
	// optimal cost to get to this dimension
	double optimal_cost;
	// corresponding matrix for this dimension,
	// obtained based on optimal cost
	// size is 3 * (width - row) * (height - column)
	double *i;
};

// height -> horizontal seam -> row
// width -> vertical seam -> column
void calculate(int width, int height, int T_width, int T_height, 
	int wanted_width, int wanted_height, struct cell_T *T) {
	int T_index = T_height * wanted_width + T_width;

	if (T_height == 0) {
		// first row -> vertical seam only
		int T_index_left = T_height * wanted_width + T_width - 1;
		double *image_left = T[T_index_left].i;
		int image_width = width - T_width + 1;
		int image_height = height - T_height;

		int *backtrack = (int *)malloc(image_height * sizeof(int));
		double optimal_cost = min_seam(image_height, image_width, image_left, 1, backtrack);

		// printf("index %d backtrack left\n", T_index);
		// for (int k = 0; k < image_height; ++k)
		// 	printf("%d ", backtrack[k]);
		// printf("\n");

		T[T_index].optimal_cost = T[T_index_left].optimal_cost + optimal_cost;
		T[T_index].i = (double *)malloc(3 * image_height * (image_width - 1) * sizeof(double));
		int k, l;
		for (k = 0; k < image_height; ++k) { // construct each row at a time
			int crr_col = 0;
			for (l = 0; l < image_width; ++l)
				if (l != backtrack[k]) { // check for elem to remove
					T[T_index].i[k*(image_width-1)+crr_col]
						= T[T_index_left].i[k*image_width + l];
					T[T_index].i[image_height*(image_width-1)+k*(image_width-1)+crr_col]
						= T[T_index_left].i[image_height*(image_width-1)+k*image_width + l];
					T[T_index].i[2*image_height*(image_width-1)+k*(image_width-1)+crr_col]
						= T[T_index_left].i[2*image_height*(image_width-1)+k*image_width + l];
					crr_col++;
				}
		}
		free(backtrack);
	} else if (T_width == 0) {
		// first column -> horizontal seam only
		int T_index_up = (T_height - 1) * wanted_width + T_width;
		double *image_up = T[T_index_up].i;
		int image_width = width - T_width;
		int image_height = height - T_height + 1;

		int *backtrack = (int *)malloc(image_width * sizeof(int));
		double optimal_cost = min_seam(image_height, image_width, image_up, 0, backtrack);

		// printf("index %d backtrack up, size %d\n", T_index, image_width);
		// for (int k = 0; k < image_width; ++k)
		// 	printf("%d ", backtrack[k]);
		// printf("\n");

		T[T_index].optimal_cost = T[T_index_up].optimal_cost + optimal_cost;
		T[T_index].i = (double *)malloc(3 * (image_height - 1) * image_width * sizeof(double));
		int k, l;
		for (k = 0; k < image_width; ++k) { // construct each column at a time
			int crr_row = 0;
			for (l = 0; l < image_height; ++l)
				if (l != backtrack[k]) { // check for elem to remove
					T[T_index].i[crr_row*image_width+k]
						= T[T_index_up].i[l*image_width + k];
					T[T_index].i[(image_height-1)*image_width+crr_row*image_width+k]
						= T[T_index_up].i[(image_height-1)*image_width+l*image_width + k];
					T[T_index].i[2*(image_height-1)*image_width+crr_row*image_width+k]
						= T[T_index_up].i[2*(image_height-1)*image_width+l*image_width + k];
					crr_row++;
				}
		}
		free(backtrack);
	} else {
		int T_index_left = T_height * wanted_width + T_width - 1;
		double *image_left = T[T_index_left].i;
		int image_left_width = width - T_width + 1;
		int image_left_height = height - T_height;

		int *backtrack_left = (int *)malloc(image_left_height * sizeof(int));
		double optimal_cost_left = min_seam(image_left_height, image_left_width, image_left, 1, backtrack_left);

		// printf("index %d backtrack left\n", T_index);
		// for (int k = 0; k < image_left_height; ++k)
		// 	printf("%d ", backtrack_left[k]);
		// printf("\n");

		int T_index_up = (T_height - 1) * wanted_width + T_width;
		double *image_up = T[T_index_up].i;
		int image_up_width = width - T_width;
		int image_up_height = height - T_height + 1;

		int *backtrack_up = (int *)malloc(image_up_width * sizeof(int));
		double optimal_cost_up = min_seam(image_up_height, image_up_width, image_up, 1, backtrack_up);

		// printf("index %d backtrack up, size %d\n", T_index, image_up_width);
		// for (int k = 0; k < image_up_width; ++k)
		// 	printf("%d ", backtrack_up[k]);
		// printf("\n");

		if (T[T_index_left].optimal_cost + optimal_cost_left <=
			T[T_index_up].optimal_cost + optimal_cost_up) {
			T[T_index].optimal_cost = T[T_index_left].optimal_cost + optimal_cost_left;
			T[T_index].i = (double *)malloc(3 * image_left_height * (image_left_width - 1) * sizeof(double));
			int k, l;
			for (k = 0; k < image_left_height; ++k) { // construct each row at a time
				int crr_col = 0;
				for (l = 0; l < image_left_width; ++l)
					if (l != backtrack_left[k]) { // check for elem to remove
						T[T_index].i[k*(image_left_width-1)+crr_col]
							= T[T_index_left].i[k*image_left_width + l];
						T[T_index].i[image_left_height*(image_left_width-1)+k*(image_left_width-1)+crr_col]
							= T[T_index_left].i[image_left_height*(image_left_width-1)+k*image_left_width + l];
						T[T_index].i[2*image_left_height*(image_left_width-1)+k*(image_left_width-1)+crr_col]
							= T[T_index_left].i[2*image_left_height*(image_left_width-1)+k*image_left_width + l];
						crr_col++;
					}
			}
		} else {
			T[T_index].optimal_cost = T[T_index_up].optimal_cost + optimal_cost_up;
			T[T_index].i = (double *)malloc(3 * (image_up_height - 1) * image_up_width * sizeof(double));
			int k, l;
			for (k = 0; k < image_up_width; ++k) { // construct each column at a time
				int crr_row = 0;
				for (l = 0; l < image_up_height; ++l)
					if (l != backtrack_up[k]) { // check for elem to remove
						T[T_index].i[crr_row*image_up_width+k]
							= T[T_index_up].i[l*image_up_width + k];
						T[T_index].i[(image_up_height-1)*image_up_width+crr_row*image_up_width+k]
							= T[T_index_up].i[(image_up_height-1)*image_up_width+l*image_up_width + k];
						T[T_index].i[2*(image_up_height-1)*image_up_width+crr_row*image_up_width+k]
							= T[T_index_up].i[2*(image_up_height-1)*image_up_width+l*image_up_width + k];
						crr_row++;
					}
			}
		}
		free(backtrack_left);
		free(backtrack_up);
	}
}

double *optimal_image(int width, int height, int wanted_width,
	int wanted_height, double *image) {
	struct cell_T *T = (struct cell_T *)malloc(wanted_width
		* wanted_height * sizeof(struct cell_T));
	T[0].optimal_cost = 0;
	T[0].i = image;

	int j, k;
	for (j = 1; j < wanted_width; ++j)
		calculate(width, height, j, 0, wanted_width, wanted_height, T);
	for (j = 1; j < wanted_height; ++j)
		for (k = 0; k < wanted_width; ++k)
			calculate(width, height, k, j, wanted_width, wanted_height, T);
	for (int i = 0; i < wanted_width * wanted_height - 1; ++i) {
		free(T[i].i);
	}
	return T[wanted_width * wanted_height - 1].i;
}