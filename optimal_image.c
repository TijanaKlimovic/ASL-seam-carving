#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <string.h>
#include "min_seam.h"


//--------------------	counter for instructions -------------------

#ifdef count_instr 
extern unsigned long long add_count;	//count the total number of add instructions
extern unsigned long long mult_count; 	//count the total number of mult instructions
#endif

//------------------------------------------------------------------


// Data structure to hold information in a T cell
struct cell_T {
	// optimal cost to get to this dimension
	int optimal_cost;
	// corresponding matrix for this dimension,
	// obtained based on optimal cost
	// size is 3 * (width - row) * (height - column)
	unsigned char *i;
};

// height -> horizontal seam -> row
// width -> vertical seam -> column
void calculate(int width, int height, int T_width, int T_height, int width_diff, struct cell_T *T) {
	int T_index = T_height * width_diff + T_width;

	#ifdef count_instr
	add_count++;
	mult_count++;
	#endif

	if (T_height == 0) {
		// first row -> vertical seam only
		int T_index_left = T_height * width_diff + T_width - 1;
		unsigned char *image_left = T[T_index_left].i;
		int image_width = width - T_width + 1;
		int image_height = height - T_height;

		#ifdef count_instr
		add_count += 5;
		mult_count++;
		#endif

		int *backtrack = (int *)malloc(image_height * sizeof(int));
		int optimal_cost = min_seam(image_height, image_width, image_left, 1, backtrack);

		T[T_index].optimal_cost = T[T_index_left].optimal_cost + optimal_cost;
		T[T_index].i = (unsigned char *)malloc(3 * image_height * (image_width - 1) * sizeof(unsigned char));

		#ifdef count_instr
		// malloc param only
		add_count++;
		mult_count += 4;

		// variable add
		add_count++;
		#endif

		int k, l;
		for (k = 0; k < image_height; ++k) { // construct each row at a time
			int crr_col = 0;
			for (l = 0; l < image_width; ++l)
				if (l != backtrack[k]) { // check for elem to remove
					// remove R 
					T[T_index].i[k*(image_width-1)*3 + crr_col*3] = T[T_index_left].i[k*image_width*3 + l*3];
					// remove G
					T[T_index].i[k*(image_width-1)*3 + crr_col*3 + 1] = T[T_index_left].i[k*image_width*3 + l*3 + 1];
					// remove B
					T[T_index].i[k*(image_width-1)*3 + crr_col*3 + 2] = T[T_index_left].i[k*image_width*3 + l*3 + 2];
					crr_col++;

					#ifdef count_instr
					// indices only
					add_count += 15;
					mult_count += 12;

					// variable add
					add_count++;
					#endif
				}
		}

		#ifdef count_instr
		// index count
		add_count += image_width * image_height;
		#endif

		free(backtrack);
	} else if (T_width == 0) {
		// first column -> horizontal seam only
		int T_index_up = (T_height - 1) * width_diff + T_width;
		unsigned char *image_up = T[T_index_up].i;
		int image_width = width - T_width;
		int image_height = height - T_height + 1;

		#ifdef count_instr
		add_count += 5;
		mult_count++;
		#endif

		int *backtrack = (int *)malloc(image_width * sizeof(int));
		int optimal_cost = min_seam(image_height, image_width, image_up, 0, backtrack);

		T[T_index].optimal_cost = T[T_index_up].optimal_cost + optimal_cost;
		T[T_index].i = (unsigned char *)malloc(3 * (image_height - 1) * image_width * sizeof(unsigned char));

		#ifdef count_instr
		// malloc param only
		add_count++;
		mult_count += 4;

		// variable add
		add_count++;
		#endif

		int k, l;
		for (k = 0; k < image_width; ++k) { // construct each column at a time
			int crr_row = 0;
			for (l = 0; l < image_height; ++l)
				if (l != backtrack[k]) { // check for elem to remove
					// remove R
					T[T_index].i[crr_row*image_width*3 + k*3] = T[T_index_up].i[l*image_width*3 + k*3];
					// remove G
					T[T_index].i[crr_row*image_width*3 + k*3 + 1] = T[T_index_up].i[l*image_width*3 + k*3 + 1];
					// remove B
					T[T_index].i[crr_row*image_width*3 + k*3 + 2] = T[T_index_up].i[l*image_width*3 + k*3 + 2];
					crr_row++;

					#ifdef count_instr
					// indices only
					add_count += 11;
					mult_count += 12;

					// variable add
					add_count++;
					#endif
				}
		}

		#ifdef count_instr
		// index count
		add_count += image_width * image_height;
		#endif

		free(backtrack);
	} else {
		int T_index_left = T_height * width_diff + T_width - 1;
		unsigned char *image_left = T[T_index_left].i;
		int image_left_width = width - T_width + 1;
		int image_left_height = height - T_height;

		#ifdef count_instr
		add_count += 5;
		mult_count++;
		#endif

		int *backtrack_left = (int *)malloc(image_left_height * sizeof(int));
		int optimal_cost_left = min_seam(image_left_height, image_left_width, image_left, 1, backtrack_left);

		#ifdef count_instr
		// malloc param only
		mult_count++;
		#endif

		int T_index_up = (T_height - 1) * width_diff + T_width;
		unsigned char *image_up = T[T_index_up].i;
		int image_up_width = width - T_width;
		int image_up_height = height - T_height + 1;

		#ifdef count_instr
		add_count += 5;
		mult_count++;
		#endif

		int *backtrack_up = (int *)malloc(image_up_width * sizeof(int));
		int optimal_cost_up = min_seam(image_up_height, image_up_width, image_up, 0, backtrack_up);

		#ifdef count_instr
		// malloc param only
		mult_count++;
		#endif

		// remove column
		if (T[T_index_left].optimal_cost + optimal_cost_left <=
			T[T_index_up].optimal_cost + optimal_cost_up) {
			T[T_index].optimal_cost = T[T_index_left].optimal_cost + optimal_cost_left;
			T[T_index].i = (unsigned char *)malloc(3 * image_left_height * (image_left_width - 1) * sizeof(unsigned char));
			
			#ifdef count_instr
			// malloc param only
			add_count++;
			mult_count += 3;

			// variable add
			add_count++;
			#endif

			int k, l;
			for (k = 0; k < image_left_height; ++k) { // construct each row at a time
				int crr_col = 0;
				for (l = 0; l < image_left_width; ++l)
					if (l != backtrack_left[k]) { // check for elem to remove
						// remove R
						T[T_index].i[k*(image_left_width-1)*3 + crr_col*3] = T[T_index_left].i[k*image_left_width*3 + l*3];
						// remove G
						T[T_index].i[k*(image_left_width-1)*3 + crr_col*3 + 1] = T[T_index_left].i[k*image_left_width*3 + l*3 + 1];
						// remove B
						T[T_index].i[k*(image_left_width-1)*3 + crr_col*3 + 2] = T[T_index_left].i[k*image_left_width*3 + l*3 + 2];
						crr_col++;
					}
			}

			#ifdef count_instr
			// index count
			add_count += image_left_height * image_left_width;
			#endif

		// remove row
		} else {
			T[T_index].optimal_cost = T[T_index_up].optimal_cost + optimal_cost_up;
			T[T_index].i = (unsigned char *)malloc(3 * (image_up_height - 1) * image_up_width * sizeof(unsigned char));
			
			#ifdef count_instr
			// malloc param only
			add_count++;
			mult_count += 3;

			// variable add
			add_count++;
			#endif

			int k, l;
			for (k = 0; k < image_up_width; ++k) { // construct each column at a time
				int crr_row = 0;
				for (l = 0; l < image_up_height; ++l)
					if (l != backtrack_up[k]) { // check for elem to remove
						// remove R
						T[T_index].i[crr_row*image_up_width*3 + k*3] = T[T_index_up].i[l*image_up_width*3 + k*3];
						// remove G
						T[T_index].i[crr_row*image_up_width*3 + k*3 + 1] = T[T_index_up].i[l*image_up_width*3 + k*3 + 1];
						// remove B
						T[T_index].i[crr_row*image_up_width*3 + k*3 + 2] = T[T_index_up].i[l*image_up_width*3 + k*3 + 2];
						crr_row++;
					}
			}

			#ifdef count_instr
			// index count
			add_count += image_up_width * image_up_height;
			#endif
		}

		#ifdef count_instr
		// if operations count
		add_count++;
		#endif

		free(backtrack_left);
		free(backtrack_up);
	}
}

unsigned char *optimal_image(int width, int height, int width_diff,
	int height_diff, unsigned char *image) {
	width_diff++;
	height_diff++;

	#ifdef count_instr
	add_count += 2;
	#endif

	struct cell_T *T = (struct cell_T *)malloc(width_diff
		* height_diff * sizeof(struct cell_T));
	T[0].optimal_cost = 0;
	T[0].i = image;

	#ifdef count_instr
	// malloc param only
	mult_count += 2;
	#endif

	if (height_diff >= 3) {
		int j, k;
		for (j = 1; j < width_diff; ++j)
			calculate(width, height, j, 0, width_diff, T);
		// fill out second row
		j = 1;
		for (k = 0; k < width_diff; ++k)
				calculate(width, height, k, j, width_diff, T);
		// fill out third row + free 1.st row without first elem
		j = 2;
		for (k = 0; k < width_diff; ++k)
				calculate(width, height, k, j, width_diff, T);
		for (int i = 1; i < width_diff; ++i) {
			free(T[i].i);
		}

		for (j = 3; j < height_diff; ++j){
			// free row 2 before
			for (int i = (j-2) * width_diff; i < (j-1) * width_diff; ++i) {
				free(T[i].i);
			}
			for (k = 0; k < width_diff; ++k)
				calculate(width, height, k, j, width_diff, T);
		}

		#ifdef count_instr
		// if operations count
		add_count += width_diff + width_diff * height_diff;
		#endif

		// copy 
		unsigned char *res = malloc(3 * (width - width_diff + 1) * (height - height_diff + 1) * sizeof(unsigned char));
		memcpy(res, T[width_diff * height_diff - 1].i, 3 * (width - width_diff + 1) * (height - height_diff + 1) * sizeof(unsigned char));

		// free last 2 rows
		for (int i = (height_diff - 2) * width_diff; i < width_diff * height_diff; ++i) {
			free(T[i].i);
		}

		#ifdef count_instr
		// malloc param only
		add_count += 4;
		mult_count += 3;

		// memcpy count all
		add_count += 5;
		mult_count += 4;

		// index count
		add_count += width_diff * height_diff - 1;
		#endif
		free(T);
		return res;

	} else {
		int j, k;
		for (j = 1; j < width_diff; ++j)
			calculate(width, height, j, 0, width_diff, T);

		for (j = 1; j < height_diff; ++j)
			for (k = 0; k < width_diff; ++k)
				calculate(width, height, k, j, width_diff, T);

		#ifdef count_instr
		// if operations count
		add_count += width_diff + width_diff * height_diff;
		#endif

		// copy 
		unsigned char *res = malloc(3 * (width - width_diff + 1) * (height - height_diff + 1) * sizeof(unsigned char));
		memcpy(res, T[width_diff * height_diff - 1].i, 3 * (width - width_diff + 1) * (height - height_diff + 1) * sizeof(unsigned char));

		// free last 2 rows
		for (int i = 1; i < width_diff * height_diff; ++i) {
			free(T[i].i);
		}

		#ifdef count_instr
		// malloc param only
		add_count += 4;
		mult_count += 3;

		// memcpy count all
		add_count += 5;
		mult_count += 4;

		// index count
		add_count += width_diff * height_diff - 1;
		#endif
		free(T);
		return res;
	}

	
}