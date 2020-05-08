#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "convolution.h"

//--------------------	counter for instructions -------------------

#ifdef count_instr 
extern unsigned long long add_count;	//count the total number of add instructions
extern unsigned long long mult_count; 	//count the total number of mult instructions
#endif

//------------------------------------------------------------------


#define MIN2(X, Y, M, IDX) if (X < Y) {M = X; IDX = 0;} else {M = Y; IDX = 1;}

#define MIN3(X, Y, Z, M, IDX) if ((X < Y) && (X < Z)) {M = X; IDX = -1;} else {MIN2(Y, Z, M, IDX)}

// #define LOG(X) X
#define LOG(X)

int min_seam(int rsize, int csize, int *img, int is_ver, int *ret_backtrack) {
	int *the_m = (int *) malloc(rsize * csize * sizeof(int));
	int *padded_img = padd0_image(rsize, csize, img);
	calc_RGB_energy(rsize + 2, csize + 2, padded_img, the_m);

	// contains index of the value from the prev row/column from where we came here
	int *backtrack = (int *) malloc(rsize * csize * sizeof(int)); //different from what we return

	// find vertical min seam
	if (is_ver) {
		for (int i = 1; i < rsize; i++) { //start from second row	
			for (int j = 0; j < csize; j++) {

				int where = i * csize + j;
				int where_before = where - csize;
				int min_idx;
				int min_val;

				// first col
				if (j == 0) {
					MIN2(the_m[where_before], 
						 the_m[where_before + 1], 
						 min_val, 
						 min_idx)

					backtrack[where] = min_idx;
				// last col
				} else if (j == csize - 1) {
					MIN2(the_m[where_before - 1],
						 the_m[where_before],
						 min_val,
						 min_idx)

					min_idx--;

					backtrack[where] = j + min_idx;
				} else {
					MIN3(the_m[where_before - 1], 
						 the_m[where_before], 
						 the_m[where_before + 1], 
						 min_val, 
						 min_idx)

					backtrack[where] = j + min_idx;
				}
				the_m[where] += min_val;
			}
		}

		//process the data to return in appropriate format
		int ret = INT_MAX;
		int direction = -1; //used in turning 2D backtrack into 1D

		//find the index of the minimum value of last row in the dp matrix
		int last_row = (rsize - 1)  * csize;
		for (int cnt = 0; cnt < csize; cnt++) {
			int current = last_row + cnt;
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}
		}

		//return the 1D backtrack (only the min seam)
		// direction -= last_start;

		for (int i = rsize - 1; i >= 0; i--) {
			ret_backtrack[i] = direction;
			direction = backtrack[last_row + direction];
			last_row -= csize;
		}

		free(the_m);
		free(backtrack);
		free(padded_img);
		
		return ret;

	} else {
	// find horizontal min seam
		for (int i = 1; i < csize; i++) { //start from second col
		
			for (int j = 0; j < rsize; j++) {

				int where = j * csize + i;
				int where_before = where - 1;
				int min_idx;
				int min_val;

				// first row
				if (j == 0) {
					MIN2(the_m[where_before], 
						 the_m[where_before + csize], 
						 min_val, 
						 min_idx)

					backtrack[where] = min_idx;

				// last row
				} else if (j == rsize - 1) {
					MIN2(the_m[where_before - csize],
						 the_m[where_before],
						 min_val,
						 min_idx)

					min_idx--;

					backtrack[where] = j + min_idx;
				} else {
					MIN3(the_m[where_before - csize], 
						 the_m[where_before], 
						 the_m[where_before + csize], 
						 min_val, 
						 min_idx)

					backtrack[where] = j + min_idx;
				}
				the_m[where] += min_val;
			}
		}
		//process the data to return in appropriate format
		int ret = INT_MAX;
		int direction = -1; //used in turning 2D backtrack into 1D

		//find the minimum of last row/col of the_m
		//set the counters to the beginning of the last row/col
		int last_col = csize - 1;
		
		for (int cnt = 0; cnt < rsize; cnt++) {
			int current = last_col + (cnt * csize);
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}
		}

		//return the 1D backtrack (only the min seam)
		// direction -= last_start;

		for (int i = csize - 1; i >= 0; i--) {
			ret_backtrack[i] = direction;
			direction = backtrack[last_col + (direction * csize)];
			last_col -= 1;
		}
		free(the_m);
		free(backtrack);
		free(padded_img);
		
		return ret;
	}
}
