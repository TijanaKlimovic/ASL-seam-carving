#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "convolution.h"

//--------------------	counter for instructions -------------------

#ifdef count_instr 
extern unsigned long long add_count;	//count the total number of add instructions
extern unsigned long long mult_count; 	//count the total number of mult instructions
#define COUNT(A, B) A += B;
#else
#define COUNT(A, B)
#endif

//------------------------------------------------------------------

#define MIN2(X, Y, M, IDX) if (X < Y) {M = X; IDX = 0;} else {M = Y; IDX = 1;}

#define MIN3(X, Y, Z, M, IDX) if ((X < Y) && (X < Z)) {M = X; IDX = -1;} else {MIN2(Y, Z, M, IDX)}

// #define LOG(X) X
#define LOG(X)

int min_seam(int rsize, int csize, unsigned char *img, int is_ver, int *ret_backtrack) {

	#ifdef count_instr        //counting adds and mults of this function
    int pointer_adds = 0;     //pointer arithmetic                      -> ADDS
    int pointer_mults = 0;    //                                        -> MULTS
    //EACH COUNT IS FOR THE COMMAND(S) ABOVE IT
	#endif

	int *the_m = (int *) malloc(rsize * csize * sizeof(int));
	COUNT(mult_count, 2)

	int *padded_img = padd0_image(rsize, csize, img); //TODO try converting in pad to uchar
	calc_RGB_energy(rsize + 2, csize + 2, padded_img, the_m);
	COUNT(add_count, 2)

	// contains index of the value from the prev row/column from where we came here
	int *backtrack = (int *) malloc(rsize * csize * sizeof(int)); //different from what we returnCOUNT(mult_count, 2)
	COUNT(mult_count, 2)

	// find vertical min seam
	if (is_ver) {
		// traverse matrix in row major order
		for (int i = 1; i < rsize; i++) { //start from second row	
			for (int j = 0; j < csize; j++) {

				int where = i * csize + j;
				COUNT(pointer_adds, 1)
				COUNT(pointer_mults, 1)

				int where_before = where - csize;
				COUNT(pointer_adds, 1)
				int min_idx;
				int min_val;

				// first col
				if (j == 0) {
					MIN2(the_m[where_before], 
						 the_m[where_before + 1], 
						 min_val, 
						 min_idx)
					COUNT(pointer_adds, 3)

					backtrack[where] = min_idx;
					COUNT(pointer_adds, 1)

				// last col
				} else if (j == csize - 1) {
					COUNT(add_count, 1) //in_lim - 1
					MIN2(the_m[where_before - 1],
						 the_m[where_before],
						 min_val,
						 min_idx)
					COUNT(pointer_adds, 3)

					min_idx--;
					COUNT(add_count, 1)

					backtrack[where] = j + min_idx;
					COUNT(pointer_adds, 1)
					COUNT(add_count, 1)

				} else {
					COUNT(add_count, 1) //in_lim - 1

					MIN3(the_m[where_before - 1], 
						 the_m[where_before], 
						 the_m[where_before + 1], 
						 min_val, 
						 min_idx)
					COUNT(pointer_adds, 5)

					backtrack[where] = j + min_idx;
					COUNT(pointer_adds, 1)
					COUNT(add_count, 1)
				}
				the_m[where] += min_val;
				COUNT(pointer_adds, 1)
				COUNT(add_count, 1)
			}
		}

		//process the data to return in appropriate format
		int ret = INT_MAX;
		int direction = -1; //used in turning 2D backtrack into 1D

		//find the index of the minimum value of last row in the dp matrix
		int last_row = (rsize - 1)  * csize;
		COUNT(add_count, 1)
		COUNT(mult_count, 1)
		for (int cnt = 0; cnt < csize; cnt++) {
			int current = last_row + cnt;
			COUNT(add_count, 1)
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}
			COUNT(pointer_adds, 1)
		}

		for (int i = rsize - 1; i >= 0; i--) {
			ret_backtrack[i] = direction;
			COUNT(pointer_adds, 1)
			direction = backtrack[last_row + direction];
			COUNT(pointer_adds, 2)
			last_row -= csize;
			COUNT(pointer_adds, 1)
		}
		COUNT(add_count, 1) //int i = out_lim - 1

		free(the_m);
		free(backtrack);
		free(padded_img);

		COUNT(add_count, pointer_adds)
		COUNT(mult_count, pointer_mults)
		
		return ret;

	} else {
	// find horizontal min seam
		for (int i = 1; i < csize; i++) { //start from second col
		
			for (int j = 0; j < rsize; j++) {

				int where = j * csize + i;
				COUNT(pointer_adds, 1)
				COUNT(pointer_mults, 1)

				int where_before = where - 1;
				COUNT(pointer_adds, 1)

				int min_idx;
				int min_val;

				// first row
				if (j == 0) {
					MIN2(the_m[where_before], 
						 the_m[where_before + csize], 
						 min_val, 
						 min_idx)
					COUNT(pointer_adds, 3)

					backtrack[where] = min_idx;
					COUNT(pointer_adds, 1)

				// last row
				} else if (j == rsize - 1) {
					COUNT(add_count, 1) //in_lim - 1

					MIN2(the_m[where_before - csize],
						 the_m[where_before],
						 min_val,
						 min_idx)
					COUNT(pointer_adds, 3)

					min_idx--;
					COUNT(add_count, 1)

					backtrack[where] = j + min_idx;
					COUNT(pointer_adds, 1)
					COUNT(add_count, 1)

				} else {
					COUNT(add_count, 1) //in_lim - 1

					MIN3(the_m[where_before - csize], 
						 the_m[where_before], 
						 the_m[where_before + csize], 
						 min_val, 
						 min_idx)
					COUNT(pointer_adds, 5)

					backtrack[where] = j + min_idx;
					COUNT(pointer_adds, 1)
					COUNT(add_count, 1)
				}
				the_m[where] += min_val;
				COUNT(pointer_adds, 1)
				COUNT(add_count, 1)
			}
		}

		//process the data to return in appropriate format
		int ret = INT_MAX;
		int direction = -1; //used in turning 2D backtrack into 1D

		//find the minimum of last row/col of the_m
		//set the counters to the beginning of the last row/col
		int last_col = csize - 1;
		COUNT(add_count, 1)
		
		for (int cnt = 0; cnt < rsize; cnt++) {
			int current = last_col + (cnt * csize);
			COUNT(add_count, 1)
			COUNT(mult_count, 1)
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}
			COUNT(pointer_adds, 1)
		}

		for (int i = csize - 1; i >= 0; i--) {
			ret_backtrack[i] = direction;
			COUNT(pointer_adds, 1)
			direction = backtrack[last_col + (direction * csize)];
			COUNT(pointer_adds, 2)
			COUNT(pointer_mults, 1)
			last_col -= 1;
			COUNT(pointer_adds, 1)
		}
		free(the_m);
		free(backtrack);
		free(padded_img);

		COUNT(add_count, pointer_adds)
		COUNT(mult_count, pointer_mults)
		
		return ret;
	}
}
