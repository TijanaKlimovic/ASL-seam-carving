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

int min_seam(int rsize, int csize, unsigned char *img, int is_ver, int *ret_backtrack) {

    #ifdef count_instr        //counting adds and mults of this function
    unsigned long long pointer_adds = 0;     //pointer arithmetic                      -> ADDS
    unsigned long long pointer_mults = 0;    //                                        -> MULTS
    #endif

	int *the_m = (int *) malloc(rsize * csize * sizeof(int));
	int *padded_img = padd0_image(rsize, csize, img); //TODO try converting in pad to uchar
	calc_RGB_energy(rsize + 2, csize + 2, padded_img, the_m);
	// contains index of the value from the prev row/column from where we came here
	int *backtrack = (int *) malloc(rsize * csize * sizeof(int)); //different from what we returnCOUNT(mult_count, 2)

	#ifdef count_instr
	// malloc count
	mult_count += 4;
	// param count
	add_count += 2;
	#endif

	// find vertical min seam
	if (is_ver) {
		// traverse matrix in row major order
		for (int i = 1; i < rsize; i++) { //start from second row	
			// first col, j == 0
			int where = i * csize;
			int where_before = where - csize;
			int min_idx;
			int min_val;

			#ifdef count_instr
			add_count++;
			mult_count++;
			#endif

			MIN2(the_m[where_before], 
				 the_m[where_before + 1], 
				 min_val, 
				 min_idx)
			backtrack[where] = min_idx;
			the_m[where] += min_val;

			#ifdef count_instr
			pointer_adds += 5;
			add_count++;
			#endif

			for (int j = 1; j < csize - 1; j++) {
				where = i * csize + j;
				where_before = where - csize;

				#ifdef count_instr
				add_count += 2;
				mult_count++;
				#endif

				MIN3(the_m[where_before - 1], 
					 the_m[where_before], 
					 the_m[where_before + 1], 
					 min_val, 
					 min_idx)
				backtrack[where] = j + min_idx;
				the_m[where] += min_val;

				#ifdef count_instr
				pointer_adds += 7;
				add_count += 2;
				#endif
			}

			// last col, j == csize - 1
			where = i * csize + csize - 1;
			where_before = where - csize;

			#ifdef count_instr
			add_count += 3;
			mult_count++;
			#endif

			MIN2(the_m[where_before - 1],
				 the_m[where_before],
				 min_val,
				 min_idx)
			min_idx--;
			backtrack[where] = csize - 1 + min_idx;
			the_m[where] += min_val;

			#ifdef count_instr
			pointer_adds += 5;
			add_count += 4;
			#endif
		}

		//process the data to return in appropriate format
		int ret = INT_MAX;
		int direction = -1; //used in turning 2D backtrack into 1D

		//find the index of the minimum value of last row in the dp matrix
		int last_row = (rsize - 1)  * csize;

		#ifdef count_instr
		add_count++;
		mult_count++;
		#endif

		for (int cnt = 0; cnt < csize; cnt++) {
			int current = last_row + cnt;
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}

			#ifdef count_instr
			// pointer count
			pointer_adds++;

			// value count
			add_count++;
			#endif
		}

		for (int i = rsize - 1; i >= 0; i--) {
			ret_backtrack[i] = direction;
			direction = backtrack[last_row + direction];
			last_row -= csize;

			#ifdef count_instr
			// pointer count
			pointer_adds += 3;

			// value count
			add_count++;
			#endif
		}

		#ifdef count_instr
		add_count += pointer_adds;
		mult_count += pointer_mults;
		#endif

		free(the_m);
		free(backtrack);
		free(padded_img);		
		return ret;
	} else {
		// find horizontal min seam
		for (int i = 1; i < csize; i++) { //start from second col
			// forst row, j == 0
			int where = i;
			int where_before = where - 1;
			int min_idx;
			int min_val;

			#ifdef count_instr
			add_count++;
			#endif

			MIN2(the_m[where_before], 
				 the_m[where_before + csize], 
				 min_val, 
				 min_idx)
			backtrack[where] = min_idx;
			the_m[where] += min_val;

			#ifdef count_instr
			pointer_adds += 5;
			add_count++;
			#endif

			for (int j = 1; j < rsize - 1; j++) {
				where = j * csize + i;
				where_before = where - 1;

				#ifdef count_instr
				add_count += 2;
				mult_count++;
				#endif

				MIN3(the_m[where_before - csize], 
					 the_m[where_before], 
					 the_m[where_before + csize], 
					 min_val, 
					 min_idx)
				backtrack[where] = j + min_idx;
				the_m[where] += min_val;

				#ifdef count_instr
				pointer_adds += 7;
				add_count += 2;
				#endif
			}

			// last row, j = rsize - 1
			where = (rsize - 1) * csize + i;
			where_before = where - 1;

			#ifdef count_instr
			add_count += 3;
			mult_count++;
			#endif

			MIN2(the_m[where_before - csize],
				 the_m[where_before],
				 min_val,
				 min_idx)
			min_idx--;
			backtrack[where] = rsize - 1 + min_idx;
			the_m[where] += min_val;

			#ifdef count_instr
			pointer_adds += 5;
			add_count += 4;
			#endif
		}

		//process the data to return in appropriate format
		int ret = INT_MAX;
		int direction = -1; //used in turning 2D backtrack into 1D

		//find the minimum of last row/col of the_m
		//set the counters to the beginning of the last row/col
		int last_col = csize - 1;

		#ifdef count_instr
		add_count++;
		#endif
		
		for (int cnt = 0; cnt < rsize; cnt++) {
			int current = last_col + (cnt * csize);
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}

			#ifdef count_instr
			// pointer count
			pointer_adds++;

			// value count
			add_count++;
			mult_count++;
			#endif
		}

		for (int i = csize - 1; i >= 0; i--) {
			ret_backtrack[i] = direction;
			direction = backtrack[last_col + (direction * csize)];
			last_col -= 1;

			#ifdef count_instr
			// pointer count
			pointer_adds += 3;
			pointer_mults++;

			// value count
			add_count++;
			#endif
		}

		#ifdef count_instr
		add_count += pointer_adds;
		mult_count += pointer_mults;
		#endif

		free(the_m);
		free(backtrack);
		free(padded_img);		
		return ret;
	}
}
