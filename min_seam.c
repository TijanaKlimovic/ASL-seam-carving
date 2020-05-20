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
	#endif

    int size = rsize * csize;
	int *the_m = (int *) malloc(size * sizeof(int));

	int *padded_img = padd0_image(rsize, csize, img); //TODO try converting in pad to uchar
	calc_RGB_energy(rsize + 2, csize + 2, padded_img, the_m);
	COUNT(add_count, 2)

	// contains index of the value from the prev row/column from where we came here
	int *backtrack = (int *) malloc(size * sizeof(int)); //different from what we returnCOUNT(mult_count, 2)
	COUNT(mult_count, 3)

	int min_idx;
	int min_val;

	int last_rsize = rsize - 1;
	int last_csize = csize - 1;
	COUNT(add_count, 2)

	// find vertical min seam
	if (is_ver) {
		// traverse matrix in row major order
		for (int i = 1; i < rsize; i++) { //start from second row	
			int where_i = i * csize;
			int where_before = where_i - csize;
			COUNT(mult_count, 1)
			COUNT(add_count, 2)
			
			// j = 0
			MIN2(the_m[where_before],
				 the_m[where_before + 1], 
				 min_val, 
				 min_idx)
			backtrack[where_i] = min_idx;
			COUNT(pointer_adds, 4)

			the_m[where_i] += min_val;
			COUNT(pointer_adds, 1)
			COUNT(add_count, 1)

			for (int j = 1; j < last_csize; j++) {
				int where = where_i + j;
				where_before = where - csize;
				COUNT(add_count, 2)

				MIN3(the_m[where_before - 1], 
					 the_m[where_before], 
					 the_m[where_before + 1], 
					 min_val, 
					 min_idx)
				backtrack[where] = j + min_idx;
				COUNT(pointer_adds, 6)
				COUNT(add_count, 1)

				the_m[where] += min_val;
				COUNT(pointer_adds, 1)
				COUNT(add_count, 1)
			}

			// j = csize - 1
			int where = where_i + last_csize;
			where_before = where - csize;
			COUNT(add_count, 2)

			MIN2(the_m[where_before - 1],
				 the_m[where_before],
				 min_val,
				 min_idx)
			min_idx--;
			backtrack[where] = last_csize + min_idx;
			COUNT(pointer_adds, 4)
			COUNT(add_count, 2)

			the_m[where] += min_val;
			COUNT(pointer_adds, 1)
			COUNT(add_count, 1)
		}

		//process the data to return in appropriate format
		int ret = INT_MAX;
		int direction = -1; //used in turning 2D backtrack into 1D

		//find the index of the minimum value of last row in the dp matrix
		int last_row = last_rsize * csize;
		COUNT(mult_count, 1)
		int cnt = 0;
		for (; cnt < csize-11; cnt+=12) {
			int cnt1 = cnt + 1;
			int cnt2 = cnt + 2;
			int cnt3 = cnt + 3;
			int cnt4 = cnt + 4;
			int cnt5 = cnt + 5;
			int cnt6 = cnt + 6;
			int cnt7 = cnt + 7;
			int cnt8 = cnt + 8;
			int cnt9 = cnt + 9;
			int cnt10 = cnt + 10;
			int cnt11 = cnt + 11;

			int current0 = last_row + cnt;
			int current1 = last_row + cnt1;
			int current2 = last_row + cnt2;
			int current3 = last_row + cnt3;
			int current4 = last_row + cnt4;
			int current5 = last_row + cnt5;
			int current6 = last_row + cnt6;
			int current7 = last_row + cnt7;
			int current8 = last_row + cnt8;
			int current9 = last_row + cnt9;
			int current10 = last_row + cnt10;
			int current11 = last_row + cnt11;

			if (the_m[current0] < ret) {
				ret = the_m[current0];
				direction = cnt;
			}
			if (the_m[current1] < ret) {
				ret = the_m[current1];
				direction = cnt1;
			}
			if (the_m[current2] < ret) {
				ret = the_m[current2];
				direction = cnt2;
			}
			if (the_m[current3] < ret) {
				ret = the_m[current3];
				direction = cnt3;
			}
			if (the_m[current4] < ret) {
				ret = the_m[current4];
				direction = cnt4;
			}
			if (the_m[current5] < ret) {
				ret = the_m[current5];
				direction = cnt5;
			}
			if (the_m[current6] < ret) {
				ret = the_m[current6];
				direction = cnt6;
			}
			if (the_m[current7] < ret) {
				ret = the_m[current7];
				direction = cnt7;
			}
			if (the_m[current8] < ret) {
				ret = the_m[current8];
				direction = cnt8;
			}
			if (the_m[current9] < ret) {
				ret = the_m[current9];
				direction = cnt9;
			}
			if (the_m[current10] < ret) {
				ret = the_m[current10];
				direction = cnt10;
			}
			if (the_m[current11] < ret) {
				ret = the_m[current11];
				direction = cnt11;
			}
			
			COUNT(add_count, 19)
			COUNT(mult_count, 12)
			COUNT(pointer_adds, 12)
		}

		while (cnt < csize) {
			int current = last_row + cnt;
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}
			cnt++;

			COUNT(add_count, 2)
			COUNT(mult_count, 1)
			COUNT(pointer_adds, 1)
		}

		for (int i = last_rsize; i >= 0; i--) {
			ret_backtrack[i] = direction;
			COUNT(pointer_adds, 1)
			direction = backtrack[last_row + direction];
			COUNT(pointer_adds, 2)
			last_row -= csize;
			COUNT(pointer_adds, 1)
		}

		free(the_m);
		free(backtrack);
		free(padded_img);

		COUNT(add_count, pointer_adds)
		COUNT(mult_count, pointer_mults)
		
		return ret;
	} else {
		// find horizontal min seam
		for (int i = 1; i < csize; i++) { //start from second col
			// j = 0
			int where = i;
			int where_before = where - 1;
			COUNT(add_count, 1)

			MIN2(the_m[where_before], 
				 the_m[where_before + csize], 
				 min_val, 
				 min_idx)
			backtrack[where] = min_idx;
			COUNT(pointer_adds, 4)

			the_m[where] += min_val;
			COUNT(pointer_adds, 1)
			COUNT(add_count, 1)

			for (int j = 1; j < last_rsize; j++) {
				where = j * csize + i;
				where_before = where - 1;
				COUNT(add_count, 2)
				COUNT(mult_count, 1)

				MIN3(the_m[where_before - csize], 
					 the_m[where_before], 
					 the_m[where_before + csize], 
					 min_val, 
					 min_idx)
				backtrack[where] = j + min_idx;
				COUNT(pointer_adds, 6)
				COUNT(add_count, 1)

				the_m[where] += min_val;
				COUNT(pointer_adds, 1)
				COUNT(add_count, 1)
			}

			// j = rsize - 1
			where = last_rsize * csize + i;
			where_before = where - 1;
			COUNT(add_count, 2)
			COUNT(mult_count, 1)

			MIN2(the_m[where_before - csize],
				 the_m[where_before],
				 min_val,
				 min_idx)
			min_idx--;
			backtrack[where] = last_rsize + min_idx;
			COUNT(pointer_adds, 4)
			COUNT(add_count, 2)

			the_m[where] += min_val;
			COUNT(pointer_adds, 1)
			COUNT(add_count, 1)
		}

		//process the data to return in appropriate format
		int ret = INT_MAX;
		int direction = -1; //used in turning 2D backtrack into 1D

		//find the minimum of last row/col of the_m
		//set the counters to the beginning of the last row/col
		int cnt = 0;
		for (; cnt < rsize-11; cnt+=12) {
			int cnt1 = cnt + 1;
			int cnt2 = cnt + 2;
			int cnt3 = cnt + 3;
			int cnt4 = cnt + 4;
			int cnt5 = cnt + 5;
			int cnt6 = cnt + 6;
			int cnt7 = cnt + 7;
			int cnt8 = cnt + 8;
			int cnt9 = cnt + 9;
			int cnt10 = cnt + 10;
			int cnt11 = cnt + 11;

			int current0 = last_csize + cnt * csize;
			int current1 = last_csize + cnt1 * csize;
			int current2 = last_csize + cnt2 * csize;
			int current3 = last_csize + cnt3 * csize;
			int current4 = last_csize + cnt4 * csize;
			int current5 = last_csize + cnt5 * csize;
			int current6 = last_csize + cnt6 * csize;
			int current7 = last_csize + cnt7 * csize;
			int current8 = last_csize + cnt8 * csize;
			int current9 = last_csize + cnt9 * csize;
			int current10 = last_csize + cnt10 * csize;
			int current11 = last_csize + cnt11 * csize;

			if (the_m[current0] < ret) {
				ret = the_m[current0];
				direction = cnt;
			}
			if (the_m[current1] < ret) {
				ret = the_m[current1];
				direction = cnt1;
			} 
			if (the_m[current2] < ret) {
				ret = the_m[current2];
				direction = cnt2;
			} 
			if (the_m[current3] < ret) {
				ret = the_m[current3];
				direction = cnt3;
			}
			if (the_m[current4] < ret) {
				ret = the_m[current4];
				direction = cnt4;
			} 
			if (the_m[current5] < ret) {
				ret = the_m[current5];
				direction = cnt5;
			} 
			if (the_m[current6] < ret) {
				ret = the_m[current6];
				direction = cnt6;
			} 
			if (the_m[current7] < ret) {
				ret = the_m[current7];
				direction = cnt7;
			} 
			if (the_m[current8] < ret) {
				ret = the_m[current8];
				direction = cnt8;
			} 
			if (the_m[current9] < ret) {
				ret = the_m[current9];
				direction = cnt9;
			} 
			if (the_m[current10] < ret) {
				ret = the_m[current10];
				direction = cnt10;
			} 
			if (the_m[current11] < ret) {
				ret = the_m[current11];
				direction = cnt11;
			} 

			COUNT(add_count, 19)
			COUNT(mult_count, 12)
			COUNT(pointer_adds, 12)
		}

		while (cnt < rsize) {
			int current = last_csize + cnt * csize;
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}
			cnt++;

			COUNT(add_count, 2)
			COUNT(mult_count, 1)
			COUNT(pointer_adds, 1)
		}

		for (int i = last_csize; i >= 0; i--) {
			ret_backtrack[i] = direction;
			COUNT(pointer_adds, 1)
			direction = backtrack[last_csize + (direction * csize)];
			COUNT(pointer_adds, 2)
			COUNT(pointer_mults, 1)
			last_csize -= 1;
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
