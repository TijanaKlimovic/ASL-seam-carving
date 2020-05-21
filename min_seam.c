#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "convolution.h"
#include "parse_img.h"

//--------------------	counter for instructions -------------------

#ifdef count_instr 
extern unsigned long long add_count;	//count the total number of add instructions
extern unsigned long long mult_count; 	//count the total number of mult instructions
#define COUNT(A, B) A += B;
#else
#define COUNT(A, B)
#endif

//------------------------------------------------------------------


// #define MIN2(X, Y, M, IDX) (X < Y) ? (M = X, IDX = 0) : (M = Y, IDX = 1);
// #define MIN21(X, Y, M, IDX) (X < Y) ? (M = X, IDX = -1) : (M = Y, IDX = 0);
// #define MIN3(X, Y, Z, M, IDX) (X < Y) && (X < Z) ? (M = X, IDX = -1) : MIN2(Y, Z, M, IDX);

#define MIN2(X, Y, M, IDX) if (X < Y) {M = X; IDX = 0;} else {M = Y; IDX = 1;}
#define MIN21(X, Y, M, IDX) if (X < Y) {M = X; IDX = -1;} else {M = Y; IDX = 0;}
#define MIN3(X, Y, Z, M, IDX) if ((X < Y) && (X < Z)) {M = X; IDX = -1;} else {MIN2(Y, Z, M, IDX)}

void rotate(int *energy, int *rotated, int h, int w) { 
    for (int i = 0; i < h; i++) { 
        for (int j = 0; j < w; j++) { 
            rotated[j * h + (h - i - 1)] = energy[i * w + j]; 
        } 
    } 
}

int min_seam(int rsize, int csize, unsigned char *img, int is_ver, int *ret_backtrack) {

	int *energy = (int *) malloc(rsize * csize * sizeof(int));
	int *padded_img = padd0_image(rsize, csize, img); //TODO try converting in pad to uchar
	calc_RGB_energy(rsize + 2, csize + 2, padded_img, energy);
	// contains index of the value from the prev row/column from where we came here
	int *backtrack = (int *) malloc(rsize * csize * sizeof(int)); //different from what we returnCOUNT(mult_count, 2)

	// if horizontal seam -> rotate +90 degrees the energy map
	int *dp;
	if (is_ver == 0) {
		dp = malloc(rsize * csize * sizeof(int));
		rotate(energy, dp, rsize, csize);
		int tmp = rsize;
		rsize = csize;
		csize = tmp;
		free(energy);
	} else {
		dp = energy;	
	}

	// find vertical min seam
	// traverse matrix in row major order
	int column_lim = csize - 1;
	int min_idx;
	int min_val;
	for (int i = 1; i < rsize; i++) { //start from second row	
		// first col, j == 0
		int row = i * csize;
		int where = row;
		int where_before = where - csize;
		
		int val1 = dp[where_before];
		int val2 = dp[where_before + 1];
		MIN2(val1, 
			 val2, 
			 min_val, 
			 min_idx);

		backtrack[where] = min_idx;
		dp[where] += min_val;

		for (int j = 1; j < column_lim; j++) {

			where = row + j;
			where_before = where - csize;

			val1 = dp[where_before - 1];
			val2 = dp[where_before];
			int val3 =  dp[where_before + 1];
			MIN3(val1, 
				 val2, 
				 val3, 
				 min_val, 
				 min_idx)

			backtrack[where] = j + min_idx;
			dp[where] += min_val;
		}
		// last row, j == csize - 1
		where = row + column_lim;
		where_before = where - csize;
		
		val1 = dp[where_before - 1];
		val2 =  dp[where_before];
		MIN21(val1,
			 val2,
			 min_val,
			 min_idx)

		backtrack[where] = column_lim + min_idx;
		dp[where] += min_val;
	}

	//process the data to return in appropriate format
	int ret = INT_MAX;
	int direction = -1; //used in turning 2D backtrack into 1D
	int row_lim = rsize - 1;


	//find the index of the minimum value of last row in the dp matrix
	int last_row = row_lim  * csize;
	for (int cnt = 0; cnt < csize; cnt++) {
		int current = last_row + cnt;
		int min = dp[current];
		if (min < ret) {
			ret = min;
			direction = cnt;
		}
	}

	//return the 1D backtrack (only the min seam)
	// direction -= last_start;

	if (is_ver == 1) {
		for (int i = row_lim; i >= 0; i--) {
			ret_backtrack[i] = direction;
			direction = backtrack[last_row + direction];
			last_row -= csize;
		}
	} else {
		//convert back indexes
		for (int i = row_lim; i >= 0; i--) {
			int d = column_lim - direction;
			ret_backtrack[i] = d;
			direction = backtrack[last_row + direction];
			last_row -= csize;
		}
	}
	free(dp);
	free(backtrack);
	free(padded_img);		
	return ret;
}
