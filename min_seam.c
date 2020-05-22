#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "convolution.h"
#include "parse_img.h"
#include <immintrin.h>

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
#define MIN21(X, Y, M, IDX) if (X < Y) {M = X; IDX = -1;} else {M = Y; IDX = 0;}
#define MIN3(X, Y, Z, M, IDX) if ((X < Y) && (X < Z)) {M = X; IDX = -1;} else {MIN2(Y, Z, M, IDX)}

void rotate(int *energy, int *rotated, int h, int w) { 
    for (int i = 0; i < h; i++) { 
    	int rotated_idx = h-i-1;
    	int energy_idx = i*w;

    	int j = 0;
        for (; j < w-7; j+=8) {
	    	int cnt1 = j+1;
	    	int cnt2 = j+2;
	    	int cnt3 = j+3;
	    	int cnt4 = j+4;
	    	int cnt5 = j+5;
	    	int cnt6 = j+6;
	    	int cnt7 = j+7;

	        rotated[rotated_idx + j * h] = energy[energy_idx + j]; 
	        rotated[rotated_idx + cnt1 * h] = energy[energy_idx + cnt1]; 
	        rotated[rotated_idx + cnt2 * h] = energy[energy_idx + cnt2]; 
	        rotated[rotated_idx + cnt3 * h] = energy[energy_idx + cnt3]; 
	        rotated[rotated_idx + cnt4 * h] = energy[energy_idx + cnt4]; 
	        rotated[rotated_idx + cnt5 * h] = energy[energy_idx + cnt5]; 
	        rotated[rotated_idx + cnt6 * h] = energy[energy_idx + cnt6]; 
	        rotated[rotated_idx + cnt7 * h] = energy[energy_idx + cnt7];
        } 

        while (j < w) {
            rotated[rotated_idx + j * h] = energy[energy_idx + j]; 
            j++;
        }
    } 
}

int min_seam(int rsize, int csize, unsigned char *img, int is_ver, int *ret_backtrack) {
	int size = rsize * csize;
	int *energy = (int *) malloc(size * sizeof(int));
	short *padded_img = padd0_image(rsize, csize, img); //TODO try converting in pad to uchar
	calc_RGB_energy(rsize + 2, csize + 2, padded_img, energy);

	// contains index of the value from the prev row/column from where we came here
	int *backtrack = (int *) malloc(size * sizeof(int)); //different from what we returnCOUNT(mult_count, 2)

	// if horizontal seam -> rotate +90 degrees the energy map
	int *dp;
	if (is_ver == 0) {
		dp = malloc(size * sizeof(int));
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
	int prev_row_idx, prev_row_idx_aux;
	for (int i = 1; i < rsize; i++) { //start from second row	
		// first column, j == 0
		int row = i * csize;
		int where = row;
		int where_before = where - csize;
		
		MIN2(dp[where_before],
			 dp[where_before + 1],
			 min_val,
			 min_idx);
		backtrack[where] = min_idx;
		dp[where] += min_val;

		// SIMD
		__m256i idx1 = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
		__m256i idx2 = _mm256_set_epi32(8, 7, 6, 5, 4, 3, 2, 1);
		__m256i idx3 = _mm256_set_epi32(9, 8, 7, 6, 5, 4, 3, 2);
		__m256i incr = _mm256_set1_epi32(8);

		int j;
		prev_row_idx_aux = (i-1) * csize - 1;
		for (j = 1; j < column_lim-7; j+=8) {
			where = row + j;
			prev_row_idx = prev_row_idx_aux + j;
			// load
			__m256i y1 = _mm256_loadu_si256((__m256i *) (dp + prev_row_idx));
			__m256i y2 = _mm256_loadu_si256((__m256i *) (dp + prev_row_idx + 1));
			__m256i y3 = _mm256_loadu_si256((__m256i *) (dp + prev_row_idx + 2));
			__m256i t1 = _mm256_loadu_si256((__m256i *) (dp + where));

			// min(y1, y2)
			__m256i mask1 = _mm256_cmpgt_epi32(y1, y2);
			__m256i min1 = _mm256_blendv_epi8(y1, y2, mask1);
			__m256i min_idx1 = _mm256_blendv_epi8(idx1, idx2, mask1);

			// min(min(y1,y2), y3)
			__m256i mask2 = _mm256_cmpgt_epi32(min1, y3);
			__m256i min2 = _mm256_blendv_epi8(min1, y3, mask2);
			__m256i min_idx2 = _mm256_blendv_epi8(min_idx1, idx3, mask2);

			// add min to current energy value
			__m256i x1 = _mm256_add_epi32(t1, min2);
			// update indexes
			idx1 = _mm256_add_epi32(idx1, incr);
			idx2 = _mm256_add_epi32(idx2, incr);
			idx3 = _mm256_add_epi32(idx3, incr);
			// // store
			_mm256_storeu_si256((__m256i *)(dp + where), x1);
			_mm256_storeu_si256((__m256i *)(backtrack + where), min_idx2);		
		}

		//rest of the row
		for (; j < column_lim; ++j) {
			where = row + j;
			where_before = where - csize;

			MIN3(dp[where_before - 1], 
				 dp[where_before], 
				 dp[where_before + 1], 
				 min_val, 
				 min_idx)
			backtrack[where] = j + min_idx;
			dp[where] += min_val;
		}

		// last column, j == csize - 1
		where = row + column_lim;
		where_before = where - csize;
		
		MIN21(dp[where_before - 1],
			 dp[where_before],
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

    __m256i incr = _mm256_set1_epi32(8);
    __m256i idx = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    __m256i minindices = idx;
    __m256i minvalues = _mm256_loadu_si256((__m256i*)(dp+last_row));

	int cnt = 8;
	for (; cnt < csize-7; cnt+=8) {
		idx = _mm256_add_epi32(idx, incr);

		__m256i values = _mm256_loadu_si256((__m256i*)(dp + last_row + cnt));
        __m256i mask = _mm256_cmpgt_epi32(values, minvalues);
        minvalues = _mm256_blendv_epi8(minvalues, values, mask);
        minindices = _mm256_blendv_epi8(minindices, idx, mask);
	}

    // find min index in vector result (in an extremely naive way)
    int32_t values_array[8];
    uint32_t indices_array[8];

    _mm256_storeu_si256((__m256i*)values_array, minvalues);
    _mm256_storeu_si256((__m256i*)indices_array, minindices);

    ret = values_array[0];
    direction = indices_array[0];
    for (int i = 1; i < 8; i++) {
        if (values_array[i] < ret) {
            ret = values_array[i];
            direction = indices_array[i];
        }
    }

	while (cnt < csize) {
		int current = last_row + cnt;
		if (dp[current] < ret) {
			ret = dp[current];
			direction = cnt;
		}
		cnt++;
	}

	//return the 1D backtrack (only the min seam)
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

// used intrinsics:
//------------------
//__m256i _mm256_cmpgt_epi32 (__m256i a, __m256i b) - comparing ints
//__m256i _mm256_load_si256 (__m256i const * mem_addr) - loading ints
// __m256i _mm256_min_epi32 (__m256i a, __m256i b) - compare a, b and store min in dst
//void _mm256_store_si256 (__m256i * mem_addr, __m256i a)
// __m256i _mm256_blend_epi32 (__m256i a, __m256i b, const int imm8)
// __m256i _mm256_set_epi32 (int e7, int e6, int e5, int e4, int e3, int e2, int e1, int e0)
// __m256i _mm256_set1_epi32 (int a)
// int _mm256_movemask_epi8 (__m256i a)
// __m256i _mm256_andnot_si256 (__m256i a, __m256i b)
// __m256i _mm256_add_epi32 (__m256i a, __m256i b)
// void _mm256_store_si256 (__m256i * mem_addr, __m256i a)
// __m256i _mm256_blendv_epi8 (__m256i a, __m256i b, __m256i mask)