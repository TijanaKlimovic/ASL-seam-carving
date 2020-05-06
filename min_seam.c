#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "convolution.h"

//--------------------	counter for instructions -------------------

#ifdef count_instr 
extern int add_count;	//count the total number of add instructions
extern int mult_count; 	//count the total number of mult instructions
#define COUNT(A, B) A += B;
#else
#define COUNT(A, B)
#endif

//------------------------------------------------------------------


#define MIN2(X, Y, M, IDX) if (X < Y) {M = X; IDX = 0;} else {M = Y; IDX = 1;} COUNT(count_ifs, 1)

#define MIN3(X, Y, Z, M, IDX) if ((X < Y) && (X < Z)) {M = X; IDX = -1;} else {MIN2(Y, Z, M, IDX)} COUNT(count_ifs, 3)

// #define LOG(X) X
#define LOG(X)


int min_seam(int rsize, int csize, int *img, int is_ver, int *ret_backtrack) {

	#ifdef count_instr        //counting adds and mults of this function
    int count_ifs = 0;        //includes explicit ifs and for loop ifs  -> ADDS
    int indexing = 0;         //includes increments of i.j,k variables  -> ADDS
    int pointer_adds = 0;     //pointer arithmetic                      -> ADDS
    int pointer_mults = 0;    //                                        -> MULTS
    //EACH COUNT IS FOR THE COMMAND(S) ABOVE IT
	#endif

	int *the_m = (int *) malloc(rsize * csize * sizeof(int));
	COUNT(mult_count, 2)

	//call Tijana's energy function to set the_m (M matrix starts as a copy of e1)
	int *padded_img = padd0_image(rsize, csize, img);

	calc_RGB_energy(rsize + 2, csize + 2, padded_img, the_m);
	COUNT(add_count, 2)

	int *backtrack = (int *) malloc(rsize * csize * sizeof(int)); //different from what we return
	COUNT(mult_count, 2)

	int out_cnt, in_cnt; //counters for loops
	int *row_ptr, *col_ptr; //for calculating where in the img we are (based on is_ver)
	int out_lim, in_lim; //limits for the counters of the loops
	int step, other_step; //step for calculating the neighbor pixel location (in the same row/col)

	if (is_ver == 1) {
		out_lim = rsize;
		in_lim = csize;
		row_ptr = &out_cnt;
		col_ptr = &in_cnt;
		step = 1;
		other_step = csize;
	} else {
		out_lim = csize;
		in_lim = rsize;
		row_ptr = &in_cnt;
		col_ptr = &out_cnt;
		step = csize;
		other_step = 1;
	}
	COUNT(count_ifs, 1)

	for (out_cnt = 1; out_cnt < out_lim; out_cnt++) { //start from second row/col
		
		for (in_cnt = 0; in_cnt < in_lim; in_cnt++) {

			//determine the location of the pixel based on is_ver
			int where = (*row_ptr) * csize + (*col_ptr);
			COUNT(pointer_adds, 1)
			COUNT(pointer_mults, 1)

			int where_before = where - other_step;
			COUNT(pointer_adds, 1)

			int min_idx;
			int min_val;
			//int min_energy;

			if (in_cnt == 0) {
				COUNT(count_ifs, 1) //if (in_cnt == 0)

				MIN2(the_m[where_before], 
					 the_m[where_before + step], 
					 min_val, 
					 min_idx)
				COUNT(pointer_adds, 3)

				backtrack[where] = min_idx;
				COUNT(pointer_adds, 1)

			} else if (in_cnt == in_lim - 1) {
				COUNT(count_ifs, 2) //if (in_cnt == 0) else if (in_cnt == in_lim - 1)
				COUNT(add_count, 1) //in_lim - 1

				MIN2(the_m[where_before - step],
					 the_m[where_before],
					 min_val,
					 min_idx)
				COUNT(pointer_adds, 3)

				min_idx--;
				COUNT(add_count, 1)

				backtrack[where] = in_cnt + min_idx;
				COUNT(pointer_adds, 1)
				COUNT(add_count, 1)
			} else {
				COUNT(count_ifs, 2) //if (in_cnt == 0) else if (in_cnt == in_lim - 1) else
				COUNT(add_count, 1) //in_lim - 1

				MIN3(the_m[where_before - step], 
					 the_m[where_before], 
					 the_m[where_before + step], 
					 min_val, 
					 min_idx)
				COUNT(pointer_adds, 5)

				backtrack[where] = in_cnt + min_idx;
				COUNT(pointer_adds, 1)
				COUNT(add_count, 1)
			}

			//min_energy = min_val;
			the_m[where] += min_val;
			COUNT(pointer_adds, 1)
			COUNT(add_count, 1)
		}

	}
	COUNT(count_ifs, out_lim + (out_lim - 1) * (in_lim + 1))
	COUNT(indexing, (out_lim - 1) + (out_lim - 1) * in_lim)

	LOG(
	printf ("[rsize] = %d, [csize] = %d\n", rsize, csize);
	
	for (int i = 0; i < rsize; i++) {
		
		for (int j = 0; j < csize; j++) {
			printf("%d ", backtrack[i * csize + j]);
		}

		printf("\n");
	}
	)
	//process the data to return in appropriate format
	int ret = INT_MAX;
	int direction = -1; //used in turning 2D backtrack into 1D

	//find the minimum of last row/col of the_m
	//set the counters to the beginning of the last row/col
	out_cnt--;
	COUNT(pointer_adds, 1)

	in_cnt = 0;

	int last_start = (*row_ptr) * csize + (*col_ptr);
	COUNT(pointer_adds, 1)
	COUNT(pointer_mults, 1)
	LOG(printf("last_start = %d\n", last_start));
	
	for (int cnt = 0; cnt < in_lim; cnt++) {
		int current = last_start + cnt * step;
		COUNT(pointer_adds, 1)
		COUNT(pointer_mults, 1)

		if (the_m[current] < ret) {

			ret = the_m[current];
			COUNT(pointer_adds, 1)

			direction = cnt;
		}
		COUNT(pointer_adds, 1)
		COUNT(count_ifs, 1)

	}
	COUNT(count_ifs, in_lim + 1)
	COUNT(indexing, in_lim)

	//return the 1D backtrack (only the min seam)
	// direction -= last_start;

	for (int i = out_lim - 1; i >= 0; i--) {
		LOG(printf("[%d] = %d\n", i, direction);)
		ret_backtrack[i] = direction;
		COUNT(pointer_adds, 1)

		LOG(printf("{%d %d %d}\n", backtrack[last_start + direction * step - step], backtrack[last_start + direction * step], backtrack[last_start + direction * step + step]);)
		direction = backtrack[last_start + direction * step];
		COUNT(pointer_adds, 2)
		COUNT(pointer_mults, 1)

		last_start -= other_step;
		COUNT(pointer_adds, 1)
	}
	COUNT(add_count, 1) //int i = out_lim - 1
	COUNT(count_ifs, out_lim)
	COUNT(indexing, out_lim - 1)

	free(the_m);
	free(backtrack);
	free(padded_img);

	COUNT(add_count, count_ifs)
	COUNT(add_count, indexing)
	COUNT(add_count, pointer_adds)
	COUNT(mult_count, pointer_mults)

	LOG(printf("DONE!\n");)
	return ret;
}

// void test() {
// 	int idx, idx2;
// 	double val, val2;
// 	MIN2(1, 10, val, idx)
// 	MIN3(-4, 9, 0, val2, idx2)
// 	printf("[%d] = %lf\n[%d] = %lf\n", idx, val, idx2, val2);
// 	double *mat = (double *) malloc(5 * 3 * sizeof(double));

// 	for (int i = 0; i < 5; i++) {

// 		for (int j = 0; j < 3; j++) {
// 			mat[i * 3 + j] = (i * 3 + j) / 2.0;
// 		}

// 	}

// 	for (int i = 0; i < 5; i++) {

// 		for (int j = 0; j < 3; j++) {
// 			printf("%lf ", mat[i * 3 + j]);
// 		}

// 		printf("\n");
// 	}

// 	// free(mat);
// }

// int main(int argc, char const *argv[]) {
// 	test();
// 	return 0;
// }