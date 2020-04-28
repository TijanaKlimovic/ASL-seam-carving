#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "convolution.h"

#define MIN2(X, Y, M, IDX) if (X < Y) {M = X; IDX = 0;} else {M = Y; IDX = 1;}

#define MIN3(X, Y, Z, M, IDX) if ((X < Y) && (X < Z)) {M = X; IDX = -1;} else {MIN2(Y, Z, M, IDX)}

// #define LOG(X) X
#define LOG(X)


double min_seam(int rsize, int csize, double *img, int is_ver, int *ret_backtrack) {
	double *the_m = (double *) malloc(rsize * csize * sizeof(double));
	//call Tijana's energy function to set the_m (M matrix starts as a copy of e1)
	double *padded_img = padd0_image(rsize, csize, img);
	calc_RGB_energy(rsize + 2, csize + 2, padded_img, the_m);
	int *backtrack = (int *) malloc(rsize * csize * sizeof(int)); //different from what we return
	int out_cnt, in_cnt; //counters for loops
	int *row_ptr, *col_ptr; //for calculating where in the img we are (based on is_ver)
	int out_lim, in_lim; //limits for the counters of the loops
	int step, other_step; //step for calculating the neighbor pixel location (in the same row/col)

	// for (int i = 0; i < rsize; i++) {

	// 	for (int j = 0; j < csize; j++) {
	// 		printf("%lf ", the_m[i * csize + j]);
	// 	}

	// 	printf("\n");
	// }
	// printf("End of Conv\n");

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

	for (out_cnt = 1; out_cnt < out_lim; out_cnt++) { //start from second row/col
		
		for (in_cnt = 0; in_cnt < in_lim; in_cnt++) {
			//determine the location of the pixel based on is_ver
			int where = (*row_ptr) * csize + (*col_ptr);
			int where_before = where - other_step;
			int min_idx;
			double min_val;
			//double min_energy;

			if (in_cnt == 0) {

				MIN2(the_m[where_before], 
					 the_m[where_before + step], 
					 min_val, 
					 min_idx)

				backtrack[where] = min_idx;
			} else if (in_cnt == in_lim - 1) {

				MIN2(the_m[where_before - step],
					 the_m[where_before],
					 min_val,
					 min_idx)

				min_idx--;
				backtrack[where] = in_cnt + min_idx;
			} else {

				MIN3(the_m[where_before - step], 
					 the_m[where_before], 
					 the_m[where_before + step], 
					 min_val, 
					 min_idx)

				backtrack[where] = in_cnt + min_idx;
			}

			//min_energy = min_val;
			the_m[where] += min_val;
		}

	}

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
	double ret = LONG_MAX;
	int direction = -1; //used in turning 2D backtrack into 1D

	//find the minimum of last row/col of the_m
	//set the counters to the beginning of the last row/col
	out_cnt--;
	in_cnt = 0;
	int last_start = (*row_ptr) * csize + (*col_ptr);
	LOG(printf("last_start = %d\n", last_start));
	
	for (int cnt = 0; cnt < in_lim; cnt++) {
		int current = last_start + cnt * step;

		if (the_m[current] < ret) {
			ret = the_m[current];
			direction = cnt;
		}

	}

	//return the 1D backtrack (only the min seam)
	// direction -= last_start;

	for (int i = out_lim - 1; i >= 0; i--) {
		LOG(printf("[%d] = %d\n", i, direction);)
		ret_backtrack[i] = direction;
		LOG(printf("{%d %d %d}\n", backtrack[last_start + direction * step - step], backtrack[last_start + direction * step], backtrack[last_start + direction * step + step]);)
		direction = backtrack[last_start + direction * step];
		last_start -= other_step;
	}

	free(the_m);
	free(backtrack);
	free(padded_img);
	LOG(printf("DONE!\n");)
	return ret;
}

void test() {
	int idx, idx2;
	double val, val2;
	MIN2(1, 10, val, idx)
	MIN3(-4, 9, 0, val2, idx2)
	printf("[%d] = %lf\n[%d] = %lf\n", idx, val, idx2, val2);
	double *mat = (double *) malloc(5 * 3 * sizeof(double));

	for (int i = 0; i < 5; i++) {

		for (int j = 0; j < 3; j++) {
			mat[i * 3 + j] = (i * 3 + j) / 2.0;
		}

	}

	for (int i = 0; i < 5; i++) {

		for (int j = 0; j < 3; j++) {
			printf("%lf ", mat[i * 3 + j]);
		}

		printf("\n");
	}

	// free(mat);
}

// int main(int argc, char const *argv[]) {
// 	test();
// 	return 0;
// }