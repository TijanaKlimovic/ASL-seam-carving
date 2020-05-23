#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "parse_img.h"
#include <immintrin.h>

//------------------------------------------------------------------
// for convolution.c
typedef union acc {
    __m256i intrin;
    short data[16];
} acc;

#define K 1
#define ABS(X) (((X)<0) ? (-(X)) : (X))

#define MIN2(X, Y, M, IDX) if (X < Y) {M = X; IDX = 0;} else {M = Y; IDX = 1;}
#define MIN21(X, Y, M, IDX) if (X < Y) {M = X; IDX = -1;} else {M = Y; IDX = 0;}
#define MIN3(X, Y, Z, M, IDX) if ((X < Y) && (X < Z)) {M = X; IDX = -1;} else {MIN2(Y, Z, M, IDX)}

int min_seam(int rsize, int csize, unsigned char *img, int is_ver, int *ret_backtrack) {

	int size = rsize * csize;
	int *energy = (int *) malloc(size * sizeof(int));
	//short *padded_img = padd0_image(rsize, csize, img); //TODO try converting in pad to uchar

	int padded_size = (rsize+2) * (csize+2) * 3;
	short* padded_img = (short*) malloc(padded_size*sizeof(short));

	//int padded_image[n+2][m+2][3];
	for(int i = 0 ; i < rsize+2 ; i++){
		int i1_idx = i*(csize+2)*3;
		int i2_idx = (i-1)*csize*3;
		for(int j = 0 ; j < csize+2 ; j++){
			if(i == 0 || j == 0 || i == rsize+1 || j == csize+1){
				//padded_image[i*(n+2)*(m+2) + (m+2)*j + k] = 0;
				padded_img[i1_idx + j*3] = 0;
				padded_img[i1_idx + j*3 + 1] = 0;
				padded_img[i1_idx + j*3 + 2] = 0;
			} else{
				padded_img[i1_idx + j*3] = (short) (img[i2_idx + (j-1)*3]);
				padded_img[i1_idx + j*3 + 1] = (short) (img[i2_idx + (j-1)*3 + 1]);
				padded_img[i1_idx + j*3 + 2] = (short) (img[i2_idx + (j-1)*3 + 2]);
			}
		}
	}	

	//calc_RGB_energy(rsize + 2, csize + 2, padded_img, energy);
{
	int n = rsize + 2;
	int m = csize + 2;
	short *padded = padded_img;

	//void calc_RGB_energy(int n, int m, short* padded, int* energy){
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    int i_limit = n - K;

    int block_width_L1 = 1141;      //working set size is 2*3*4(m+2) + m*4 < C
    int width_limit_L1 = m - K - block_width_L1 + 1;

    int jj, jj_old;
    int block_width_L1_9 = block_width_L1 - 9;
    for(jj = 1; jj < width_limit_L1; jj+= block_width_L1){
        int j_limit = jj + block_width_L1_9;

        for(int i = 1 ; i < i_limit ; i++){
            int j;
            int i1m = (i - 1) * m * 3;
            int i2m = i * m * 3;
            int i3m = (i + 1) * m * 3;
            int im = (i - 1) * (m - 2);
            for(j = jj ; j < j_limit ; j += 10){
                short *row0 = padded + i1m + (j - 1) * 3;
                short *row1 = padded + i2m * 3 + (j - 1) * 3;
                short *row2 = padded + i3m + (j - 1) * 3;
                int imj = im + j;

                //loads are all unaligned because we're dealing with shorts
                __m256i r1 = _mm256_loadu_si256((__m256i *) (row0 + 3));
                __m256i r6 = _mm256_loadu_si256((__m256i *) (row2 + 3));
                r1 = _mm256_sub_epi16(r6, r1);
                r1 = _mm256_add_epi16(r1, r1);

                __m256i r3 = _mm256_loadu_si256((__m256i *) row1);
                __m256i r4 = _mm256_loadu_si256((__m256i *) (row1 + 6));
                r3 = _mm256_sub_epi16(r4, r3);
                r3 = _mm256_add_epi16(r3, r3);

                __m256i r0 = _mm256_loadu_si256((__m256i *) row0);
                acc r2;
                r2.intrin = _mm256_loadu_si256((__m256i *) (row0 + 6));
                __m256i r5 = _mm256_loadu_si256((__m256i *) row2);
                r4 = _mm256_sub_epi16(r2.intrin, r0);
                r6 = _mm256_sub_epi16(r5, r0);

                acc r7;
                r7.intrin = _mm256_loadu_si256((__m256i *) (row2 + 6));
                r0 = _mm256_sub_epi16(r7.intrin, r5);
                r5 = _mm256_sub_epi16(r7.intrin, r2.intrin);

                r2.intrin = _mm256_add_epi16(r4, r3);
                r2.intrin = _mm256_add_epi16(r2.intrin, r0);
                r7.intrin = _mm256_add_epi16(r1, r6);
                r7.intrin = _mm256_add_epi16(r7.intrin, r5);

                __m256i r9 = _mm256_loadu_si256((__m256i *) (row0 + 18));
                __m256i r14 = _mm256_loadu_si256((__m256i *) (row2 + 18));
                r9 = _mm256_sub_epi16(r14, r9);
                r9 = _mm256_add_epi16(r9, r9);

                __m256i r11 = _mm256_loadu_si256((__m256i *) (row1 + 15));
                __m256i r12 = _mm256_loadu_si256((__m256i *) (row1 + 21));
                r11 = _mm256_sub_epi16(r12, r11);
                r11 = _mm256_add_epi16(r11, r11);

                __m256i r8 = _mm256_loadu_si256((__m256i *) (row0 + 15));
                acc r10;
                r10.intrin = _mm256_loadu_si256((__m256i *) (row0 + 21));
                __m256i r13 = _mm256_loadu_si256((__m256i *) (row2 + 15));
                r12 = _mm256_sub_epi16(r10.intrin, r8);
                r14 = _mm256_sub_epi16(r13, r8);

                acc r15;
                r15.intrin = _mm256_loadu_si256((__m256i *) (row2 + 21));
                r8 = _mm256_sub_epi16(r15.intrin, r13);
                r13 = _mm256_sub_epi16(r15.intrin, r10.intrin);

                r10.intrin = _mm256_add_epi16(r12, r11);
                r10.intrin = _mm256_add_epi16(r10.intrin, r8);
                r15.intrin = _mm256_add_epi16(r9, r14);
                r15.intrin = _mm256_add_epi16(r15.intrin, r13);

                r2.intrin = _mm256_abs_epi16(r2.intrin);
                r7.intrin = _mm256_abs_epi16(r7.intrin);
                r10.intrin = _mm256_abs_epi16(r10.intrin);
                r15.intrin = _mm256_abs_epi16(r15.intrin);

                short result0 =  r2.data[ 0] +  r2.data[ 1] +  r2.data[ 2] +  r7.data[ 0] +  r7.data[ 1] +  r7.data[ 2];
                short result1 =  r2.data[ 3] +  r2.data[ 4] +  r2.data[ 5] +  r7.data[ 3] +  r7.data[ 4] +  r7.data[ 5];
                short result2 =  r2.data[ 6] +  r2.data[ 7] +  r2.data[ 8] +  r7.data[ 6] +  r7.data[ 7] +  r7.data[ 8];
                short result3 =  r2.data[ 9] +  r2.data[10] +  r2.data[11] +  r7.data[ 9] +  r7.data[10] +  r7.data[11];
                short result4 =  r2.data[12] +  r2.data[13] +  r2.data[14] +  r7.data[12] +  r7.data[13] +  r7.data[14];

                short result5 = r10.data[ 0] + r10.data[ 1] + r10.data[ 2] + r15.data[ 0] + r15.data[ 1] + r15.data[ 2];
                short result6 = r10.data[ 3] + r10.data[ 4] + r10.data[ 5] + r15.data[ 3] + r15.data[ 4] + r15.data[ 5];
                short result7 = r10.data[ 6] + r10.data[ 7] + r10.data[ 8] + r15.data[ 6] + r15.data[ 7] + r15.data[ 8];
                short result8 = r10.data[ 9] + r10.data[10] + r10.data[11] + r15.data[ 9] + r15.data[10] + r15.data[11];
                short result9 = r10.data[12] + r10.data[13] + r10.data[14] + r15.data[12] + r15.data[13] + r15.data[14];

                *(energy + imj - 1) = (int) result0;
                *(energy + imj    ) = (int) result1;
                *(energy + imj + 1) = (int) result2;
                *(energy + imj + 2) = (int) result3;
                *(energy + imj + 3) = (int) result4;
                *(energy + imj + 4) = (int) result5;
                *(energy + imj + 5) = (int) result6;
                *(energy + imj + 6) = (int) result7;
                *(energy + imj + 7) = (int) result8;
                *(energy + imj + 8) = (int) result9;
            }
            
            for(; j < jj + block_width_L1 ; j++) {
                int acc1;
                int acc2;
                int acc3;
                int acc4;
                int acc5;
                int acc6;

                int j1 = (j - 1) * 3;
                int j2 = j * 3;
                int j3 = (j + 1) * 3;
                int *energy_pos = energy + im + (j-1);

                // channel R
                //H_y
                acc1 = -(padded[i1m  + j1] + ((padded[i1m + j2]) << 1));
                acc2 = padded[i3m + j1] - padded[i1m + j3];
                acc3 = ((padded[i3m + j2]) << 1) + padded[i3m + j3];
                *(energy_pos) = (int) ABS(acc1 + acc2 + acc3);
                //H_x
                acc4 = padded[i1m + j3] - padded[i1m + j1];
                acc5 = (padded[i2m + j3] - padded[i2m + j1]) << 1;
                acc6 = padded[i3m + j3] - padded[i3m + j1];
                *(energy_pos) += (int) ABS(acc4 + acc5 + acc6);

                // channel G
                //H_y
                int k = 1;
                acc1 = -(padded[i1m  + j1 + k] + ((padded[i1m + j2 + k]) << 1));
                acc2 = padded[i3m + j1 + k] - padded[i1m + j3 + k];
                acc3 = ((padded[i3m + j2 + k]) << 1) + padded[i3m + j3 + k];
                *(energy_pos) += (int) ABS(acc1 + acc2 + acc3);
                //H_x
                acc4 = padded[i1m + j3 + k] - padded[i1m + j1 + k];
                acc5 = (padded[i2m + j3 + k] - padded[i2m + j1 + k]) << 1;
                acc6 = padded[i3m + j3 + k] - padded[i3m + j1 + k];
                *(energy_pos) += (int) ABS(acc4 + acc5 + acc6);

                // channel B
                //H_y
                k = 2;
                acc1 = -(padded[i1m  + j1 + k] + ((padded[i1m + j2 + k]) << 1));
                acc2 = padded[i3m + j1 + k] - padded[i1m + j3 + k];
                acc3 = ((padded[i3m + j2 + k]) << 1) + padded[i3m + j3 + k];
                *(energy_pos) += (int) ABS(acc1 + acc2 + acc3);
                //H_x
                acc4 = padded[i1m + j3 + k] - padded[i1m + j1 + k];
                acc5 = (padded[i2m + j3 + k] - padded[i2m + j1 + k]) << 1;
                acc6 = padded[i3m + j3 + k] - padded[i3m + j1 + k];
                *(energy_pos) += (int) ABS(acc4 + acc5 + acc6);
            }

        }

    }

    int jj_limit =  m - K - 9;
    jj_old = jj;


    for (int i = 1; i < i_limit; i++) { //single level reg block calculation 
    	int i1m = (i - 1) * m * 3;
        int i2m = i * m * 3;
        int i3m = (i + 1) * m * 3;
        int im = (i - 1) * (m - 2);
        for (jj = jj_old; jj < jj_limit; jj += 10) {
            short *row0 = padded + i1m + (jj - 1) * 3;
            short *row1 = padded + i2m + (jj - 1) * 3;
            short *row2 = padded + i3m + (jj - 1) * 3;
            int imj = im + jj;

            //loads are all unaligned because we're dealing with shorts
            __m256i r1 = _mm256_loadu_si256((__m256i *) (row0 + 3));
            __m256i r6 = _mm256_loadu_si256((__m256i *) (row2 + 3));
            r1 = _mm256_sub_epi16(r6, r1);
            r1 = _mm256_add_epi16(r1, r1);

            __m256i r3 = _mm256_loadu_si256((__m256i *) row1);
            __m256i r4 = _mm256_loadu_si256((__m256i *) (row1 + 6));
            r3 = _mm256_sub_epi16(r4, r3);
            r3 = _mm256_add_epi16(r3, r3);

            __m256i r0 = _mm256_loadu_si256((__m256i *) row0);
            acc r2;
            r2.intrin = _mm256_loadu_si256((__m256i *) (row0 + 6));
            __m256i r5 = _mm256_loadu_si256((__m256i *) row2);
            r4 = _mm256_sub_epi16(r2.intrin, r0);
            r6 = _mm256_sub_epi16(r5, r0);

            acc r7;
            r7.intrin = _mm256_loadu_si256((__m256i *) (row2 + 6));
            r0 = _mm256_sub_epi16(r7.intrin, r5);
            r5 = _mm256_sub_epi16(r7.intrin, r2.intrin);

            r2.intrin = _mm256_add_epi16(r4, r3);
            r2.intrin = _mm256_add_epi16(r2.intrin, r0);
            r7.intrin = _mm256_add_epi16(r1, r6);
            r7.intrin = _mm256_add_epi16(r7.intrin, r5);

            __m256i r9 = _mm256_loadu_si256((__m256i *) (row0 + 18));
            __m256i r14 = _mm256_loadu_si256((__m256i *) (row2 + 18));
            r9 = _mm256_sub_epi16(r14, r9);
            r9 = _mm256_add_epi16(r9, r9);

            __m256i r11 = _mm256_loadu_si256((__m256i *) (row1 + 15));
            __m256i r12 = _mm256_loadu_si256((__m256i *) (row1 + 21));
            r11 = _mm256_sub_epi16(r12, r11);
            r11 = _mm256_add_epi16(r11, r11);

            __m256i r8 = _mm256_loadu_si256((__m256i *) (row0 + 15));
            acc r10;
            r10.intrin = _mm256_loadu_si256((__m256i *) (row0 + 21));
            __m256i r13 = _mm256_loadu_si256((__m256i *) (row2 + 15));
            r12 = _mm256_sub_epi16(r10.intrin, r8);
            r14 = _mm256_sub_epi16(r13, r8);

            acc r15;
            r15.intrin = _mm256_loadu_si256((__m256i *) (row2 + 21));
            r8 = _mm256_sub_epi16(r15.intrin, r13);
            r13 = _mm256_sub_epi16(r15.intrin, r10.intrin);

            r10.intrin = _mm256_add_epi16(r12, r11);
            r10.intrin = _mm256_add_epi16(r10.intrin, r8);
            r15.intrin = _mm256_add_epi16(r9, r14);
            r15.intrin = _mm256_add_epi16(r15.intrin, r13);

            r2.intrin = _mm256_abs_epi16(r2.intrin);
            r7.intrin = _mm256_abs_epi16(r7.intrin);
            r10.intrin = _mm256_abs_epi16(r10.intrin);
            r15.intrin = _mm256_abs_epi16(r15.intrin);

            short result0 =  r2.data[ 0] +  r2.data[ 1] +  r2.data[ 2] +  r7.data[ 0] +  r7.data[ 1] +  r7.data[ 2];
            short result1 =  r2.data[ 3] +  r2.data[ 4] +  r2.data[ 5] +  r7.data[ 3] +  r7.data[ 4] +  r7.data[ 5];
            short result2 =  r2.data[ 6] +  r2.data[ 7] +  r2.data[ 8] +  r7.data[ 6] +  r7.data[ 7] +  r7.data[ 8];
            short result3 =  r2.data[ 9] +  r2.data[10] +  r2.data[11] +  r7.data[ 9] +  r7.data[10] +  r7.data[11];
            short result4 =  r2.data[12] +  r2.data[13] +  r2.data[14] +  r7.data[12] +  r7.data[13] +  r7.data[14];

            short result5 = r10.data[ 0] + r10.data[ 1] + r10.data[ 2] + r15.data[ 0] + r15.data[ 1] + r15.data[ 2];
            short result6 = r10.data[ 3] + r10.data[ 4] + r10.data[ 5] + r15.data[ 3] + r15.data[ 4] + r15.data[ 5];
            short result7 = r10.data[ 6] + r10.data[ 7] + r10.data[ 8] + r15.data[ 6] + r15.data[ 7] + r15.data[ 8];
            short result8 = r10.data[ 9] + r10.data[10] + r10.data[11] + r15.data[ 9] + r15.data[10] + r15.data[11];
            short result9 = r10.data[12] + r10.data[13] + r10.data[14] + r15.data[12] + r15.data[13] + r15.data[14];

            *(energy + imj - 1) = (int) result0;
            *(energy + imj    ) = (int) result1;
            *(energy + imj + 1) = (int) result2;
            *(energy + imj + 2) = (int) result3;
            *(energy + imj + 3) = (int) result4;
            *(energy + imj + 4) = (int) result5;
            *(energy + imj + 5) = (int) result6;
            *(energy + imj + 6) = (int) result7;
            *(energy + imj + 7) = (int) result8;
            *(energy + imj + 8) = (int) result9;
        }

        //jj_old = jj;

        for(; jj < m - K ; jj++) {
            int acc1;
            int acc2;
            int acc3;
            int acc4;
            int acc5;
            int acc6;

            int j1 = (jj - 1) * 3;
            int j2 = jj * 3;
            int j3 = (jj + 1) * 3;
            int *energy_pos = energy + im + (jj-1);

            // channel R
            //H_y
            acc1 = -(padded[i1m + j1] + ((padded[i1m + j2]) << 1));
            acc2 = padded[i3m + j1] - padded[i1m + j3];
            acc3 = ((padded[i3m + j2]) << 1) + padded[i3m + j3];
            *(energy_pos) = (int) ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = padded[i1m + j3] - padded[i1m + j1];
            acc5 = (padded[i2m + j3] - padded[i2m + j1]) << 1;
            acc6 = padded[i3m + j3] - padded[i3m + j1];
            *(energy_pos) += (int) ABS(acc4 + acc5 + acc6);

            // channel G
            //H_y
            int k = 1;
            acc1 = -(padded[i1m  + j1 + k] + ((padded[i1m + j2 + k]) << 1));
            acc2 = padded[i3m + j1 + k] - padded[i1m + j3 + k];
            acc3 = ((padded[i3m + j2 + k]) << 1) + padded[i3m + j3 + k];
            *(energy_pos) += (int) ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = padded[i1m + j3 + k] - padded[i1m + j1 + k];
            acc5 = (padded[i2m + j3 + k] - padded[i2m + j1 + k]) << 1;
            acc6 = padded[i3m + j3 + k] - padded[i3m + j1 + k];
            *(energy_pos) += (int) ABS(acc4 + acc5 + acc6);

            // channel B
            //H_y
            k = 2;
            acc1 = -(padded[i1m  + j1 + k] + ((padded[i1m + j2 + k]) << 1));
            acc2 = padded[i3m + j1 + k] - padded[i1m + j3 + k];
            acc3 = ((padded[i3m + j2 + k]) << 1) + padded[i3m + j3 + k];
            *(energy_pos) += (int) ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = padded[i1m + j3 + k] - padded[i1m + j1 + k];
            acc5 = (padded[i2m + j3 + k] - padded[i2m + j1 + k]) << 1;
            acc6 = padded[i3m + j3 + k] - padded[i3m + j1 + k];
            *(energy_pos) += (int) ABS(acc4 + acc5 + acc6);
        }

    }
}

	// contains index of the value from the prev row/column from where we came here
	int *backtrack = (int *) malloc(size * sizeof(int)); //different from what we returnCOUNT(mult_count, 2)

	// if horizontal seam -> rotate +90 degrees the energy map
	int *dp;
	if (is_ver == 0) {
		dp = malloc(size * sizeof(int));
	    for (int i = 0; i < rsize; i++) { 
	    	int dp_idx = rsize - i - 1;
	    	int energy_idx = i * csize;
	        for (int j = 0; j < csize; j++) { 
	            dp[j * rsize + dp_idx] = energy[energy_idx + j]; 
	        } 
	    } 

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
	int prev_row_idx;
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
		int prev_row_idx_i = (i-1) * csize - 1;
		for (j = 1; j < column_lim-8; j+=8) {
			where = row + j;
			prev_row_idx = prev_row_idx_i + j;
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
    __m256i minindices = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    __m256i minvalues = _mm256_loadu_si256((__m256i*)(dp+last_row));

	int cnt = 8;
	for (; cnt < csize-7; cnt+=8) {
		idx = _mm256_add_epi32(idx, incr);
		__m256i values = _mm256_loadu_si256((__m256i*)(dp + last_row + cnt));
        __m256i mask = _mm256_cmpgt_epi32(minvalues, values);
        minvalues = _mm256_blendv_epi8(minvalues, values, mask);
        minindices = _mm256_blendv_epi8(minindices, idx, mask);
	}

    // find min index in vector result (in an extremely naive way)
    int32_t values_array[8];
    int32_t indices_array[8];

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

    if (csize < 8) {
    	ret = INT_MAX;
    	cnt = 0;
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