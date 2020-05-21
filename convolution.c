#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse_img.h"
//#include "count.h"
#include <immintrin.h>

#define K 1
#define ABS(X) (((X)<0) ? (-(X)) : (X))
// #define debug   //uncomment for debugging

//--------------------  counter for instructions -------------------

#ifdef count_instr 
extern unsigned long long add_count; //count the total number of add instructions
extern unsigned long long mult_count;  //count the total number of mult instructions
#endif

typedef union acc {
    __m256i intrin;
    short data[16];
} acc;

//------------------------------------------------------------------

// int debug = 1;
//assuming that preprocessing is made of 0 padding 
// Given n rows, m columns of channel padded of some image and the kernel H computes partial gradient corresponding to H given
//padded is of size 3 x n x m
void calc_RGB_energy(int n, int m, short* padded, int* energy){
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    #ifdef count_instr        //counting adds and mults of this function
    unsigned long long count_ifs = 0;        //includes explicit ifs and for loop ifs  -> ADDS
    unsigned long long indexing = 0;         //includes increments of i.j,k variables  -> ADDS
    unsigned long long pointer_adds = 0;     //pointer arithmetic                      -> ADDS
    unsigned long long pointer_mults = 0;    //                                        -> MULTS
    #endif
    int i_limit = n - K;

    int block_width_L1 = 1065;      //working set size is 4(m+2)+m
    int width_limit_L1 = m - K - block_width_L1 + 1;

    int jj, jj_old;
    for(jj = 1; jj < width_limit_L1; jj+= block_width_L1){
        int j_limit = jj + block_width_L1 - 9;

        for(int i = 1 ; i < i_limit ; i++){
            int j;
            for(j = jj ; j < j_limit ; j += 10){
                short *row0 = padded + (i - 1) * m * 3 + (j - 1) * 3;
                short *row1 = padded + (i    ) * m * 3 + (j - 1) * 3;
                short *row2 = padded + (i + 1) * m * 3 + (j - 1) * 3;

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

                *(energy + (i - 1) * (m - 2) + j - 1) = (int) result0;
                *(energy + (i - 1) * (m - 2) + j    ) = (int) result1;
                *(energy + (i - 1) * (m - 2) + j + 1) = (int) result2;
                *(energy + (i - 1) * (m - 2) + j + 2) = (int) result3;
                *(energy + (i - 1) * (m - 2) + j + 3) = (int) result4;
                *(energy + (i - 1) * (m - 2) + j + 4) = (int) result5;
                *(energy + (i - 1) * (m - 2) + j + 5) = (int) result6;
                *(energy + (i - 1) * (m - 2) + j + 6) = (int) result7;
                *(energy + (i - 1) * (m - 2) + j + 7) = (int) result8;
                *(energy + (i - 1) * (m - 2) + j + 8) = (int) result9;
            }
            
            for(; j < jj + block_width_L1 ; j++) {
                int acc1;
                int acc2;
                int acc3;
                int acc4;
                int acc5;
                int acc6;
                // channel R
                //H_y
                int k = 0;
                acc1 = -(padded[(i - 1) * m * 3  + (j - 1)*3 + k] + ((padded[(i - 1) * m *3 + j * 3 + k]) << 1));
                acc2 = padded[(i + 1) * m * 3 + (j - 1) * 3 + k] - padded[(i - 1) * m * 3 + (j + 1) * 3 + k];
                acc3 = ((padded[(i + 1) * m * 3 + j*3 + k]) << 1) + padded[(i + 1) * m * 3 + (j + 1) * 3 + k];
                *(energy + (i-1)*(m-2) + (j-1)) = (int) ABS(acc1 + acc2 + acc3);
                //H_x
                acc4 = padded[(i - 1) * m * 3 + (j + 1) * 3 + k] - padded[(i - 1) * m * 3 + (j - 1) * 3 + k];
                acc5 = (padded[i * m  * 3+ (j + 1) * 3 + k] - padded[i * m * 3 + (j - 1) * 3 + k]) << 1;
                acc6 = padded[(i + 1) * m * 3 + (j + 1) * 3 + k] - padded[(i + 1) * m * 3 + (j - 1) * 3 + k];
                *(energy + (i-1)*(m-2) + (j-1)) += (int) ABS(acc4 + acc5 + acc6);

                // channel G
                //H_y
                k = 1;
                acc1 = -(padded[(i - 1) * m * 3  + (j - 1)*3 + k] + ((padded[(i - 1) * m *3 + j * 3 + k]) << 1));
                acc2 = padded[(i + 1) * m * 3 + (j - 1) * 3 + k] - padded[(i - 1) * m * 3 + (j + 1) * 3 + k];
                acc3 = ((padded[(i + 1) * m * 3 + j*3 + k]) << 1) + padded[(i + 1) * m * 3 + (j + 1) * 3 + k];
                *(energy + (i-1)*(m-2) + (j-1)) += (int) ABS(acc1 + acc2 + acc3);
                //H_x
                acc4 = padded[(i - 1) * m * 3 + (j + 1) * 3 + k] - padded[(i - 1) * m * 3 + (j - 1) * 3 + k];
                acc5 = (padded[i * m  * 3+ (j + 1) * 3 + k] - padded[i * m * 3 + (j - 1) * 3 + k]) << 1;
                acc6 = padded[(i + 1) * m * 3 + (j + 1) * 3 + k] - padded[(i + 1) * m * 3 + (j - 1) * 3 + k];
                *(energy + (i-1)*(m-2) + (j-1)) += (int) ABS(acc4 + acc5 + acc6);

                // channel B
                //H_y
                k = 2;
                acc1 = -(padded[(i - 1) * m * 3  + (j - 1)*3 + k] + ((padded[(i - 1) * m *3 + j * 3 + k]) << 1));
                acc2 = padded[(i + 1) * m * 3 + (j - 1) * 3 + k] - padded[(i - 1) * m * 3 + (j + 1) * 3 + k];
                acc3 = ((padded[(i + 1) * m * 3 + j*3 + k]) << 1) + padded[(i + 1) * m * 3 + (j + 1) * 3 + k];
                *(energy + (i-1)*(m-2) + (j-1)) += (int) ABS(acc1 + acc2 + acc3);
                //H_x
                acc4 = padded[(i - 1) * m * 3 + (j + 1) * 3 + k] - padded[(i - 1) * m * 3 + (j - 1) * 3 + k];
                acc5 = (padded[i * m  * 3+ (j + 1) * 3 + k] - padded[i * m * 3 + (j - 1) * 3 + k]) << 1;
                acc6 = padded[(i + 1) * m * 3 + (j + 1) * 3 + k] - padded[(i + 1) * m * 3 + (j - 1) * 3 + k];
                *(energy + (i-1)*(m-2) + (j-1)) += (int) ABS(acc4 + acc5 + acc6);
            }

        }

    }

    int jj_limit =  m - K - 9;
    jj_old = jj;


    for (int i = 1; i < i_limit; i++) { //single level reg block calculation 
        for (jj = jj_old; jj < jj_limit; jj += 10) {
            short *row0 = padded + (i - 1) * m * 3 + (jj - 1) * 3;
            short *row1 = padded + (i    ) * m * 3 + (jj - 1) * 3;
            short *row2 = padded + (i + 1) * m * 3 + (jj - 1) * 3;

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

            *(energy + (i - 1) * (m - 2) + jj - 1) = (int) result0;
            *(energy + (i - 1) * (m - 2) + jj    ) = (int) result1;
            *(energy + (i - 1) * (m - 2) + jj + 1) = (int) result2;
            *(energy + (i - 1) * (m - 2) + jj + 2) = (int) result3;
            *(energy + (i - 1) * (m - 2) + jj + 3) = (int) result4;
            *(energy + (i - 1) * (m - 2) + jj + 4) = (int) result5;
            *(energy + (i - 1) * (m - 2) + jj + 5) = (int) result6;
            *(energy + (i - 1) * (m - 2) + jj + 6) = (int) result7;
            *(energy + (i - 1) * (m - 2) + jj + 7) = (int) result8;
            *(energy + (i - 1) * (m - 2) + jj + 8) = (int) result9;
        }

        //jj_old = jj;

        for(; jj < m - K ; jj++) {
            int acc1;
            int acc2;
            int acc3;
            int acc4;
            int acc5;
            int acc6;
            // channel R
            //H_y
            int k = 0;
            acc1 = -(padded[(i - 1) * m * 3  + (jj - 1)*3 + k] + ((padded[(i - 1) * m *3 + jj * 3 + k]) << 1));
            acc2 = padded[(i + 1) * m * 3 + (jj - 1) * 3 + k] - padded[(i - 1) * m * 3 + (jj + 1) * 3 + k];
            acc3 = ((padded[(i + 1) * m * 3 + jj*3 + k]) << 1) + padded[(i + 1) * m * 3 + (jj + 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (jj-1)) = (int) ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = padded[(i - 1) * m * 3 + (jj + 1) * 3 + k] - padded[(i - 1) * m * 3 + (jj - 1) * 3 + k];
            acc5 = (padded[i * m  * 3+ (jj + 1) * 3 + k] - padded[i * m * 3 + (jj - 1) * 3 + k]) << 1;
            acc6 = padded[(i + 1) * m * 3 + (jj + 1) * 3 + k] - padded[(i + 1) * m * 3 + (jj - 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (jj-1)) += (int) ABS(acc4 + acc5 + acc6);

            // channel G
            //H_y
            k = 1;
            acc1 = -(padded[(i - 1) * m * 3  + (jj - 1)*3 + k] + ((padded[(i - 1) * m *3 + jj * 3 + k]) << 1));
            acc2 = padded[(i + 1) * m * 3 + (jj - 1) * 3 + k] - padded[(i - 1) * m * 3 + (jj + 1) * 3 + k];
            acc3 = ((padded[(i + 1) * m * 3 + jj*3 + k]) << 1) + padded[(i + 1) * m * 3 + (jj + 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (jj-1)) += (int) ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = padded[(i - 1) * m * 3 + (jj + 1) * 3 + k] - padded[(i - 1) * m * 3 + (jj - 1) * 3 + k];
            acc5 = (padded[i * m  * 3+ (jj + 1) * 3 + k] - padded[i * m * 3 + (jj - 1) * 3 + k]) << 1;
            acc6 = padded[(i + 1) * m * 3 + (jj + 1) * 3 + k] - padded[(i + 1) * m * 3 + (jj - 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (jj-1)) += (int) ABS(acc4 + acc5 + acc6);

            // channel B
            //H_y
            k = 2;
            acc1 = -(padded[(i - 1) * m * 3  + (jj - 1)*3 + k] + ((padded[(i - 1) * m *3 + jj * 3 + k]) << 1));
            acc2 = padded[(i + 1) * m * 3 + (jj - 1) * 3 + k] - padded[(i - 1) * m * 3 + (jj + 1) * 3 + k];
            acc3 = ((padded[(i + 1) * m * 3 + jj*3 + k]) << 1) + padded[(i + 1) * m * 3 + (jj + 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (jj-1)) += (int) ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = padded[(i - 1) * m * 3 + (jj + 1) * 3 + k] - padded[(i - 1) * m * 3 + (jj - 1) * 3 + k];
            acc5 = (padded[i * m  * 3+ (jj + 1) * 3 + k] - padded[i * m * 3 + (jj - 1) * 3 + k]) << 1;
            acc6 = padded[(i + 1) * m * 3 + (jj + 1) * 3 + k] - padded[(i + 1) * m * 3 + (jj - 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (jj-1)) += (int) ABS(acc4 + acc5 + acc6);
        }

    }
    


    #ifdef debug 
    char *fname = "energy_map.png";
    save_as_grayscale_image(fname, m-2, n-2, energy);
    printf("Saved first energy map as %s\n", fname);
    // debug = 1;
    #endif

    #ifdef count_instr 
    count_ifs += n-1 + (n-2)*(m-1) + (n-2)*(m-2)*4 + (n-2)*(m-2)*3*4; //count lines 46-52
    indexing += n-2 + (n-2)*(m-2) + (n-2)*(m-2)*3 + (n-2)*(m-2)*3*3;
    pointer_adds += (n-2)*(m-2)*3*3*(2 + 2 + 3);
    pointer_mults += (n-2)*(m-2)*3*3*2;
    add_count += (n-2)*(m-2)*3*3;                              //count directly the add and mult in line 50
    mult_count += (n-2)*(m-2)*3*3;

    //count total
    add_count += count_ifs + indexing + pointer_adds; 
    mult_count += pointer_mults;
    printf("NO ADDS paddedOR calc_energy IS: %llu \n", add_count); 
    printf("NO MULTS paddedOR calc_energy IS: %llu \n", mult_count); 
    #endif
}


//first method called before anything else to create a 0 frame aronund original image
//make sure to free the returned pointer after the first seam is found and removed 
//repeat for each seam  
short* padd0_image(int n, int m, unsigned char* channels){

  int size = (n+2) * (m+2) * 3;
  short* padded_image = (short*) malloc( size*sizeof(short));

  //int padded_image[n+2][m+2][3];
  for(int i = 0 ; i < n+2 ; i++){
    for(int j = 0 ; j < m+2 ; j++){
      for(int k = 0 ; k < 3 ; k++){
        //if the column is 0 or m+1 or the row is 0 or n+1 we set 0 otherwise copy the value 
        if(i == 0 || j == 0 || i == n+1 || j == m+1){
          //padded_image[i*(n+2)*(m+2) + (m+2)*j + k] = 0;
          padded_image[i*(m+2)*3 + j*3 + k] = 0;
        } else{
          //printf("in pad its %f\n", channels[i*n*m + m*j + k]);
          //padded_image[i*(n+2)*(m+2) + (m+2)*j + k] = channels[i*n*m + m*(j-1) + k-1];
          padded_image[i*(m+2)*3 + j*3 + k] = (short) (channels[(i-1)*m*3 + (j-1)*3 + k]);
        }
    }
  }
}
  return padded_image;
}

