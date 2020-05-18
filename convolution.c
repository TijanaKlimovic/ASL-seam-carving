#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse_img.h"
// #include "count.h"

#define K 1
#define ABS(X)(((X) < 0) ? (-(X)) : (X))
// #define debug   //uncomment for debugging

//--------------------  counter for instructions -------------------

#ifdef count_instr
extern unsigned long long add_count; //count the total number of add instructions
extern unsigned long long mult_count; //count the total number of mult instructions
#endif

//------------------------------------------------------------------

// int debug = 1;
//assuming that preprocessing is made of 0 padding 
// Given n rows, m columns of channel F of some image and the kernel H computes partial gradient corresponding to H given
//F is of size 3 x n x m
void calc_energy_1(int n, int m, int * F, int * part_grad) {
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    #ifdef count_instr //counting adds and mults of this function
    unsigned long long count_ifs = 0; //includes explicit ifs and for loop ifs  -> ADDS
    unsigned long long indexing = 0; //includes increments of i.j,k variables  -> ADDS
    unsigned long long pointer_adds = 0; //pointer arithmetic                      -> ADDS
    unsigned long long pointer_mults = 0; //                                        -> MULTS
    #endif

    int block_height_L3 = n - K;
    int block_height_L2 = n - K;
    int block_height_L1 = n - K;
    int block_height_r = n - K;

    int block_width_L3 = 99998; //currently going channel by channel maybe later we can parallelise it
    int block_width_L2 = 12798;
    int block_width_L1 = 1400;//1598;
    int block_width_r = 1;  //working set size is 4(m+2)+m

    int width_limit_L3 = m - K - block_width_L3;
    int width_limit_L2;
    int width_limit_L1 = m - K - block_width_L1;
    int width_limit_r;

    int j, i;
    for (j = 1; j < width_limit_L1; j = j + block_width_L1) {
        for (i = 1; i < n - K; i += block_height_L1) {

            width_limit_r = j + block_width_L1 - block_width_r; //take care in case block_width_L2 doesnt divide width_limit_L2
            int j_L3, i_L3;

            for (j_L3 = j; j_L3 < j + block_width_L1; j_L3 += block_width_r) {
                for (i_L3 = i; i_L3 < i + block_height_L1; i_L3 += block_height_r) {
                    for (int ii = i_L3; ii < i_L3 + block_height_r; ii++) {
                        for (int jj = j_L3; jj < j_L3 + block_width_r; jj++) {

                            int acc1;
                            int acc2;
                            int acc3;
                            int acc4;
                            int acc5;
                            int acc6;
                            //fixed kernels 
                            // int H_y[3][3] = {
                            //   {-1,-2,-1},
                            //   {0,0,0},
                            //   {1,2,1}};

                            // int H_x[3][3] = {
                            //   {-1,0,1},
                            //   {-2,0,2},
                            //   {-1,0,1}};

                            int nw, n, ne, w, e, sw, s, se;
                            nw = F[(ii - 1) * m + (jj - 1)];
                            n = F[(ii - 1) * m + jj];
                            ne = F[(ii - 1) * m + jj + 1];
                            w = F[ii * m + jj - 1];
                            e = F[ii * m + jj + 1];
                            sw = F[(ii + 1) * m + jj - 1];
                            s = F[(ii + 1) * m + jj];
                            se = F[(ii + 1) * m + jj + 1];
                            //H_y
                            acc1 = (s - n) << 1;
                            acc2 = se - nw;
                            acc3 = sw - ne;

                            // acc1 = -(nw + (n << 1));
                            // acc2 = sw - ne;
                            // acc3 = (s << 1) + se;
                            //H_x

                            acc4 = ne - nw;
                            acc5 = (e - w) << 1;
                            acc6 = se - sw;
                            *(part_grad + ii * m + jj) = ABS(acc1 + acc2 + acc3) + ABS(acc4 + acc5 + acc6);
                            #ifdef count_instr //count line 55
                            mult_count++;
                            pointer_adds += 2 * 2; //assuming worse case that the after having the pointer arithmetic done its not saved but redone
                            pointer_mults += 2;
                            #endif
                        }
                    }
                }

                #ifdef count_instr //count line 54
                count_ifs++; //when not taken check must be made too
                pointer_adds += 2; //need to dereference first to compare 
                pointer_mults++;
                #endif
            }

            //bad non cooparating js
            for (; j_L3 < j + block_width_L1; j_L3++) {
                for (i_L3 = 1; i_L3 < i + block_height_L1; i_L3++) {
                    int acc1;
                    int acc2;
                    int acc3;
                    int acc4;
                    int acc5;
                    int acc6;
                    int nw, n, ne, w, e, sw, s, se;

                    nw = F[(i_L3 - 1) * m + (j_L3 - 1)];
                    n = F[(i_L3 - 1) * m + j_L3];
                    ne = F[(i_L3 - 1) * m + j_L3 + 1];
                    w = F[i_L3 * m + j_L3 - 1];
                    e = F[i_L3 * m + j_L3 + 1];
                    sw = F[(i_L3 + 1) * m + j_L3 - 1];
                    s = F[(i_L3 + 1) * m + j_L3];
                    se = F[(i_L3 + 1) * m + j_L3 + 1];

                    acc1 = (s - n) << 1;
                    acc2 = se - nw;
                    acc3 = sw - ne;

                    acc4 = ne - nw;
                    acc5 = (e - w) << 1;
                    acc6 = se - sw;

                    *(part_grad + i_L3 * m + j_L3) = ABS(acc1 + acc2 + acc3) + ABS(acc4 + acc5 + acc6);
                }
            }

        }
    }

    //bad non cooparating js
    for (; j < m - K; j++) {
        for (i = 1; i < n - K; i++) {
            int acc1;
            int acc2;
            int acc3;
            int acc4;
            int acc5;
            int acc6;
            //H_y
            acc1 = -(F[(i - 1) * m + (j - 1)] + ((F[(i - 1) * m + j]) << 1));
            acc2 = F[(i + 1) * m + (j - 1)] - F[(i - 1) * m + j + 1];
            acc3 = ((F[(i + 1) * m + j]) << 1) + F[(i + 1) * m + j + 1];
            *(part_grad + i * m + j) = ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = F[(i - 1) * m + j + 1] - F[(i - 1) * m + (j - 1)];
            acc5 = (F[i * m + j + 1] - F[i * m + j - 1]) << 1;
            acc6 = F[(i + 1) * m + j + 1] - F[(i + 1) * m + j - 1];
            *(part_grad + i * m + j) += ABS(acc4 + acc5 + acc6);

        }
    }

    #ifdef count_instr
    count_ifs += n - 1 + (n - 2) * (m - 1) + (n - 2) * (m - 2) * 4 + (n - 2) * (m - 2) * 3 * 4; //count lines 46-52
    indexing += n - 2 + (n - 2) * (m - 2) + (n - 2) * (m - 2) * 3 + (n - 2) * (m - 2) * 3 * 3;
    pointer_adds += (n - 2) * (m - 2) * 3 * 3 * (2 + 2 + 3);
    pointer_mults += (n - 2) * (m - 2) * 3 * 3 * 2;
    add_count += (n - 2) * (m - 2) * 3 * 3; //count directly the add and mult in line 50
    mult_count += (n - 2) * (m - 2) * 3 * 3;

    //count total
    add_count += count_ifs + indexing + pointer_adds;
    mult_count += pointer_mults;
    printf("NO ADDS FOR calc_energy IS: %llu \n", add_count);
    printf("NO MULTS FOR calc_energy IS: %llu \n", mult_count);
    #endif
}


//version with no accumulators only accounting for the 
void calc_energy(int n, int m, int * F, int * part_grad) {
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    #ifdef count_instr //counting adds and mults of this function
    unsigned long long count_ifs = 0; //includes explicit ifs and for loop ifs  -> ADDS
    unsigned long long indexing = 0; //includes increments of i.j,k variables  -> ADDS
    unsigned long long pointer_adds = 0; //pointer arithmetic                      -> ADDS
    unsigned long long pointer_mults = 0; //                                        -> MULTS
    #endif

    //int block_height_L3 = n - K;
    //int block_height_L2 = n - K;
    int block_height_L1 = n - K;
    int block_height_R = n - K;

    //int block_width_L3 = 166665; //currently going channel by channel maybe later we can parallelise it 166667
    //int block_width_L2 = 21332; // 21334
    int block_width_L1 = 1598;
    int block_width_R = 2;  //working set size is 4(m+2)+m

    //int width_limit_L3 = m - K - block_width_L3;
    //int width_limit_L2;
    int width_limit_L1 = m - K - block_width_L1 + 1;
    int width_limit_R ; //m - K - block_width_R + 1

    int j_old;

    int j, i;
    for (j = 1; j < width_limit_L1; j = j + block_width_L1) {
        width_limit_R = j + block_width_L1 - block_width_R + 1; //take care in case block_width_L2 doesnt divide width_limit_L2
        int j_L3;

        for (int ii = 1; ii < n - K; ii++) { //single level 2 block calculation 
            for (j_L3 = j; j_L3 < width_limit_R; j_L3 += block_width_R) { //single level 3 block calculation
            
                for (int jj = j_L3; jj < j_L3 + block_width_R; jj++) {
                    int acc;
                    //H_y
                    *(part_grad + ii * m + jj) = (-(F[(ii - 1) * m + (jj - 1)] + ((F[(ii - 1) * m + jj]) << 1))) + (F[(ii + 1) * m + (jj - 1)] - F[(ii - 1) * m + jj + 1]) + (((F[(ii + 1) * m + jj]) << 1) + F[(ii + 1) * m + jj + 1]);
                    *(part_grad + ii * m + jj) = ABS(*(part_grad + ii * m + jj));
                    //H_x
                    acc = (F[(ii - 1) * m + jj + 1] - F[(ii - 1) * m + (jj - 1)]) + (((F[ii * m + jj + 1] - F[ii * m + jj - 1]) << 1)) + (F[(ii + 1) * m + jj + 1] - F[(ii + 1) * m + jj - 1]);
                    *(part_grad + ii * m + jj) += ABS(acc);
                    #ifdef count_instr //count line 55
                    mult_count++;
                    pointer_adds += 2 * 2; //assuming worse case that the after having the pointer arithmetic done its not saved but redone
                    pointer_mults += 2;
                    #endif
                }
            }

            //bad non cooparating jjs
            for (; j_L3 < j + block_width_L1; j_L3++) {
                int acc;
                //H_y
                *(part_grad + ii * m + j_L3) = (-(F[(ii - 1) * m + (j_L3 - 1)] + ((F[(ii - 1) * m + j_L3]) << 1))) + (F[(ii + 1) * m + (j_L3 - 1)] - F[(ii - 1) * m + j_L3 + 1]) + (((F[(ii + 1) * m + j_L3]) << 1) + F[(ii + 1) * m + j_L3 + 1]);
                *(part_grad + ii * m + j_L3) = ABS(*(part_grad + ii * m + j_L3));
                //H_x
                acc = (F[(ii - 1) * m + j_L3 + 1] - F[(ii - 1) * m + (j_L3 - 1)]) + (((F[ii * m + j_L3 + 1] - F[ii * m + j_L3 - 1]) << 1)) + (F[(ii + 1) * m + j_L3 + 1] - F[(ii + 1) * m + j_L3 - 1]);
                *(part_grad + ii * m + j_L3) += ABS(acc);
            }
            
        }

        #ifdef count_instr //count line 54
        count_ifs++; //when not taken check must be made too
        pointer_adds += 2; //need to dereference first to compare 
        pointer_mults++;
        #endif

    }

    //bad non cooparating js
    int j_limit =  m - K - block_width_R + 1;
    j_old = j;

    for (int ii = 1; ii < n - K; ii++) { //single level 2 block calculation 
        for (j = j_old; j < j_limit; j = j + block_width_R) {
                for (int jj = j; jj < j + block_width_R; jj++) {
                    int acc;
                    //H_y
                    *(part_grad + ii * m + jj) = (-(F[(ii - 1) * m + (jj - 1)] + ((F[(ii - 1) * m + jj]) << 1))) + (F[(ii + 1) * m + (jj - 1)] - F[(ii - 1) * m + jj + 1]) + (((F[(ii + 1) * m + jj]) << 1) + F[(ii + 1) * m + jj + 1]);
                    *(part_grad + ii * m + jj) = ABS(*(part_grad + ii * m + jj));
                    //H_x
                    acc = (F[(ii - 1) * m + jj + 1] - F[(ii - 1) * m + (jj - 1)]) + (((F[ii * m + jj + 1] - F[ii * m + jj - 1]) << 1)) + (F[(ii + 1) * m + jj + 1] - F[(ii + 1) * m + jj - 1]);
                    *(part_grad + ii * m + jj) += ABS(acc);
                    #ifdef count_instr //count line 55
                    mult_count++;
                    pointer_adds += 2 * 2; //assuming worse case that the after having the pointer arithmetic done its not saved but redone
                    pointer_mults += 2;
                    #endif
                }
            }        
    }

    j_old = j;

    for (int ii = 1; ii < n - K; ii++) { //single level 2 block calculation 
         for (j = j_old; j < m - K; j++) {
                    int acc;
                    //H_y
                    *(part_grad + ii * m + j) = -((F[(ii - 1) * m + (j - 1)] + ((F[(ii - 1) * m + j]) << 1))) + (F[(ii + 1) * m + (j - 1)] - F[(ii - 1) * m + j + 1]) + (((F[(ii + 1) * m + j]) << 1) + F[(ii + 1) * m + j + 1]);
                    *(part_grad + ii * m + j) = ABS(*(part_grad + ii * m + j));
                    //H_x
                    acc = (F[(ii - 1) * m + j + 1] - F[(ii - 1) * m + (j - 1)]) + (((F[ii * m + j+ 1] - F[ii * m + j - 1]) << 1)) + (F[(ii + 1) * m + j + 1] - F[(ii + 1) * m + j - 1]);
                    *(part_grad + ii * m + j) += ABS(acc);
                    #ifdef count_instr //count line 55
                    mult_count++;
                    pointer_adds += 2 * 2; //assuming worse case that the after having the pointer arithmetic done its not saved but redone
                    pointer_mults += 2;
                    #endif
            }        
    }

    #ifdef count_instr
    count_ifs += n - 1 + (n - 2) * (m - 1) + (n - 2) * (m - 2) * 4 + (n - 2) * (m - 2) * 3 * 4; //count lines 46-52
    indexing += n - 2 + (n - 2) * (m - 2) + (n - 2) * (m - 2) * 3 + (n - 2) * (m - 2) * 3 * 3;
    pointer_adds += (n - 2) * (m - 2) * 3 * 3 * (2 + 2 + 3);
    pointer_mults += (n - 2) * (m - 2) * 3 * 3 * 2;
    add_count += (n - 2) * (m - 2) * 3 * 3; //count directly the add and mult in line 50
    mult_count += (n - 2) * (m - 2) * 3 * 3;

    //count total
    add_count += count_ifs + indexing + pointer_adds;
    mult_count += pointer_mults;
    printf("NO ADDS FOR calc_energy IS: %llu \n", add_count);
    printf("NO MULTS FOR calc_energy IS: %llu \n", mult_count);
    #endif
}

//version with accumulators and scalar replacement
void calc_energy_acc(int n, int m, int * F, int * part_grad) {
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    #ifdef count_instr //counting adds and mults of this function
    unsigned long long count_ifs = 0; //includes explicit ifs and for loop ifs  -> ADDS
    unsigned long long indexing = 0; //includes increments of i.j,k variables  -> ADDS
    unsigned long long pointer_adds = 0; //pointer arithmetic                      -> ADDS
    unsigned long long pointer_mults = 0; //                                        -> MULTS
    #endif

    //int block_height_L3 = n - K;
    //int block_height_L2 = n - K;
    int block_height_L1 = n - K;
    int block_height_R = n - K;

    //int block_width_L3 = 166665; //currently going channel by channel maybe later we can parallelise it 166667
    //int block_width_L2 = 21332; // 21334
    int block_width_L1 = 1598;
    int block_width_R = 1;  //working set size is 4(m+2)+m

    //int width_limit_L3 = m - K - block_width_L3;
    //int width_limit_L2;
    int width_limit_L1 = m - K - block_width_L1 + 1;
    int width_limit_R ; //m - K - block_width_R + 1

    int j, i;
    for (j = 1; j < width_limit_L1; j = j + block_width_L1) {
        width_limit_R = j + block_width_L1 - block_width_R + 1; //take care in case block_width_L2 doesnt divide width_limit_L2
        int j_L3;

        for (int ii = 1; ii < n - K; ii++) { //single level 2 block calculation 
            for (j_L3 = j; j_L3 < width_limit_R; j_L3 += block_width_R) { //single level 3 block calculation
            
                for (int jj = j_L3; jj < j_L3 + block_width_R; jj++) {
                    int acc1;
                    int acc2;
                    int acc3;
                    int acc4;
                    int acc5;
                    int acc6;
                    int nw, n, ne, w, e, sw, s, se;
                    nw = F[(ii - 1) * m + (jj - 1)];
                    n = F[(ii - 1) * m + jj];
                    ne = F[(ii - 1) * m + jj + 1];
                    w = F[ii * m + jj - 1];
                    e = F[ii * m + jj + 1];
                    sw = F[(ii + 1) * m + jj - 1];
                    s = F[(ii + 1) * m + jj];
                    se = F[(ii + 1) * m + jj + 1];
                    //H_y
                    acc1 = (s - n) << 1;
                    acc2 = se - nw;
                    acc3 = sw - ne;

                    // acc1 = -(nw + (n << 1));
                    // acc2 = sw - ne;
                    // acc3 = (s << 1) + se;
                    //H_x

                    acc4 = ne - nw;
                    acc5 = (e - w) << 1;
                    acc6 = se - sw;
                    *(part_grad + ii * m + jj) = ABS(acc1 + acc2 + acc3) + ABS(acc4 + acc5 + acc6);
                }
            }
        }

        #ifdef count_instr //count line 54
        count_ifs++; //when not taken check must be made too
        pointer_adds += 2; //need to dereference first to compare 
        pointer_mults++;
        #endif

        //bad non cooparating jjs
        for (; j_L3 < j + block_width_L1; j_L3++) {
            for (int i_L3 = 1; i_L3 < n - K; i_L3++) {
                int acc1;
                int acc2;
                int acc3;
                int acc4;
                int acc5;
                int acc6;
                int nw, n, ne, w, e, sw, s, se;
                nw = F[(i_L3 - 1) * m + (j_L3 - 1)];
                n = F[(i_L3 - 1) * m + j_L3];
                ne = F[(i_L3 - 1) * m + j_L3 + 1];
                w = F[i_L3 * m + j_L3 - 1];
                e = F[i_L3 * m + j_L3 + 1];
                sw = F[(i_L3 + 1) * m + j_L3 - 1];
                s = F[(i_L3 + 1) * m + j_L3];
                se = F[(i_L3 + 1) * m + j_L3 + 1];
                //H_y
                acc1 = (s - n) << 1;
                acc2 = se - nw;
                acc3 = sw - ne;

                // acc1 = -(nw + (n << 1));
                // acc2 = sw - ne;
                // acc3 = (s << 1) + se;
                //H_x

                acc4 = ne - nw;
                acc5 = (e - w) << 1;
                acc6 = se - sw;
                *(part_grad + i_L3 * m + j_L3) = ABS(acc1 + acc2 + acc3) + ABS(acc4 + acc5 + acc6);
            }
        }
    }

    //bad non cooparating js
    int j_limit =  m - K - block_width_R + 1;

    for (int ii = 1; ii < n - K; ii++) { //single level 2 block calculation 
        for (; j < j_limit; j = j + block_width_R) {
                for (int jj = j; jj < j + block_width_R; jj++) {
                    int acc1;
                    int acc2;
                    int acc3;
                    int acc4;
                    int acc5;
                    int acc6;
                    int nw, n, ne, w, e, sw, s, se;
                    nw = F[(ii - 1) * m + (jj - 1)];
                    n = F[(ii - 1) * m + jj];
                    ne = F[(ii - 1) * m + jj + 1];
                    w = F[ii * m + jj - 1];
                    e = F[ii * m + jj + 1];
                    sw = F[(ii + 1) * m + jj - 1];
                    s = F[(ii + 1) * m + jj];
                    se = F[(ii + 1) * m + jj + 1];
                    //H_y
                    acc1 = (s - n) << 1;
                    acc2 = se - nw;
                    acc3 = sw - ne;

                    // acc1 = -(nw + (n << 1));
                    // acc2 = sw - ne;
                    // acc3 = (s << 1) + se;
                    //H_x

                    acc4 = ne - nw;
                    acc5 = (e - w) << 1;
                    acc6 = se - sw;
                    *(part_grad + ii * m + jj) = ABS(acc1 + acc2 + acc3) + ABS(acc4 + acc5 + acc6);
                }
            }        
    }
    for (int ii = 1; ii < n - K; ii++) { //single level 2 block calculation 
         for (; j < m - K; j++) {
                int acc1;
                int acc2;
                int acc3;
                int acc4;
                int acc5;
                int acc6;
                int nw, n, ne, w, e, sw, s, se;
                nw = F[(ii - 1) * m + (j - 1)];
                n = F[(ii - 1) * m + j];
                ne = F[(ii - 1) * m + j + 1];
                w = F[ii * m + j - 1];
                e = F[ii * m + j + 1];
                sw = F[(ii + 1) * m + j - 1];
                s = F[(ii + 1) * m + j];
                se = F[(ii + 1) * m + j + 1];
                //H_y
                acc1 = (s - n) << 1;
                acc2 = se - nw;
                acc3 = sw - ne;

                // acc1 = -(nw + (n << 1));
                // acc2 = sw - ne;
                // acc3 = (s << 1) + se;
                //H_x

                acc4 = ne - nw;
                acc5 = (e - w) << 1;
                acc6 = se - sw;
                *(part_grad + ii * m + j) = ABS(acc1 + acc2 + acc3) + ABS(acc4 + acc5 + acc6);
            }        
    }

    #ifdef count_instr
    count_ifs += n - 1 + (n - 2) * (m - 1) + (n - 2) * (m - 2) * 4 + (n - 2) * (m - 2) * 3 * 4; //count lines 46-52
    indexing += n - 2 + (n - 2) * (m - 2) + (n - 2) * (m - 2) * 3 + (n - 2) * (m - 2) * 3 * 3;
    pointer_adds += (n - 2) * (m - 2) * 3 * 3 * (2 + 2 + 3);
    pointer_mults += (n - 2) * (m - 2) * 3 * 3 * 2;
    add_count += (n - 2) * (m - 2) * 3 * 3; //count directly the add and mult in line 50
    mult_count += (n - 2) * (m - 2) * 3 * 3;

    //count total
    add_count += count_ifs + indexing + pointer_adds;
    mult_count += pointer_mults;
    printf("NO ADDS FOR calc_energy IS: %llu \n", add_count);
    printf("NO MULTS FOR calc_energy IS: %llu \n", mult_count);
    #endif
}

//assumes channels is padded with 0 of size 3 x n x m
//assumes result is of size n-2 x m-2 
//calculates the cumulative sum over the individual channel energies
//returns the energy map result over color image of size  n-2 x m-2 
void calc_RGB_energy(int n, int m, int * channels, int * result) {

    #ifdef count_instr //counting adds and mults of this function
    unsigned long long count_ifs = 0; //includes explicit ifs and for loop ifs  -> ADDS
    unsigned long long indexing = 0; //includes increments of i.j,k variables  -> ADDS
    unsigned long long pointer_adds = 0; //pointer arithmetic                      -> ADDS
    unsigned long long pointer_mults = 0; //                                        -> MULTS
    #endif

    //fixed kernels 
    // int H_y[3][3] = {
    //   {-1,-2,-1},
    //   {0,0,0},
    //   {1,2,1}};

    // int H_x[3][3] = {
    //   {-1,0,1},
    //   {-2,0,2},
    //   {-1,0,1}};

    int size = 3 * n * m;

    int * partial = (int * ) malloc(size * sizeof(int));

    //calculate the parital derivatives 
    for (int i = 0; i < 3; i++) {
        //pass the ith channel for energy calculation
        calc_energy(n, m, channels + n * m * i, partial + n * m * i);
    }

    #ifdef count_instr //counts lines 120-124
    count_ifs += 4;
    indexing += 3;
    pointer_adds += 3 * 4;
    pointer_mults += 3 * 8;
    #endif

    for (int i = 0; i < n - 2; i++) {
        for (int j = 0; j < m - 2; j++) {
            result[(m - 2) * i + j] = 0;
        }
    }

    #ifdef count_instr //counts lines 134-138
    count_ifs += n - 1 + (n - 2) * (m - 1);
    indexing += n - 2 + (n - 2) * (m - 2);
    pointer_adds += (n - 2) * (m - 2) * 2;
    pointer_mults += (n - 2) * (m - 2);
    #endif

    //calculate the total 3d energy 

    for (int i = 0; i < 3; i++) {
        for (int j = 1; j < n - 1; j++) {
            for (int k = 1; k < m - 1; k++) {
                //add elementwise along the z axis 
                *(result + (m - 2) * (j - 1) + k - 1) += * (partial + i * m * n + j * m + k);
            }
        }
    }

    #ifdef count_instr //counts lines 134-138
    count_ifs += n - 1 + (n - 2) * (m - 1) + (n - 2) * (m - 2) * 4;
    indexing += n - 2 + (n - 2) * (m - 2) + (n - 2) * (m - 2) * 3;
    pointer_adds += (n - 2) * (m - 2) * 3 * (5 + 3 + 3);
    pointer_mults += (n - 2) * (m - 2) * 3 * 7;
    add_count += (n - 2) * (m - 2) * 3 * 2; //count directly the adds in line 154

    //count total
    add_count += count_ifs + indexing + pointer_adds;
    mult_count += pointer_mults;
    printf("NO ADDS FOR calc_energy IS: %llu \n", add_count);
    printf("NO MULTS FOR calc_energy IS: %llu \n", mult_count);

    #endif

    //save img
    free(partial);

    //unsigned char *energy_map = NULL;
    #ifdef debug
    char * fname = "energy_map.png";
    save_as_grayscale_image(fname, m - 2, n - 2, result);
    printf("Saved first energy map as %s\n", fname);
    // debug = 1;
    #endif
}

//first method called before anything else to create a 0 frame aronund original image
//make sure to free the returned pointer after the first seam is found and removed 
//repeat for each seam  
int * padd0_image(int n, int m, int * channels) {

    int size = 3 * (n + 2) * (m + 2);
    int * padded_image = (int * ) malloc(size * sizeof(int));

    //int padded_image[3][n+2][m+2];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < n + 2; j++) {
            for (int k = 0; k < m + 2; k++) {
                //if the column is 0 or m+1 or the row is 0 or n+1 we set 0 otherwise copy the value 
                if (j == 0 || k == 0 || j == n + 1 || k == m + 1) {
                    padded_image[i * (n + 2) * (m + 2) + (m + 2) * j + k] = 0;
                } else {
                    //printf("in pad its %f\n", channels[i*n*m + m*j + k]);
                    padded_image[i * (n + 2) * (m + 2) + (m + 2) * j + k] = channels[i * n * m + m * (j - 1) + k - 1];
                }
            }
        }
    }
    return padded_image;
}

// ========================== TESTING FUNCTIONS =================================

void test_computation() {
    int n = 4;
    int m = 4;
    int result[4];

    //prepadded with 0 to test the main convolution boi 
    int test_array[3 * 4 * 4] = {
        0,
        0,
        0,
        0,
        0,
        2,
        2,
        0,
        0,
        50,
        100,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        2,
        2,
        0,
        0,
        50,
        100,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        2,
        2,
        0,
        0,
        50,
        100,
        0,
        0,
        0,
        0,
        0
    };

    calc_RGB_energy(n, m, test_array, result); //n and m passed are the ones that are increased by 2 due to padding with 0 
    for (int i = 0; i < n - 2; i++) {
        for (int j = 0; j < m - 2; j++) {
            printf("%d\n", result[i * (m - 2) + j]);
        }
    }
}

void test_padding() {
    int n = 1;
    int m = 2;
    int channels[6] = {
        100,
        50,
        100,
        50,
        100,
        50
    };

    int * result = padd0_image(n, m, channels);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < n + 2; j++) {
            for (int k = 0; k < m + 2; k++) {
                printf("%d ", result[i * (n + 2) * (m + 2) + (m + 2) * j + k]);
            }
            printf("\n");
        }
        printf("\n\n");
    }
}

// ========================== TESTING FUNCTIONS =================================

// int main()
// {
//     test_computation();
//     //test_padding();

// }