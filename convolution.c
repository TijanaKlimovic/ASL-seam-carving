#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse_img.h"
// #include "count.h"

#define K 1

#ifdef count_instr 
#define ABS(X) (((X) < 0) ? ((-(X)), count_ifs++) : (X))
#else 
#define ABS(X) (((X) < 0) ? (-(X)) : (X))
#endif

// #define debug   //uncomment for debugging

//--------------------  counter for instructions -------------------

#ifdef count_instr
extern unsigned long long add_count; //count the total number of add instructions
extern unsigned long long mult_count; //count the total number of mult instructions
#endif

//------------------------------------------------------------------

//version with no accumulators only accounting for the 
void calc_energy(int n, int m, int * F, int * part_grad) {
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    #ifdef count_instr //counting adds and mults of this function
    unsigned long long count_ifs = 0; //includes explicit ifs and for loop ifs  -> ADDS
    unsigned long long indexing = 0; //includes increments of i.j,k variables   -> ADDS
    unsigned long long pointer_adds = 0; //pointer arithmetic                   -> ADDS
    unsigned long long pointer_mults = 0; //                                    -> MULTS
    unsigned long long shift_count = 0;
    #endif

    int block_height_L1 = n - K;
    int block_height_R = n - K;
    #ifdef count_instr
    add_count += 2;
    #endif

    int block_width_L1 = 1598;
    int block_width_R = 2;  //working set size is 4(m+2)+m

    int width_limit_L1 = m - K - block_width_L1 + 1;
    #ifdef count_instr
    add_count += 3;
    #endif

    int width_limit_R;

    int j_old;
    int ii_limit = n - K;
    #ifdef count_instr
    add_count++;
    #endif

    int j, i;
    for (j = 1; j < width_limit_L1; j = j + block_width_L1) {
        width_limit_R = j + block_width_L1 - block_width_R + 1; //take care in case block_width_L2 doesnt divide width_limit_L2
        int j_L3;

        for (int ii = 1; ii < ii_limit; ii++) { //single level 2 block calculation 
            int acc_save = F[(ii - 1) * m + j - 1] + (F[(ii * m) + j - 1] << 1) + F[(ii + 1) * m + j - 1];
            int acc_save2 = F[(ii - 1) * m + j] + (F[ii * m + j] << 1) + F[(ii + 1) * m + j];

            #ifdef count_instr
            pointer_adds += 19;
            pointer_mults += 6;
            add_count += 4;
            shift_count += 2;
            #endif

            for (j_L3 = j; j_L3 < width_limit_R; j_L3 += block_width_R) { //single level 3 block calculation
                int acc;
                int acc2;

                int nw = F[(ii - 1) * m + (j_L3 - 1)];
                int n = F[(ii - 1) * m + j_L3];
                int ne = F[(ii - 1) * m + j_L3 + 1];
                int se = F[(ii + 1) * m + j_L3 + 1];
                int sw = F[(ii + 1) * m + j_L3 - 1];
                int s = F[(ii + 1) * m + j_L3];
                int e = F[ii * m + j_L3 + 1];

                int nee = F[(ii - 1) * m + j_L3 + 2]; //LOL
                int see = F[(ii + 1) * m + j_L3 + 2];
                int ee = F[ii * m + j_L3 + 2];

                int dst = *(part_grad + ii * m + j_L3);
                int dst2 = *(part_grad + ii * m + j_L3 + 1);

                #ifdef count_instr
                pointer_adds += 41;
                pointer_mults += 12;
                #endif

                int new_acc_save = (ne + se) + (e << 1); //2a 1s
                int new_acc_save2 = (nee + see) + (ee << 1); //2a 1s
                //4a 2s
                //H_y
                dst = (sw - ne) + (((s) << 1) + se) - (nw + ((n) << 1)); //5a 2s
                dst = ABS(dst);
                //9a 4s
                //H_x
                acc = new_acc_save - acc_save; //1a
                dst += ABS(acc); //1a
                //11a 4s
                //H_y
                dst2 = (s - nee) + (((se) << 1) + see) - (n + ((ne) << 1)); //5a 2s
                dst2 = ABS(dst2);
                //16a 6s
                //H_x
                acc2 = new_acc_save2 - acc_save2; //1a
                dst2 += ABS(acc2); //1a
                //18a 6s
                #ifdef count_instr
                add_count += 18;
                shift_count += 6;
                #endif

                acc_save = new_acc_save;
                acc_save2 = new_acc_save2;
                
                *(part_grad + ii * m + j_L3) = dst;
                *(part_grad + ii * m + j_L3 + 1) = dst2;

                #ifdef count_instr
                pointer_adds += 5;
                pointer_mults += 2;
                #endif
            }
            #ifdef count_instr
            indexing += (width_limit_R - j) / block_width_R;
            count_ifs += (width_limit_R - j) / block_width_R + 1;
            #endif


            #ifdef count_instr
            indexing += j + block_width_L1 - j_L3;
            count_ifs += j + block_width_L1 - j_L3 + 1;
            #endif
            //bad non cooparating jjs
            for (; j_L3 < j + block_width_L1; j_L3++) {
                int acc;

                int nw = F[(ii - 1) * m + (j_L3 - 1)];
                int n = F[(ii - 1) * m + j_L3];
                int ne = F[(ii - 1) * m + j_L3 + 1];
                int se = F[(ii + 1) * m + j_L3 + 1];
                int sw = F[(ii + 1) * m + (j_L3 - 1)];
                int s = F[(ii + 1) * m + j_L3];
                int e = F[ii * m + j_L3 + 1];

                #ifdef count_instr
                pointer_adds += 25;
                pointer_mults += 7;
                #endif

                int dst = *(part_grad + ii * m + j_L3);

                #ifdef count_instr
                pointer_adds += 2;
                pointer_mults++;
                #endif

                int new_acc_save = (ne + se) + (e << 1);

                //H_y
                dst = (sw - ne) + (((s) << 1) + se) - (nw + ((n) << 1));
                dst = ABS(dst);
                //H_x
                acc = new_acc_save - acc_save;
                dst += ABS(acc);

                *(part_grad + ii * m + j_L3) = dst;
                #ifdef count_instr
                add_count += 8;
                shift_count += 3;
                pointer_adds += 2;
                pointer_mults++;
                #endif
            }
            
        }
        #ifdef count_instr 
        indexing += ii_limit - 1;
        count_ifs += ii_limit;
        #endif

    }
    #ifdef count_instr
    indexing += (width_limit_L1 - 1) / block_width_L1 + 1;
    count_ifs += (width_limit_L1 - 1) / block_width_L1 + 2;
    #endif

    //bad non cooparating js
    int j_limit =  m - K - block_width_R + 1;
    #ifdef count_instr
    add_count += 3;
    #endif
    j_old = j;


    for (int ii = 1; ii < ii_limit; ii++) { //single level reg block calculation 
        int acc_save = F[(ii - 1) * m + j_old - 1] + (F[(ii * m) + j_old - 1] << 1) + F[(ii + 1) * m + j_old - 1];
        int acc_save2 = F[(ii - 1) * m + j_old] + (F[ii * m + j_old] << 1) + F[(ii + 1) * m + j_old];
        #ifdef count_instr
        pointer_adds += 19;
        pointer_mults += 6;
        add_count += 4;
        shift_count += 2;
        #endif

        for (j = j_old; j < j_limit; j = j + block_width_R) {
            int acc;
            int acc2;

            int nw = F[(ii - 1) * m + (j - 1)];
            int n = F[(ii - 1) * m + j];
            int ne = F[(ii - 1) * m + j + 1];
            int se = F[(ii + 1) * m + j + 1];
            int sw = F[(ii + 1) * m + j - 1];
            int s = F[(ii + 1) * m + j];
            int e = F[ii * m + j + 1];

            int nee = F[(ii - 1) * m + j + 2]; //LOL
            int see = F[(ii + 1) * m + j + 2];
            int ee = F[ii * m + j + 2];

            int dst = *(part_grad + ii * m + j);
            int dst2 = *(part_grad + ii * m + j + 1);

            #ifdef count_instr
            pointer_adds += 41;
            pointer_mults += 12;
            #endif

            int new_acc_save = (ne + se) + (e << 1);
            int new_acc_save2 = (nee + see) + (ee << 1);

            //H_y
            dst = (sw - ne) + (((s) << 1) + se) - (nw + ((n) << 1));
            dst = ABS(dst);
            //H_x
            // acc = (ne - nw) + (((e - w) << 1)) + (se - sw);
            acc = new_acc_save - acc_save;
            dst += ABS(acc);

            //H_y
            dst2 = (s - nee) + (((se) << 1) + see) - (n + ((ne) << 1));
            dst2 = ABS(dst2);
            //H_x
            // acc2 = (nee - n) + (((ee - center) << 1)) + (see - s);
            acc2 = new_acc_save2 - acc_save2;
            dst2 += ABS(acc2);

            #ifdef count_instr
            add_count += 18;
            shift_count += 6;
            #endif

            acc_save = new_acc_save;
            acc_save2 = new_acc_save2;
            
            *(part_grad + ii * m + j) = dst;
            *(part_grad + ii * m + j + 1) = dst2;

            #ifdef count_instr
            pointer_adds += 5;
            pointer_mults += 2;
            #endif
        }
        #ifdef count_instr
        indexing += (j_limit - j_old) / block_width_R;
        count_ifs += (j_limit - j_old) / block_width_R + 1;
        #endif
    }
    #ifdef count_instr
    indexing += (ii_limit - 1);
    count_ifs += ii_limit;
    #endif

    j_old = j;
    j_limit = m - K;
    #ifdef count_instr
    add_count++;
    #endif

    for (int ii = 1; ii < ii_limit; ii++) { //single level 2 block calculation 
        for (j = j_old; j < j_limit; j++) {
            int acc;

            int nw = F[(ii - 1) * m + (j - 1)];
            int n = F[(ii - 1) * m + j];
            int ne = F[(ii - 1) * m + j + 1];
            int se = F[(ii + 1) * m + j + 1];
            int sw = F[(ii + 1) * m + (j - 1)];
            int s = F[(ii + 1) * m + j];
            int e = F[ii * m + j + 1];
            int w = F[ii * m + j - 1];

            #ifdef count_instr
            pointer_adds += 28;
            pointer_mults += 8;
            #endif

            int dst = *(part_grad + ii * m + j);

            #ifdef count_instr
            pointer_adds += 2;
            pointer_mults++;
            #endif

            //H_y
            dst = (sw - nw) + ((s - n) << 1) + (se - ne);
            // dst = (sw - ne) + (((s) << 1) + se) - (nw + ((n) << 1));
            dst = ABS(dst);
            //H_x
            acc = (ne - nw) + (((e - w) << 1)) + (se - sw);
            dst += ABS(acc);

            *(part_grad + ii * m + j) = dst;

            #ifdef count_instr
            add_count += 11;
            shift_count += 2;
            pointer_adds += 2;
            pointer_mults++;
            #endif
        }
        #ifdef count_instr
        indexing += j_limit - j_old;
        count_ifs += j_limit - j_old + 1;
        #endif    
    }
    #ifdef count_instr
    indexing += ii_limit - 1;
    count_ifs += ii_limit;
    #endif

    #ifdef count_instr
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
