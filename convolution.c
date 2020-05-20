#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse_img.h"
// #include "count.h"

#define K 1
#ifdef count_instr 
#define ABS(X)(((X) < 0) ? ((-(X)), count_ifs++) : (X))
#else 
#define ABS(X)(((X) < 0) ? (-(X)) : (X))
#endif
//#define debug   //uncomment for debugging

//--------------------  counter for instructions -------------------

#ifdef count_instr
extern unsigned long long add_count; //count the total number of add instructions
extern unsigned long long mult_count; //count the total number of mult instructions
#endif

//------------------------------------------------------------------

void calc_energy(int n, int m, int * F, int * part_grad) {
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    #ifdef count_instr //counting adds and mults of this function
    unsigned long long count_ifs = 0; //includes explicit ifs and for loop ifs -> ADDS
    unsigned long long indexing = 0; //includes increments of i.j,k variables  -> ADDS
    unsigned long long pointer_adds = 0; //pointer arithmetic                  -> ADDS
    unsigned long long pointer_mults = 0; //                                   -> MULTS
    unsigned long long shifts = 0;
    #endif

    int block_height_L1 = n - K;
    int block_height_R = n - K;
 
    int block_width_L1 = 1598;      //working set size is 4(m+2)+m
    int block_width_R = 2;          //working set size is 3(m+2)+m

    int width_limit_L1 = m - K - block_width_L1 + 1;
    int width_limit_R ; //m - K - block_width_R + 1

    int j, j_old;
    int ii_limit = n - K;

    #ifdef count_instr 
    add_count += 6;
    #endif

    for (j = 1; j < width_limit_L1; j = j + block_width_L1) {
        width_limit_R = j + block_width_L1 - block_width_R + 1; //take care in case block_width_L2 doesnt divide width_limit_L2
        int j_L3;

        for (int ii = 1; ii < ii_limit; ii++) { //single level 2 block calculation 
            for (j_L3 = j; j_L3 < width_limit_R; j_L3 += block_width_R) { //single level 3 block calculation
                int acc;
                int acc2;

                int nw = F[(ii - 1) * m + (j_L3 - 1)];  //3
                int n = F[(ii - 1) * m + j_L3]; //2
                int ne = F[(ii - 1) * m + j_L3 + 1]; //3
                int se = F[(ii + 1) * m + j_L3 + 1]; //3 
                int sw = F[(ii + 1) * m + j_L3 - 1]; //3
                int s = F[(ii + 1) * m + j_L3]; //2
                int e = F[ii * m + j_L3 + 1]; //2
                int w = F[ii * m + j_L3 - 1]; //2

                int nee = F[(ii - 1) * m + j_L3 + 2]; //LOL 3
                int see = F[(ii + 1) * m + j_L3 + 2]; //3
                int ee = F[ii * m + j_L3 + 2]; //2
                int center = F[ii * m + j_L3]; //1

                int dst = *(part_grad + ii * m + j_L3); //2
                int dst2 = *(part_grad + ii * m + j_L3 + 1);    //3

                //H_y
                dst = (-(nw + ((n) << 1))) + (sw - ne) + (((s) << 1) + se);
                dst = ABS(dst);
                //H_x
                acc = (ne - nw) + (((e - w) << 1)) + (se - sw);
                dst += ABS(acc);

                //H_y
                dst2 = (-(n + ((ne) << 1))) + (s - nee) + (((se) << 1) + see);
                dst2 = ABS(dst2);
                //H_x
                acc2 = (nee - n) + (((ee - center) << 1)) + (see - s);
                dst2 += ABS(acc2);
                
                *(part_grad + ii * m + j_L3) = dst;
                *(part_grad + ii * m + j_L3 + 1) = dst2;

                #ifdef count_instr 
                indexing += 0;
                pointer_adds += 0; 
                pointer_mults += 0;
                add_count += 
                mult_count += 
                shifts += 
                #endif

            }

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
                int w = F[ii * m + j_L3 - 1];

                int dst = *(part_grad + ii * m + j_L3);

                //H_y
                dst = (-(nw + ((n) << 1))) + (sw - ne) + (((s) << 1) + se);
                dst = ABS(dst);
                //H_x
                acc = (ne - nw) + (((e - w) << 1)) + (se - sw);
                dst += ABS(acc);

                *(part_grad + ii * m + j_L3) = dst;

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

    for (int ii = 1; ii < ii_limit; ii++) { //single level reg block calculation 
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
            int w = F[ii * m + j - 1];

            int nee = F[(ii - 1) * m + j + 2]; //LOL
            int see = F[(ii + 1) * m + j + 2];
            int ee = F[ii * m + j + 2];
            int center = F[ii * m + j];

            int dst = *(part_grad + ii * m + j);
            int dst2 = *(part_grad + ii * m + j + 1);

            //H_y
            dst = (-(nw + ((n) << 1))) + (sw - ne) + (((s) << 1) + se);
            dst = ABS(dst);
            //H_x
            acc = (ne - nw) + (((e - w) << 1)) + (se - sw);
            dst += ABS(acc);

            //H_y
            dst2 = (-(n + ((ne) << 1))) + (s - nee) + (((se) << 1) + see);
            dst2 = ABS(dst2);
            //H_x
            acc2 = (nee - n) + (((ee - center) << 1)) + (see - s);
            dst2 += ABS(acc2);
            
            *(part_grad + ii * m + j) = dst;
            *(part_grad + ii * m + j + 1) = dst2;
        }        
    }

    j_old = j;

    for (int ii = 1; ii < ii_limit; ii++) { //single level 2 block calculation 
         for (j = j_old; j < m - K; j++) {
                int acc;

                int nw = F[(ii - 1) * m + (j - 1)];
                int n = F[(ii - 1) * m + j];
                int ne = F[(ii - 1) * m + j + 1];
                int se = F[(ii + 1) * m + j + 1];
                int sw = F[(ii + 1) * m + (j - 1)];
                int s = F[(ii + 1) * m + j];
                int e = F[ii * m + j + 1];
                int w = F[ii * m + j - 1];

                int dst = *(part_grad + ii * m + j);

                //H_y
                dst = (-(nw + ((n) << 1))) + (sw - ne) + (((s) << 1) + se);
                dst = ABS(dst);
                //H_x
                acc = (ne - nw) + (((e - w) << 1)) + (se - sw);
                dst += ABS(acc);

                *(part_grad + ii * m + j) = dst;
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