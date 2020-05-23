#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse_img.h"

#define K 1
#ifdef count_instr 
#define ABS(X)(((X) < 0) ? (mult_count++, (-(X))) : (X))
#else 
#define ABS(X)(((X) < 0) ? (-(X)) : (X))
#endif

//#define debug   //uncomment for debugging

//--------------------  counter for instructions -------------------

#ifdef count_instr 
extern unsigned long long add_count; //count the total number of add instructions
extern unsigned long long mult_count;  //count the total number of mult instructions
#endif

//------------------------------------------------------------------

// int debug = 1;
// assuming that preprocessing is made of 0 padding 
// Given n rows, m columns of channel F of some image and the kernel H computes partial gradient corresponding to H given
// F is of size 3 x n x m
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
    add_count += pointer_adds;
    mult_count += pointer_mults;
    printf("NO ADDS FOR calc_energy IS: %llu \n", add_count);
    printf("NO MULTS FOR calc_energy IS: %llu \n", mult_count);
    #endif
}


//assumes channels is padded with 0 of size 3 x n x m
//assumes result is of size n-2 x m-2 
//calculates the cumulative sum over the individual channel energies
//returns the energy map result over color image of size  n-2 x m-2 
void calc_RGB_energy(int n, int m, int* channels, int* result){
    
    #ifdef count_instr        //counting adds and mults of this function
    unsigned long long pointer_adds = 0;     //pointer arithmetic                      -> ADDS
    unsigned long long pointer_mults = 0;    //                                        -> MULTS
    #endif

    int* partial = (int*) malloc(3*n*m*sizeof(int));

    //calculate the parital derivatives 
    //pass the ith channel for energy calculation
    int step1 = n*m, step2 = 2*n*m;
    calc_energy(n,m,channels, partial);
    calc_energy(n,m,channels + step1, partial + step1);
    calc_energy(n,m,channels + step2, partial + step2);

    #ifdef count_instr //counts 89-96        
    pointer_adds += 3*4;     
    pointer_mults += 3*8;
    #endif

    //calculate the total 3d energy 
    int k;
    int *result_step = result-1, *partial_step = partial;
    int *result_step1 = result-1, *partial_step1 = partial+step1;
    int *result_step2 = result-1, *partial_step2 = partial+step2;
    // first channel
    for(int j = 1 ; j < n-1 ; j++) {
        partial_step += m;
        for(k = 1 ; k < m-16 ; k += 16) {
            //add elementwise along the z axis 
            *(result_step+k) = *(partial_step+k);
            *(result_step+k+1) = *(partial_step+k+1);
            *(result_step+k+2) = *(partial_step+k+2);
            *(result_step+k+3) = *(partial_step+k+3);
            *(result_step+k+4) = *(partial_step+k+4);
            *(result_step+k+5) = *(partial_step+k+5);
            *(result_step+k+6) = *(partial_step+k+6);
            *(result_step+k+7) = *(partial_step+k+7);
            *(result_step+k+8) = *(partial_step+k+8);
            *(result_step+k+9) = *(partial_step+k+9);
            *(result_step+k+10) = *(partial_step+k+10);
            *(result_step+k+11) = *(partial_step+k+11);
            *(result_step+k+12) = *(partial_step+k+12);
            *(result_step+k+13) = *(partial_step+k+13);
            *(result_step+k+14) = *(partial_step+k+14);
            *(result_step+k+15) = *(partial_step+k+15);
        }
        while(k < m-1) {
            *(result_step+k) = *(partial_step+k);
            k++;
        }
        result_step += m-2;
    }
    // second channel
    for(int j = 1 ; j < n-1 ; j++) {
        partial_step1 += m;
        for(k = 1 ; k < m-16 ; k += 16) {
            //add elementwise along the z axis 
            *(result_step1+k) += *(partial_step1+k);
            *(result_step1+k+1) += *(partial_step1+k+1);
            *(result_step1+k+2) += *(partial_step1+k+2);
            *(result_step1+k+3) += *(partial_step1+k+3);
            *(result_step1+k+4) += *(partial_step1+k+4);
            *(result_step1+k+5) += *(partial_step1+k+5);
            *(result_step1+k+6) += *(partial_step1+k+6);
            *(result_step1+k+7) += *(partial_step1+k+7);
            *(result_step1+k+8) += *(partial_step1+k+8);
            *(result_step1+k+9) += *(partial_step1+k+9);
            *(result_step1+k+10) += *(partial_step1+k+10);
            *(result_step1+k+11) += *(partial_step1+k+11);
            *(result_step1+k+12) += *(partial_step1+k+12);
            *(result_step1+k+13) += *(partial_step1+k+13);
            *(result_step1+k+14) += *(partial_step1+k+14);
            *(result_step1+k+15) += *(partial_step1+k+15);
        }
        while(k < m-1) {
            *(result_step1+k) += *(partial_step1+k);
            k++;
        }
        result_step1 += m-2;
    }
    // third channel
    for(int j = 1 ; j < n-1 ; j++) {
        partial_step2 += m;
        for(k = 1 ; k < m-16 ; k += 16) {
            //add elementwise along the z axis 
            *(result_step2+k) += *(partial_step2+k);
            *(result_step2+k+1) += *(partial_step2+k+1);
            *(result_step2+k+2) += *(partial_step2+k+2);
            *(result_step2+k+3) += *(partial_step2+k+3);
            *(result_step2+k+4) += *(partial_step2+k+4);
            *(result_step2+k+5) += *(partial_step2+k+5);
            *(result_step2+k+6) += *(partial_step2+k+6);
            *(result_step2+k+7) += *(partial_step2+k+7);
            *(result_step2+k+8) += *(partial_step2+k+8);
            *(result_step2+k+9) += *(partial_step2+k+9);
            *(result_step2+k+10) += *(partial_step2+k+10);
            *(result_step2+k+11) += *(partial_step2+k+11);
            *(result_step2+k+12) += *(partial_step2+k+12);
            *(result_step2+k+13) += *(partial_step2+k+13);
            *(result_step2+k+14) += *(partial_step2+k+14);
            *(result_step2+k+15) += *(partial_step2+k+15);
        }
        while(k < m-1) {
            *(result_step2+k) += *(partial_step2+k);
            k++;
        }
        result_step2 += m-2;
    }

    #ifdef count_instr // count 216
    add_count += 2*3*(n-2)*(m-2);
    pointer_adds += 8*3*(n-2)*(m-2);
    pointer_mults += 4*3*(n-2)*(m-2);
    #endif

    #ifdef count_instr
    //count total
    add_count += pointer_adds; 
    mult_count += pointer_mults;
    #endif

    //save img
    free(partial);

    //unsigned char *energy_map = NULL;
    #ifdef debug 
      char *fname = "energy_map.png";
      save_as_grayscale_image(fname, m-2, n-2, result);
      printf("Saved first energy map as %s\n", fname);
      // debug = 1;
    #endif
}



//first method called before anything else to create a 0 frame aronund original image
//make sure to free the returned pointer after the first seam is found and removed 
//repeat for each seam  
int* padd0_image(int n, int m, unsigned char* channels){

  int size = 3*(n+2)*(m+2);
  int* padded_image = (int*) malloc( size*sizeof(int));

  //int padded_image[3][n+2][m+2];
  for(int i = 0 ; i < 3 ; i++){
    for(int j = 0 ; j < n+2 ; j++){
      for(int k = 0 ; k < m+2 ; k++){
        //if the column is 0 or m+1 or the row is 0 or n+1 we set 0 otherwise copy the value 
        if(j == 0 || k == 0 || j == n+1 || k == m+1){
          padded_image[i*(n+2)*(m+2) + (m+2)*j + k] = 0;
        }else{
          //printf("in pad its %f\n", channels[i*n*m + m*j + k]);
          padded_image[i*(n+2)*(m+2) + (m+2)*j + k] = (int)channels[i*n*m + m*(j-1) + k-1];
        }
      }
    }
  }
  return padded_image;
}

