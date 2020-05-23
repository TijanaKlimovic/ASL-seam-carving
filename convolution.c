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
    unsigned long long count_ifs = 0; //includes explicit ifs and for loop ifs -> ADDS
    unsigned long long indexing = 0; //includes increments of i.j,k variables  -> ADDS
    unsigned long long pointer_adds = 0; //pointer arithmetic                  -> ADDS
    unsigned long long pointer_mults = 0; //                                   -> MULTS
    unsigned long long shifts = 0;
    #endif
 
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

        #ifdef count_instr 
        add_count += 3; //width_limit_R
        #endif

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

                //total 20 adds

                int nee = F[(ii - 1) * m + j_L3 + 2]; //LOL 3
                int see = F[(ii + 1) * m + j_L3 + 2]; //3
                int ee = F[ii * m + j_L3 + 2]; //2
                int center = F[ii * m + j_L3]; //1

                int dst = *(part_grad + ii * m + j_L3); //2
                int dst2 = *(part_grad + ii * m + j_L3 + 1);    //3

                //H_y
                dst = (-(nw + ((n) << 1))) + (sw - ne) + (((s) << 1) + se); //5
                dst = ABS(dst);
                //H_x
                acc = (ne - nw) + (((e - w) << 1)) + (se - sw); //5
                dst += ABS(acc); //1

                //H_y
                dst2 = (-(n + ((ne) << 1))) + (s - nee) + (((se) << 1) + see); //5
                dst2 = ABS(dst2);
                //H_x
                acc2 = (nee - n) + (((ee - center) << 1)) + (see - s); //5
                dst2 += ABS(acc2); //1
                
                *(part_grad + ii * m + j_L3) = dst; //2
                *(part_grad + ii * m + j_L3 + 1) = dst2; //3

                #ifdef count_instr 
                count_ifs += 5; //4 from abs and 1 for comparing in for loop
                indexing += 1;
                pointer_adds += 39; 
                pointer_mults += 16;
                add_count += 22;
                mult_count += 2; //-1 multiplications 
                shifts += 6;
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

                #ifdef count_instr 
                count_ifs += 3; //abs and for loop condition
                indexing += 1;  //jl INCREMENTS
                pointer_adds += 24; 
                pointer_mults += 10;
                add_count += 11;
                mult_count += 1;
                shifts += 3;
                #endif

            }

        #ifdef count_instr 
        count_ifs+= 2; //when not taken check must be made too and for loop condition for ii
        indexing += 1;  //ii increments 
        #endif
        }
        #ifdef count_instr 
        count_ifs+= 2; //when not taken check must be made too and for loop condition for j
        indexing += 1;  //j increments 
        #endif
    }

    //bad non cooparating js
    int j_limit =  m - K - block_width_R + 1;
    j_old = j;

    #ifdef count_instr 
    count_ifs++; //js last for loop comparison 
    add_count += 3; //j_limit
    #endif

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

            //20 pointer adds and 8 mults

            int nee = F[(ii - 1) * m + j + 2]; //LOL
            int see = F[(ii + 1) * m + j + 2];
            int ee = F[ii * m + j + 2];
            int center = F[ii * m + j];

            int dst = *(part_grad + ii * m + j);
            int dst2 = *(part_grad + ii * m + j + 1);

            //14 pointer adds and 6 mults

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

            #ifdef count_instr 
            count_ifs += 5; //4 from abs and 1 for comparing in for loop
            indexing += 1;
            pointer_adds += 39; 
            pointer_mults += 16;
            add_count += 22;
            mult_count += 2; //-1 multiplications 
            shifts += 6;
            #endif 
        }        
        #ifdef count_instr 
        count_ifs += 2; //one for j that isnt taken and 1 for ii for condition 
        indexing += 1; //ii increments
        #endif 
    }

    #ifdef count_instr 
    count_ifs ++; //one for ii not taken
    #endif 

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

            #ifdef count_instr 
            count_ifs += 3; //abs and for loop condition
            indexing += 1;  //j INCREMENTS
            pointer_adds += 24; 
            pointer_mults += 10;
            add_count += 11;
            mult_count += 1;
            shifts += 3;
            #endif
        }     
        #ifdef count_instr 
        count_ifs+= 2; //when not taken check must be made too and for loop condition for ii
        indexing += 1;  //ii increments 
        #endif   
    }

    #ifdef count_instr
    //count total
    count_ifs++;
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

