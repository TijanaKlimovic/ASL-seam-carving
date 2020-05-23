#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse_img.h"

#define K 1
#define ABS(X) (((X)<0) ? (-(X)) : (X))
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
void calc_energy(int n, int m, int* F, int* part_grad){
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    #ifdef count_instr        //counting adds and mults of this function
    unsigned long long pointer_adds = 0;     //pointer arithmetic                      -> ADDS
    unsigned long long pointer_mults = 0;    //                                        -> MULTS
    #endif

    for(int i = 1 ; i < n-K ; i++){
        for(int j = 1 ; j < m-K ; j++){
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
            *(part_grad + i*m + j) = ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = F[(i - 1) * m + j + 1] - F[(i - 1) * m + (j - 1)];
            acc5 = (F[i * m + j + 1] - F[i * m + j - 1]) << 1;
            acc6 = F[(i + 1) * m + j + 1] - F[(i + 1) * m + j - 1];
            *(part_grad + i*m + j) += ABS(acc4 + acc5 + acc6);

            #ifdef count_instr //count 41-49
            pointer_adds += 7; pointer_mults += 2; // 41
            pointer_adds += 8; pointer_mults += 2; // 42
            pointer_adds += 7; pointer_mults += 2; // 43
            pointer_adds += 8; pointer_mults += 2; // 46
            pointer_adds += 6; pointer_mults += 2; // 47
            pointer_adds += 8; pointer_mults += 2; // 48

            // values adds
            add_count += 13;
            #endif
        }
    }

    #ifdef count_instr 
    //count total
    add_count += pointer_adds; 
    mult_count += pointer_mults;
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

    #ifdef count_instr // counts 83
    mult_count += 3;
    #endif

    //calculate the parital derivatives 
    for(int i = 0 ; i < 3 ; i ++){
        //pass the ith channel for energy calculation
        calc_energy(n,m,channels + n*m*i, partial + n*m*i);
    }

    #ifdef count_instr //counts 90-93
    pointer_adds += 3*2;
    pointer_mults += 3*2; // n*m*i counted once for each iteration
    #endif

    for (int i = 0; i < n - 2; i++) {
        for (int j = 0; j < m - 2; j++) {
            result[(m - 2) * i + j] = 0;
        }
    }

    //calculate the total 3d energy
    for(int i = 0 ; i < 3 ; i++) {
        for(int j = 1 ; j < n-1 ; j++) {
            for(int k = 1 ; k < m-1 ; k++) {
                //add elementwise along the z axis 
                *(result+(m-2)*(j-1)+k-1) += *(partial + i*m*n + j*m + k);

                #ifdef count_instr // count 111
                add_count += 2;
                pointer_adds += 8;
                pointer_mults += 4;
                #endif
            }
        }
    }

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

