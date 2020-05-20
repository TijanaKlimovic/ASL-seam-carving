#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse_img.h"
//#include "count.h"

#define K 1
#define ABS(X) (((X)<0) ? (-(X)) : (X))
// #define debug   //uncomment for debugging

//--------------------  counter for instructions -------------------

#ifdef count_instr 
extern unsigned long long add_count; //count the total number of add instructions
extern unsigned long long mult_count;  //count the total number of mult instructions
#endif

//------------------------------------------------------------------

//assumes channels is padded with 0 of size 3 x n x m
//assumes result is of size n-2 x m-2 
//calculates the cumulative sum over the individual channel energies
void calc_RGB_energy(int n, int m, int* padded, int* energy ){
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

            // channel R
            //H_y
            int k = 0;
            acc1 = -(padded[(i - 1) * m * 3  + (j - 1)*3 + k]
                + ((padded[(i - 1) * m *3 + j * 3 + k]) << 1));
            acc2 = padded[(i + 1) * m * 3 + (j - 1) * 3 + k]
                - padded[(i - 1) * m * 3 + (j + 1) * 3 + k];
            acc3 = ((padded[(i + 1) * m * 3 + j*3 + k]) << 1)
                + padded[(i + 1) * m * 3 + (j + 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (j-1)) = ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = padded[(i - 1) * m * 3 + (j + 1) * 3 + k]
                - padded[(i - 1) * m * 3 + (j - 1) * 3 + k];
            acc5 = (padded[i * m  * 3+ (j + 1) * 3 + k]
                - padded[i * m * 3 + (j - 1) * 3 + k]) << 1;
            acc6 = padded[(i + 1) * m * 3 + (j + 1) * 3 + k]
                - padded[(i + 1) * m * 3 + (j - 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (j-1)) += ABS(acc4 + acc5 + acc6);

            // channel G
            //H_y
            k = 1;
            acc1 = -(padded[(i - 1) * m * 3  + (j - 1)*3 + k]
                + ((padded[(i - 1) * m *3 + j * 3 + k]) << 1));
            acc2 = padded[(i + 1) * m * 3 + (j - 1) * 3 + k]
                - padded[(i - 1) * m * 3 + (j + 1) * 3 + k];
            acc3 = ((padded[(i + 1) * m * 3 + j*3 + k]) << 1)
                + padded[(i + 1) * m * 3 + (j + 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (j-1)) += ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = padded[(i - 1) * m * 3 + (j + 1) * 3 + k]
                - padded[(i - 1) * m * 3 + (j - 1) * 3 + k];
            acc5 = (padded[i * m  * 3+ (j + 1) * 3 + k]
                - padded[i * m * 3 + (j - 1) * 3 + k]) << 1;
            acc6 = padded[(i + 1) * m * 3 + (j + 1) * 3 + k]
                - padded[(i + 1) * m * 3 + (j - 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (j-1)) += ABS(acc4 + acc5 + acc6);

            // channel B
            //H_y
            k = 2;
            acc1 = -(padded[(i - 1) * m * 3  + (j - 1)*3 + k]
                + ((padded[(i - 1) * m *3 + j * 3 + k]) << 1));
            acc2 = padded[(i + 1) * m * 3 + (j - 1) * 3 + k]
                - padded[(i - 1) * m * 3 + (j + 1) * 3 + k];
            acc3 = ((padded[(i + 1) * m * 3 + j*3 + k]) << 1)
                + padded[(i + 1) * m * 3 + (j + 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (j-1)) += ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = padded[(i - 1) * m * 3 + (j + 1) * 3 + k]
                - padded[(i - 1) * m * 3 + (j - 1) * 3 + k];
            acc5 = (padded[i * m  * 3+ (j + 1) * 3 + k]
                - padded[i * m * 3 + (j - 1) * 3 + k]) << 1;
            acc6 = padded[(i + 1) * m * 3 + (j + 1) * 3 + k]
                - padded[(i + 1) * m * 3 + (j - 1) * 3 + k];
            *(energy + (i-1)*(m-2) + (j-1)) += ABS(acc4 + acc5 + acc6);

            #ifdef count_instr
            pointer_adds += 3*2; // i-1 and i+1
            pointer_adds += 3*2; // j-1 and j+1
            pointer_mults += 3*10; // all
            pointer_adds += 3*36; // all padded[...]
            add_count += 17; // enery value adds
            pointer_adds += 3*3; // energy pointer
            #endif
        }
    }

    #ifdef count_instr 
    //count total
    add_count += pointer_adds; 
    mult_count += pointer_mults;
    #endif
}


//first method called before anything else to create a 0 frame aronund original image
//make sure to free the returned pointer after the first seam is found and removed 
//repeat for each seam  
int* padd0_image(int n, int m, unsigned char* channels){

  int size = (n+2) * (m+2) * 3;
  int* padded_image = (int*) malloc( size*sizeof(int));

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
          padded_image[i*(m+2)*3 + j*3 + k] = channels[(i-1)*m*3 + (j-1)*3 + k];
        }
    }
  }
}
  return padded_image;
}

