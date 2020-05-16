#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse_img.h"
#include "count.h"

#define ABS(X) (((X)<0) ? (-(X)) : (X))
// #define debug   //uncomment for debugging

//--------------------  counter for instructions -------------------

#ifdef count_instr 
extern unsigned long long add_count; //count the total number of add instructions
extern unsigned long long mult_count;  //count the total number of mult instructions
#endif

//------------------------------------------------------------------

// int debug = 1;
//assuming that preprocessing is made of 0 padding 
// Given n rows, m columns of channel F of some image and the kernel H computes partial gradient corresponding to H given
//F is of size 3 x n x m
void calc_energy(int n, int m, int* F, int* part_grad ){
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    #ifdef count_instr        //counting adds and mults of this function
    unsigned long long count_ifs = 0;        //includes explicit ifs and for loop ifs  -> ADDS
    unsigned long long indexing = 0;         //includes increments of i.j,k variables  -> ADDS
    unsigned long long pointer_adds = 0;     //pointer arithmetic                      -> ADDS
    unsigned long long pointer_mults = 0;    //                                        -> MULTS
    #endif

    int *Fstep = F, *Fstep1 = F+m, *Fstep2 = F+2*m, step = m;
    for(int i = 1 ; i < n-1 ; i++){
        for(int j = 1 ; j < m-1 ; j++){
            int acc1;
            int acc2;
            int acc3;
            int acc4;
            int acc5;
            int acc6;
            int acc_total;

            int accF1 = *(Fstep + j - 1);
            int accF2 = *(Fstep + j);
            int accF3 = *(Fstep + j + 1);
            int accF4 = *(Fstep1 + j - 1);
            int accF5 = *(Fstep1 + j + 1);
            int accF6 = *(Fstep2 + j - 1);
            int accF7 = *(Fstep2 + j);
            int accF8 = *(Fstep2 + j + 1);

            //H_y
            acc1 = -(accF1 + (accF2 << 1));
            acc2 = accF6 - accF3;
            acc3 = (accF7 << 1) + accF8;
            acc_total = ABS(acc1 + acc2 + acc3);
            //H_x
            acc4 = accF3 - accF1;
            acc5 = (accF5 - accF4) << 1;
            acc6 = accF8 - accF6;
            acc_total += ABS(acc4 + acc5 + acc6);

            *(part_grad + step + j) = acc_total;

            #ifdef count_instr      //count line 55
            mult_count ++;
            pointer_adds += 2*2;    //assuming worse case that the after having the pointer arithmetic done its not saved but redone
            pointer_mults += 2;
            #endif
        }

        Fstep += m;
        Fstep1 += m;
        Fstep2 += m;
        step += m;

        #ifdef count_instr  //count line 54
        count_ifs ++;       //when not taken check must be made too
        pointer_adds += 2;  //need to dereference first to compare 
        pointer_mults ++; 
        #endif
    }

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
    unsigned long long count_ifs = 0;        //includes explicit ifs and for loop ifs  -> ADDS
    unsigned long long indexing = 0;         //includes increments of i.j,k variables  -> ADDS
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

    #ifdef count_instr        //counts lines 120-124
    count_ifs += 4;          
    indexing += 3;         
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

    #ifdef count_instr                                       //counts lines 134-138
    count_ifs += n-1 + (n-2)*(m-1) + (n-2)*(m-2)*4;          
    indexing += n-2 + (n-2)*(m-2) + (n-2)*(m-2)*3;         
    pointer_adds += (n-2)*(m-2)*3*(5 + 3 + 3);     
    pointer_mults += (n-2)*(m-2)*3*7;  
    add_count +=  (n-2)*(m-2)*3*2;                           //count directly the adds in line 154

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
    char *fname = "energy_map.png";
    save_as_grayscale_image(fname, m-2, n-2, result);
    printf("Saved first energy map as %s\n", fname);
    // debug = 1;
  #endif
}

//first method called before anything else to create a 0 frame aronund original image
//make sure to free the returned pointer after the first seam is found and removed 
//repeat for each seam  
int* padd0_image(int n, int m, int* channels){

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
          padded_image[i*(n+2)*(m+2) + (m+2)*j + k] = channels[i*n*m + m*(j-1) + k-1];
        }
      }
    }
  }
  return padded_image;
}


// ========================== TESTING FUNCTIONS =================================

void test_computation(){
    int n = 4;
    int m = 4; 
    int result[4]; 

    //prepadded with 0 to test the main convolution boi 
    int test_array[3*4*4] = {0,0,0,0,  0,2,2,0,  0,50,100,0,  0,0,0,0,         0,0,0,0,  0,2,2,0,  0,50,100,0,  0,0,0,0,          0,0,0,0,  0,2,2,0,  0,50,100,0,  0,0,0,0 };

    calc_RGB_energy(n,m,test_array,result); //n and m passed are the ones that are increased by 2 due to padding with 0 
    for(int i = 0 ; i < n-2 ; i++){ 
      for(int j = 0 ; j < m-2 ; j++){
        printf("%d\n", result[i*(m-2)+j]);
      }
    }
}

void test_padding(){
  int n =1;
  int m =2;
  int channels[6] = {100, 50 , 100, 50 ,100 ,50 };

  int* result = padd0_image(n,m,channels);
  for(int i = 0 ; i < 3 ; i++){
    for(int j = 0 ; j < n+2 ; j++){
      for(int k = 0 ; k < m+2 ; k++){
        printf("%d ", result[i*(n+2)*(m+2) + (m+2)*j + k]);
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
