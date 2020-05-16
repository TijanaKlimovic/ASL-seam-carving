#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse_img.h"

#define K 1
//#define debug   //uncomment for debugging

//--------------------  counter for instructions -------------------

#ifdef count_instr 
extern unsigned long long add_count; //count the total number of add instructions
extern unsigned long long mult_count;  //count the total number of mult instructions
#endif

//------------------------------------------------------------------

int debug = 1;
//assuming that preprocessing is made of 0 padding 
// Given n rows, m columns of channel F of some image and the kernel H computes partial gradient corresponding to H given

void calc_energy(int n, int m, int* F, int* part_grad, int H[3][3] ){
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    #ifdef count_instr        //counting adds and mults of this function
    unsigned long long count_ifs = 0;        //includes explicit ifs and for loop ifs  -> ADDS
    unsigned long long indexing = 0;         //includes increments of i.j,k variables  -> ADDS
    unsigned long long pointer_adds = 0;     //pointer arithmetic                      -> ADDS
    unsigned long long pointer_mults = 0;    //                                        -> MULTS
    #endif


    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            *(part_grad + i*m + j) = 0;
        }
    }

    #ifdef count_instr            //count lines 33-37
    count_ifs += n+1 + n*(m+1);
    indexing += n + n*m;
    pointer_mults += m*n;        //assuming perfect prediciton
    pointer_adds += 2*m*n;       //assuming perfect prediciton
    #endif

    for(int i = 1 ; i < n-K ; i++){
        for(int j = 1 ; j < m-K ; j++){
            for(int u = -K ; u <= K; u++){
                for(int v = -K ; v <= K ; v++){
                    *(part_grad + i*m + j) += H[u+K][v+K]*F[(i+u)*m + j+v];
                }
            }

            //calculate absolute value of each element in partial derivative of channel F 
            if(*(part_grad + i*m + j) < 0){
              *(part_grad + i*m + j) = (-1) * (*(part_grad + i*m + j));

              #ifdef count_instr      //count line 55
              mult_count ++;
              pointer_adds += 2*2;    //assuming worse case that the after having the pointer arithmetic done its not saved but redone
              pointer_mults += 2;
              #endif
            }

            #ifdef count_instr  //count line 54
            count_ifs ++;       //when not taken check must be made too
            pointer_adds += 2;  //need to dereference first to compare 
            pointer_mults ++; 
            #endif
        }
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
    //printf("NO ADDS FOR calc_energy IS: %llu \n", add_count); 
    //printf("NO MULTS FOR calc_energy IS: %llu \n", mult_count); 
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

  //fixed kernels 
  int H_y[3][3] = {
    {-1,-2,-1},
    {0,0,0},
    {1,2,1}};

  int H_x[3][3] = {
    {-1,0,1},
    {-2,0,2},
    {-1,0,1}};

    int size = 3*n*m ;

    int* partial_x = (int*) malloc( size*sizeof(int));
    int* partial_y = (int*) malloc( size*sizeof(int));

    //calculate the parital derivatives 
    for(int i = 0 ; i < 3 ; i ++){
      //pass the ith channel for energy calculation
       calc_energy(n,m,channels + n*m*i, partial_x + n*m*i, H_x);
       calc_energy(n,m,channels + n*m*i, partial_y + n*m*i, H_y);
    }


    #ifdef count_instr        //counts lines 120-124
    count_ifs += 4;          
    indexing += 3;         
    pointer_adds += 3*4;     
    pointer_mults += 3*8;    
    #endif

    for (int i = 0; i < n - 2; i++) {
        for (int j = 0; j < m - 2; j++) {
            result[(m - 2) * i + j] = 0;
        }
    }


    #ifdef count_instr                                      //counts lines 134-138
    count_ifs += n-1 + (n-2)*(m-1);          
    indexing += n-2 + (n-2)*(m-2);         
    pointer_adds += (n-2)*(m-2)*2;     
    pointer_mults += (n-2)*(m-2);    
    #endif

    //calculate the total 3d energy 
    
      for(int j = 1 ; j < n-1 ; j ++){
        for(int k = 1 ; k < m-1 ; k++){
          for(int i = 0 ; i < 3 ; i ++){
            //add elementwise along the z axis 
          *(result+(m-2)*(j-1)+k-1) += *(partial_x + i*m*n + j*m + k) + *(partial_y + i*m*n + j*m + k);
        }
      } 
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
    //printf("NO ADDS FOR calc_energy IS: %llu \n", add_count); 
    //printf("NO MULTS FOR calc_energy IS: %llu \n", mult_count); 

    #endif

    //save img
    free(partial_x);
    free(partial_y);

  //unsigned char *energy_map = NULL;
  #ifdef debug 
    char *fname = "energy_map.png";
    save_as_grayscale_image(fname, m-2, n-2, result);
    printf("Saved first energy map as %s\n", fname);
    debug = 1;
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
