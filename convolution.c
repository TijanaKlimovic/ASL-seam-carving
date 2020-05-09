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

//assuming that preprocessing is made of 0 padding 
// Given n rows, m columns of channel F of some image and the kernel H computes partial gradient corresponding to H given

void calc_energy(int n, int m, int* F, int* part_grad, int H[3][3] ){
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel


    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            *(part_grad + i*m + j) = 0;
        }
    }

    #ifdef count_instr
    // index count
    add_count += (n-1)*(m-1);

    // pointer count
    add_count += 2*(n-1)*(m-1);
    mult_count += (n-1)*(m-1);
    #endif

    for(int i = 1 ; i < n-K ; i++){
        for(int j = 1 ; j < m-K ; j++){
            for(int u = -K ; u <= K; u++){
                for(int v = -K ; v <= K ; v++){
                    *(part_grad + i*m + j) += H[u+K][v+K]*F[(i+u)*m + j+v];

                    #ifdef count_instr
                    // operation count
                    add_count++;
                    mult_count++;

                    // pointer count
                    add_count += 7;
                    mult_count += 2;
                    #endif
                }
            }

            #ifdef count_instr
            // index count (u and v)
            add_count += 9;
            #endif

            //calculate absolute value of each element in partial derivative of channel F 
            if(*(part_grad + i*m + j) < 0){
              *(part_grad + i*m + j) = (-1) * (*(part_grad + i*m + j));

              #ifdef count_instr
              // operation count
              mult_count ++;

              // pointer count
              add_count += 4;
              mult_count += 2;
              #endif
            }

            #ifdef count_instr
            // if count
            add_count += 2;
            mult_count++;
            #endif
        }
    }

    #ifdef count_instr
    // index count (i and j)
    add_count += (n-2)*(m-2);
    #endif
}

//assumes channels is padded with 0 of size 3 x n x m
//assumes result is of size n-2 x m-2 
//calculates the cumulative sum over the individual channel energies
//returns the energy map result over color image of size  n-2 x m-2 
void calc_RGB_energy(int n, int m, int* channels, int* result){
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

  #ifdef count_instr
  mult_count += 4;
  #endif

  //calculate the parital derivatives
  for(int i = 0 ; i < 3 ; i ++){
    //pass the ith channel for energy calculation
     calc_energy(n,m,channels + n*m*i, partial_x + n*m*i, H_x);
     calc_energy(n,m,channels + n*m*i, partial_y + n*m*i, H_y);
  }


  #ifdef count_instr
  // index count
  add_count += 3;

  // pointer count
  add_count += 12;
  mult_count += 24;
  #endif

  for (int i = 0; i < n - 2; i++) {
      for (int j = 0; j < m - 2; j++) {
          result[(m - 2) * i + j] = 0;
      }
  }

  #ifdef count_instr
  // index count
  add_count += (n-2)*(m-2);

  // pointer count
  add_count += 2*(n-2)*(m-2);
  mult_count += (n-2)*(m-2);    
  #endif

  //calculate the total 3d energy 
  
  for(int j = 1 ; j < n-1 ; j ++){
    for(int k = 1 ; k < m-1 ; k++){
      for(int i = 0 ; i < 3 ; i ++){
        //add elementwise along the z axis 
        *(result+(m-2)*(j-1)+k-1) += *(partial_x + i*m*n + j*m + k)
                                  + *(partial_y + i*m*n + j*m + k);

        #ifdef count_instr
        // pointer count
        add_count += 11;
        mult_count += 7;

        // operation count
        add_count += 2;
        #endif
      }
    } 
  }

  #ifdef count_instr
  // index count
  add_count += 3*(n-2)*(m-2);

  printf("%llu %llu\n", add_count, mult_count); 
  #endif

  //save img
  free(partial_x);
  free(partial_y);
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
