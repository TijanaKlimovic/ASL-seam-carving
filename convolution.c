#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse_img.h"

//--------------------  counter for instructions -------------------

#ifdef count_instr 
extern unsigned long long add_count; //count the total number of add instructions
extern unsigned long long mult_count;  //count the total number of mult instructions
#endif

//------------------------------------------------------------------

// Assuming that preprocessing is made of 0 padding 
// Given n rows, m columns of channel F of some image and the kernel H computes partial gradient corresponding to H given

void calc_energy(int n, int m, int* F, int* part_grad, int H[9] ){
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    int res, f_index_0, f_index_1, f_index_2;
    for(int i = 1 ; i < n-1 ; i++){
        f_index_0 = (i-1)*m;
        f_index_1 = i*m;
        f_index_2 = (i+1)*m;

        #ifdef count_instr
        add_count += 2;
        mult_count += 3;
        #endif

        for(int j = 1 ; j < m-1 ; j++){
            res = H[0]*F[f_index_0 + j-1]
                  + H[1]*F[f_index_0 + j]
                  + H[2]*F[f_index_0 + j+1]
                  + H[3]*F[f_index_1 + j-1]
                  + H[4]*F[f_index_1 + j]
                  + H[5]*F[f_index_1 + j+1]
                  + H[6]*F[f_index_2 + j-1]
                  + H[7]*F[f_index_2 + j]
                  + H[8]*F[f_index_2 + j+1];

            //calculate absolute value of each element in partial derivative of channel F 
            if(res < 0){
              res = (-1) * res;

              #ifdef count_instr
              mult_count++;
              #endif
            }

            *(part_grad + i*m + j) = res;

            #ifdef count_instr
            // pointer count
            add_count += 17;
            mult_count += 1;

            // operation count
            add_count += 8;
            mult_count += 9;
            #endif
        }
    }

    #ifdef count_instr
    // index count
    add_count += (n-2)*(m-2);
    #endif
}

//assumes channels is padded with 0 of size 3 x n x m
//assumes result is of size n-2 x m-2 
//calculates the cumulative sum over the individual channel energies
//returns the energy map result over color image of size  n-2 x m-2 
void calc_RGB_energy(int n, int m, int* channels, int* result){
  //fixed kernels 
  int H_y[9] = { -1, -2, -1, 0, 0, 0, 1, 2, 1};
  int H_x[9] = { -1, 0, 1, -2, 0, 2, -1, 0, 1};

  int size = 3*n*m ;
  int* partial_x = (int*) malloc(size*sizeof(int));
  int* partial_y = (int*) malloc(size*sizeof(int));

  int step_1 = n*m, step_2 = 2*n*m;

  calc_energy(n,m,channels, partial_x, H_x);
  calc_energy(n,m,channels, partial_y, H_y);
  calc_energy(n,m,channels + step_1, partial_x + step_1, H_x);
  calc_energy(n,m,channels + step_1, partial_y + step_1, H_y);
  calc_energy(n,m,channels + step_2, partial_x + step_2, H_x);
  calc_energy(n,m,channels + step_2, partial_y + step_2, H_y);

  #ifdef count_instr
  add_count += 8;
  mult_count += 5;
  #endif

  //calculate the total 3d energy 
  
  int *partial_x_index, *partial_y_index;
  for(int j = 1 ; j < n-1 ; j ++){
    for(int k = 1 ; k < m-1 ; k++){
      partial_x_index = partial_x + j*m + k;
      partial_y_index = partial_y + j*m + k;
      *(result+(m-2)*(j-1)+k-1) = *(partial_x_index) + *(partial_y_index)
                                + *(partial_x_index + step_1) + *(partial_y_index + step_1)
                                + *(partial_x_index + step_2) + *(partial_y_index + step_2);

      #ifdef count_instr
      // pointer count
      add_count += 9;
      mult_count++;

      // operation count
      add_count += 9;
      mult_count += 2;
      #endif
    } 
  }

  #ifdef count_instr
  // index count
  add_count += (n-2)*(m-2);
  printf("NO ADDS FOR calc_energy IS: %llu ", add_count); 
  printf("NO MULTS FOR calc_energy IS: %llu ", mult_count);
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
