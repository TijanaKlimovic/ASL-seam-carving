#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse_img.h"

int debug = 0;
//assuming that preprocessing is made of 0 padding 
// Given n rows, m columns of channel F of some image and the kernel H computes partial gradient corresponding to H given

void calc_energy(int n, int m, int k, double* F , double* part_grad , double H[3][3] ){
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            *(part_grad + i*m + j) = 0;
        }
    }

    for(int i = 1 ; i < n-k ; i++){
        for(int j = 1 ; j < m-k ; j++){
            for(int u = -k ; u <= k; u++){
                for(int v = -k ; v <= k ; v++){
                    *(part_grad + i*m + j) += H[u+k][v+k]*F[(i+u)*m + j+v];
                }
            }
          //calculate absolute value of each element in partial derivative of channel F 
          *(part_grad + i*m + j) = fabs(*(part_grad + i*m + j));
        }
    }
}

//assumes channels is padded with 0 of size 3 x n x m
//assumes result is of size n-2 x m-2 
//calculates the cumulative sum over the individual channel energies
//returns the energy map result over color image of size  n-2 x m-2 
void calc_RGB_energy(int n, int m, double* channels, double* result){
    //fixed kernels 
  double H_y[3][3] = {
    {-1,-2,-1},
    {0,0,0},
    {1,2,1}};

  double H_x[3][3] = {
    {-1,0,1},
    {-2,0,2},
    {-1,0,1}};
    
    int k= 1;

    //double partial_x[3][n][m]; //3d parital derivative in x
    //double partial_y[3][n][m]; //3d partial derivative in y 

    int size = 3*n*m ;

    double* partial_x = (double*) malloc( size*sizeof(double));
    double* partial_y = (double*) malloc( size*sizeof(double));
    /*
    if(partial_x == NULL){
      printf("%s\n", "too much FOR X !");
    }
    */
    /*
    if(partial_y == NULL){
      printf("%s\n", "too much FOR Y !");
    }*/

    //calculate the parital derivatives 
    for(int i = 0 ; i < 3 ; i ++){
      //pass the ith channel for energy calculation
       calc_energy(n,m,k,channels + n*m*i, partial_x + n*m*i, H_x);
       calc_energy(n,m,k,channels + n*m*i, partial_y + n*m*i, H_y);
    }

    for (int i = 0; i < n - 2; i++) {
        for (int j = 0; j < m - 2; j++) {
            result[(m - 2) * i + j] = 0;
        }
    }

    //calculate the total 3d energy 
    
      for(int j = 1 ; j < n-1 ; j ++){
        for(int k = 1 ; k < m-1 ; k++){
          for(int i = 0 ; i < 3 ; i ++){
            //add elementwise along the z axis 
          *(result+(m-2)*(j-1)+k-1) += *(partial_x + i*m*n + j*m + k) + *(partial_y + i*m*n + j*m + k);
        }
      } 
    }
  // save img
    free(partial_x);
    free(partial_y);

  //unsigned char *energy_map = NULL;
  if (!debug) {
    char *fname = "energy_map.png";
    save_as_grayscale_image(fname, m-2, n-2, result);
    printf("Saved first energy map as %s\n", fname);
    debug = 1;
  }
}



//first method called before anything else to create a 0 frame aronund original image
//make sure to free the returned pointer after the first seam is found and removed 
//repeat for each seam  
double* padd0_image(int n, int m, double* channels){

  int size = 3*(n+2)*(m+2);
  double* padded_image = (double*) malloc( size*sizeof(double));

  //double padded_image[3][n+2][m+2];
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
    double result[4]; 

    //prepadded with 0 to test the main convolution boi 
    double test_array[3*4*4] = {0,0,0,0,  0,2,2,0,  0,50,100,0,  0,0,0,0,         0,0,0,0,  0,2,2,0,  0,50,100,0,  0,0,0,0,          0,0,0,0,  0,2,2,0,  0,50,100,0,  0,0,0,0 };

    calc_RGB_energy(n,m,test_array,result); //n and m passed are the ones that are increased by 2 due to padding with 0 
    for(int i = 0 ; i < n-2 ; i++){ 
      for(int j = 0 ; j < m-2 ; j++){
        printf("%f\n", result[i*(m-2)+j]);
      }
    }
}

void test_padding(){
  int n =1;
  int m =2;
  double channels[6] = {100, 50 , 100, 50 ,100 ,50 };

  double* result = padd0_image(n,m,channels);
  for(int i = 0 ; i < 3 ; i++){
    for(int j = 0 ; j < n+2 ; j++){
      for(int k = 0 ; k < m+2 ; k++){
        printf("%f ", result[i*(n+2)*(m+2) + (m+2)*j + k]);
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
