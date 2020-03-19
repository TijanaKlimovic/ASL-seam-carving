#include <stdio.h>
#include <stdlib.h>

//assuming that preprocessing is made of 0 padding 
// Given n rows, m columns of channel F of some image and the kernel H computes partial gradient corresponding to H given

void calc_energy(int n, int m, int k, double* F , double part_grad[n][m] , double H[3][3] ){
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel
    for(int i = 1 ; i < n-k ; i++){
        for(int j = 1 ; j < m-k ; j++){
            for(int u = -k ; u <= k; u++){
                for(int v = -k ; v <= k ; v++){
                    part_grad[i][j] += H[u+k][v+k]*F[(i+u)*m + j+v];
                }
            }
          //calculate absolute value of each element in partial derivative of channel F 
          part_grad[i][j] = abs(part_grad[i][j]);
        }
    }
}

//assumes channels is padded with 0
//calculates the cumulative sum over the individual channel energies
//returns the energy map result over color image of size 3 x n x m 
void calc_RGB_energy(int n, int m, double* channels, double* result){

    //fixed kernels 
    double H_y[3][3] = {
    {-1,2,-1},
    {0,0,0},
    {1,2,1}};

  double H_x[3][3] = {
    {-1,0,1},
    {-2,0,2},
    {-1,0,1}};
    
    int k= 1;

    double partial_x[3][n][m]; //3d parital derivative in x
    double partial_y[3][n][m]; //3d partial derivative in y 

    //calculate the parital derivatives 
    for(int i = 0 ; i < 3 ; i ++){
       calc_energy(n,m,k,channels + n*m*i, partial_x[i], H_x);
       calc_energy(n,m,k,channels + n*m*i, partial_y[i], H_y);
    }

    //calculate the total 3d energy 
    
      for(int j = 0 ; j < n ; j ++){
        for(int k = 0 ; k < m ; k++){
          for(int i = 0 ; i < 3 ; i ++){
            //add elementwise along the z axis 
          *(result+m*j+k) += partial_x[i][j][k] + partial_y[i][j][k];
        }
      } 
    }
}


//first method called before anything else to create a 0 frame aronund original image
//maybe we can have this body be placed into the methods before we want to calculate the energy.  
double* padd0_image(int n, int m, double* channels){
  double padded_image[3][n+2][m+2];
  for(int i = 0 ; i < 3 ; i++){
    for(int j = 0 ; j < n+2 ; j++){
      for(int k = 0 ; k < m+2 ; k++){
        //if the column is 0 or m+1 or the row is 0 or n+1 we set 0 otherwise copy the value 
        if(i == 0 || j == 0 || j == n+1 || i == m+1){
          padded_image[i][j][k] = 0;
        }else
        padded_image[i][j][k] = channels[i*n*m + m*j + k];
      }
    }
  }
  return padded_image;
}

int main()
{
    return 0;
}
