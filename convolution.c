#include <stdio.h>
#include <stdlib.h>

//size of neighbourhood;
int k=1;


//testing RGB matrices 
//red matrix 
double R[4][4] = {
    {0,0,0,0},
    {0, 11, 12, 0},
    {0, 15, 16, 0},
    {0,0,0,0}
}; 

//green matrix
double G[4][4] = {
    {0,0,0,0},
    {0, 11, 12, 0},
    {0, 15, 16, 0},
    {0,0,0,0}
}; 

//blue matrix
double B[4][4] = {
    {0,0,0,0},
    {0, 11, 12, 0},
    {0, 15, 16, 0},
    {0,0,0,0}
}; 

double H_y[3][3] = {
  {1,2,1},
  {0,0,0},
  {-1,-2,-1}
};

double H_x[3][3] = {
  {-1,0,1},
  {-2,0,2},
  {-1,0,1}
};


//assuming that preprocessing is made of 0 padding 
//n rows, m columns of channel F
double** calc_energy(int n, int m, int k, double F[n][m]){
    double part_grad_x[n][m];
    double part_grad_y[n][m];
    double** energy;
    
    energy = malloc(n*sizeof(*energy));
    for(int i = 0; i < n; ++i)
        energy[i] = malloc(m*sizeof(**energy));
    
    //start at 1 and end at n-1/m-1 to avoid padding
    // i,j are the current pixel
    for(int i = 1 ; i < n-k ; i++){
        for(int j = 1 ; j < m-k ; j++){
            for(int u = -k ; u <= k; u++){
                for(int v = -k ; v <= k ; v++){
                    part_grad_x[i][j] += H_x[u+k][v+k]*F[i+u][j+v];
                    part_grad_x[i][j] += H_y[u+k][v+k]*F[i+u][j+v];
                }
            }
            //calculate energy of pixel i,j in channel F
            energy[i][j] = abs(part_grad_x[i][j]) + abs(part_grad_y[i][j]);
        }
    }
    return energy; 
}

//calculates the cumulative sum over the individual channel energies
//each channel is padded with a 0  frame 
//retuns the energy map over color image of size n x m 
double** calc_RGB_energy(int n, int m, int k, double channels[3][n][m]){
    double** channel_energy;
    //energy matrix
    double** result;
    
    result = malloc(n*sizeof(*result));
    for(int i = 0; i < n; ++i)
        result[i] = malloc(m*sizeof(**result));
     
    for(int i = 0 ; i < 3 ; i ++){
       channel_energy = calc_energy(n,m,k,channels[i]);
       for(int j = 0 ; j < n ; j++){
           for(int z = 0 ; z < m ; z++){
               result[j][z] += channel_energy[j][z];
           }
       }
       free(channel_energy);
    }
    return result; 
}

int main()
{
    //need to free the array returned by calc_RGB_energy
    return 0;
}
