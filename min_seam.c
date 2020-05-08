#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "convolution.h"

//--------------------	counter for instructions -------------------

#ifdef count_instr 
extern unsigned long long add_count;	//count the total number of add instructions
extern unsigned long long mult_count; 	//count the total number of mult instructions
#endif

//------------------------------------------------------------------


#define MIN2(X, Y, M, IDX) if (X < Y) {M = X; IDX = 0;} else {M = Y; IDX = 1;}

#define MIN3(X, Y, Z, M, IDX) if ((X < Y) && (X < Z)) {M = X; IDX = -1;} else {MIN2(Y, Z, M, IDX)}

// #define LOG(X) X
#define LOG(X)

int min_seam(int height, int width, int *img, int is_ver, int *ret_backtrack) {
	int *the_m = (int *) malloc(height * width * sizeof(int));

	int size = 3*(height+2)*(width+2);
	int* padded_img = (int*) malloc(size*sizeof(int));

	//int padded_img[3][n+2][m+2];
	for(int i = 0 ; i < 3 ; i++) {
	   	for(int j = 0 ; j < height+2 ; j++) {
	     	for(int k = 0 ; k < width+2 ; k++) {
	        //if the column is 0 or m+1 or the row is 0 or n+1 we set 0 otherwise copy the value 
	        if(j == 0 || k == 0 || j == height+1 || k == width+1) {
	        	padded_img[i*(height+2)*(width+2) + (width+2)*j + k] = 0;
	        } else {
	        	padded_img[i*(height+2)*(width+2) + (width+2)*j + k] = img[i*height*width + width*(j-1) + k-1];
	        }
	      	}
	    }
	}

	int padded_h = height + 2;
	int padded_w = width + 2;
	  
	//Sobel filters 
	int H_y[3][3] = {
	    {-1,-2,-1},
	    {0,0,0},
	    {1,2,1}};

	int H_x[3][3] = {
	    {-1,0,1},
	    {-2,0,2},
	    {-1,0,1}};

	size = 3 * padded_h * padded_w ;

	int* partial_x = (int*) malloc( size*sizeof(int));
	int* partial_y = (int*) malloc( size*sizeof(int));

    //calculate the parital derivatives 
    for(int i = 0 ; i < 3 ; i ++){
      	//pass the ith channel for energy calculation
       	int *channel = padded_img + padded_h * padded_w * i;
       	int *part_grad_x = partial_x + padded_h * padded_w * i;

	    //start at 1 and end at n-1/m-1 to avoid padding
	    // i,j are the current pixel
	    for (int i = 0; i < padded_h; i++) {
	        for (int j = 0; j < padded_w; j++) {
	            *(part_grad_x + i*padded_w + j) = 0;
	        }
	    }

	    for(int i = 1 ; i < padded_h-1 ; i++){
	        for(int j = 1 ; j < padded_w-1 ; j++){
	            for(int u = -1 ; u <= 1; u++){
	                for(int v = -1 ; v <= 1 ; v++){
	                    *(part_grad_x + i*padded_w + j) += H_x[u+1][v+1]*channel[(i+u)*padded_w + j+v];
	                }
	            }

	            //calculate absolute value of each element in partial derivative of channel F 
	            if(*(part_grad_x + i*padded_w + j) < 0){
	              *(part_grad_x + i*padded_w + j) = (-1) * (*(part_grad_x + i*padded_w + j));
	            }
	        }
	    }

       int *part_grad_y = partial_y + padded_h * padded_w * i;

	    //start at 1 and end at n-1/m-1 to avoid padding
	    // i,j are the current pixel
	    for (int i = 0; i < padded_h; i++) {
	        for (int j = 0; j < padded_w; j++) {
	            *(part_grad_y + i*padded_w + j) = 0;
	        }
	    }

	    for(int i = 1 ; i < padded_h-1 ; i++){
	        for(int j = 1 ; j < padded_w-1 ; j++){
	            for(int u = -1 ; u <= 1; u++){
	                for(int v = -1 ; v <= 1 ; v++){
	                    *(part_grad_y + i*padded_w + j) += H_y[u+1][v+1]*channel[(i+u)*padded_w + j+v];
	                }
	            }

	            //calculate absolute value of each element in partial derivative of channel F 
	            if(*(part_grad_y + i*padded_w + j) < 0){
	              *(part_grad_y + i*padded_w + j) = (-1) * (*(part_grad_y + i*padded_w + j));
	            }
	        }
	    }
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            the_m[width * i + j] = 0;
        }
    }

    //calculate the total 3d energy 
      for(int j = 1 ; j < padded_h-1 ; j ++){
        for(int k = 1 ; k < padded_w-1 ; k++){
          for(int i = 0 ; i < 3 ; i ++){
            //add elementwise along the z axis 
          *(the_m+ width*(j-1)+k-1) += *(partial_x + i*padded_h*padded_w + j*padded_w + k) + *(partial_y + i*padded_h*padded_w + j*padded_w + k);
        }
      } 
    }

    free(partial_x);
    free(partial_y);

	// contains index of the value from the prev row/column from where we came here
	int *backtrack = (int *) malloc(height * width * sizeof(int)); //different from what we return

	// find vertical min seam
	if (is_ver) {
		for (int i = 1; i < height; i++) { //start from second row	
			for (int j = 0; j < width; j++) {

				int where = i * width + j;
				int where_before = where - width;
				int min_idx;
				int min_val;

				// first col
				if (j == 0) {
					MIN2(the_m[where_before], 
						 the_m[where_before + 1], 
						 min_val, 
						 min_idx)

					backtrack[where] = min_idx;
				// last col
				} else if (j == width - 1) {
					MIN2(the_m[where_before - 1],
						 the_m[where_before],
						 min_val,
						 min_idx)

					min_idx--;

					backtrack[where] = j + min_idx;
				} else {
					MIN3(the_m[where_before - 1], 
						 the_m[where_before], 
						 the_m[where_before + 1], 
						 min_val, 
						 min_idx)

					backtrack[where] = j + min_idx;
				}
				the_m[where] += min_val;
			}
		}

		//process the data to return in appropriate format
		int ret = INT_MAX;
		int direction = -1; //used in turning 2D backtrack into 1D

		//find the index of the minimum value of last row in the dp matrix
		int last_row = (height - 1)  * width;
		for (int cnt = 0; cnt < width; cnt++) {
			int current = last_row + cnt;
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}
		}

		//return the 1D backtrack (only the min seam)
		// direction -= last_start;

		for (int i = height - 1; i >= 0; i--) {
			ret_backtrack[i] = direction;
			direction = backtrack[last_row + direction];
			last_row -= width;
		}

		free(the_m);
		free(backtrack);
		free(padded_img);
		
		return ret;

	} else {
	// find horizontal min seam
		for (int i = 1; i < width; i++) { //start from second col
		
			for (int j = 0; j < height; j++) {

				int where = j * width + i;
				int where_before = where - 1;
				int min_idx;
				int min_val;

				// first row
				if (j == 0) {
					MIN2(the_m[where_before], 
						 the_m[where_before + width], 
						 min_val, 
						 min_idx)

					backtrack[where] = min_idx;

				// last row
				} else if (j == height - 1) {
					MIN2(the_m[where_before - width],
						 the_m[where_before],
						 min_val,
						 min_idx)

					min_idx--;

					backtrack[where] = j + min_idx;
				} else {
					MIN3(the_m[where_before - width], 
						 the_m[where_before], 
						 the_m[where_before + width], 
						 min_val, 
						 min_idx)

					backtrack[where] = j + min_idx;
				}
				the_m[where] += min_val;
			}
		}
		//process the data to return in appropriate format
		int ret = INT_MAX;
		int direction = -1; //used in turning 2D backtrack into 1D

		//find the minimum of last row/col of the_m
		//set the counters to the beginning of the last row/col
		int last_col = width - 1;
		
		for (int cnt = 0; cnt < height; cnt++) {
			int current = last_col + (cnt * width);
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}
		}

		//return the 1D backtrack (only the min seam)
		// direction -= last_start;

		for (int i = width - 1; i >= 0; i--) {
			ret_backtrack[i] = direction;
			direction = backtrack[last_col + (direction * width)];
			last_col -= 1;
		}
		free(the_m);
		free(backtrack);
		free(padded_img);
		
		return ret;
	}
}
