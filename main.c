//#error Please comment out the next two lines under linux, then comment this error
//#include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler
//#include <windows.h> // Include if under windows

#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#include "tsc_x86.h"

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
//#define FREQUENCY 3.6e9
#define CALIBRATE

//--------------- necessary for benchmarking---------------

#include <stdio.h>
#include "parse_img.h"

//--------------- necessary for run---------------

#include <string.h>

//--------------- extern count variables --------------

#include "count.h"

#ifdef count_instr 
unsigned long long add_count = 0;	//count the total number of add instructions
unsigned long long mult_count = 0; //count the total number of mult instructions
#endif

//---------------------------------------------------------

#define MIN2(X, Y, M, IDX) if (X < Y) {M = X; IDX = 0;} else {M = Y; IDX = 1;}

#define MIN3(X, Y, Z, M, IDX) if ((X < Y) && (X < Z)) {M = X; IDX = -1;} else {MIN2(Y, Z, M, IDX)}

void init_tsc() {
	; // no need to initialize anything for x86
}

myInt64 start_tsc(void) {
    tsc_counter start;
    CPUID();
    RDTSC(start);
    return COUNTER_VAL(start);
}

myInt64 stop_tsc(myInt64 start) {
	tsc_counter end;
	RDTSC(end);
	CPUID();
	return COUNTER_VAL(end) - start;
}

// Data structure to hold information in a T cell
struct cell_T {
	// optimal cost to get to this dimension
	int optimal_cost;
	// corresponding matrix for this dimension,
	// obtained based on optimal cost
	// size is 3 * (width - row) * (height - column)
	int *i;
};

int run(int width, int height, int *image, const char *output_file_name, int width_diff, int height_diff ) {
	width_diff++;
	height_diff++;

	struct cell_T *T = (struct cell_T *)malloc(width_diff * height_diff * sizeof(struct cell_T));
	T[0].optimal_cost = 0;
	T[0].i = image;

	// calc first row -> vertical seam only
	int j;
	for (j = 1; j < width_diff; ++j) {
		int T_index_left = j - 1;
		int *image_left = T[T_index_left].i;
		int image_width = width - j + 1;

		int *backtrack_indexes = (int *)malloc(height * sizeof(int));
	
		int optimal_cost;
		int h = height;
		int w = image_width;
		int *img = image_left;

		


		int *the_m = (int *) malloc(h * w * sizeof(int));

		int size = 3*(h+2)*(w+2);
		int* padded_img = (int*) malloc(size*sizeof(int));

		//int padded_img[3][n+2][m+2];
		for(int i = 0 ; i < 3 ; i++) {
		   	for(int j = 0 ; j < h+2 ; j++) {
		     	for(int k = 0 ; k < w+2 ; k++) {
		        //if the column is 0 or m+1 or the row is 0 or n+1 we set 0 otherwise copy the value 
		        if(j == 0 || k == 0 || j == h+1 || k == w+1) {
		        	padded_img[i*(h+2)*(w+2) + (w+2)*j + k] = 0;
		        } else {
		        	padded_img[i*(h+2)*(w+2) + (w+2)*j + k] = img[i*h*w + w*(j-1) + k-1];
		        }
		      	}
		    }
		}

		int padded_h = h + 2;
		int padded_w = w + 2;
		  
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

	    for (int i = 0; i < h; i++) {
	        for (int j = 0; j < w; j++) {
	            the_m[w * i + j] = 0;
	        }
	    }

	    //calculate the total 3d energy 
	      for(int j = 1 ; j < padded_h-1 ; j ++){
	        for(int k = 1 ; k < padded_w-1 ; k++){
	          for(int i = 0 ; i < 3 ; i ++){
	            //add elementwise along the z axis 
	          *(the_m+ w*(j-1)+k-1) += *(partial_x + i*padded_h*padded_w + j*padded_w + k) + *(partial_y + i*padded_h*padded_w + j*padded_w + k);
	        }
	      } 
	    }

	    free(partial_x);
	    free(partial_y);

		// contains index of the value from the prev row/column from where we came here
		int *backtrack = (int *) malloc(h * w * sizeof(int)); //different from what we return

		// find vertical min seam
		for (int i = 1; i < h; i++) { //start from second row	
			for (int j = 0; j < w; j++) {

				int where = i * w + j;
				int where_before = where - w;
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
				} else if (j == w - 1) {
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
		int last_row = (h - 1)  * w;
		for (int cnt = 0; cnt < w; cnt++) {
			int current = last_row + cnt;
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}
		}

		optimal_cost = ret;

		//return the 1D backtrack (only the min seam)
		// direction -= last_start;

		for (int i = h - 1; i >= 0; i--) {
			backtrack_indexes[i] = direction;
			direction = backtrack[last_row + direction];
			last_row -= w;
		}

		free(the_m);
		free(backtrack);
		free(padded_img);




		T[j].optimal_cost = T[T_index_left].optimal_cost + optimal_cost;
		T[j].i = (int *)malloc(3 * height * (image_width - 1) * sizeof(int));

		int k, l;
		for (k = 0; k < height; ++k) { // construct each row at a time
			int crr_col = 0;
			for (l = 0; l < image_width; ++l)
				if (l != backtrack_indexes[k]) { // check for elem to remove
					T[j].i[k*(image_width-1)+crr_col] = T[T_index_left].i[k*image_width + l];
					T[j].i[height*(image_width-1)+k*(image_width-1)+crr_col] = T[T_index_left].i[height*image_width+k*image_width + l];
					T[j].i[2*height*(image_width-1)+k*(image_width-1)+crr_col] = T[T_index_left].i[2*height*image_width+k*image_width + l];
					crr_col++;
				}
		}
		free(backtrack_indexes);
	}
		
	// going through the rest of the dp matrix
	int i;
	for (i = 1; i < height_diff; ++i) {
		// first column -> horizontal seam only
		int T_index = i * width_diff;
		int T_index_up = (i - 1) * width_diff;
		int *image_up = T[T_index_up].i;
		int image_height = height - i + 1;
	
		int *backtrack_indexes = (int *)malloc(width * sizeof(int));
	
		int optimal_cost;
		int h = image_height;
		int w = width;
		int *img = image_up;

		int *the_m = (int *) malloc(h * w * sizeof(int));

		int size = 3*(h+2)*(w+2);
		int* padded_img = (int*) malloc(size*sizeof(int));

		//int padded_img[3][n+2][m+2];
		for(int i = 0 ; i < 3 ; i++) {
		   	for(int j = 0 ; j < h+2 ; j++) {
		     	for(int k = 0 ; k < w+2 ; k++) {
		        //if the column is 0 or m+1 or the row is 0 or n+1 we set 0 otherwise copy the value 
		        if(j == 0 || k == 0 || j == h+1 || k == w+1) {
		        	padded_img[i*(h+2)*(w+2) + (w+2)*j + k] = 0;
		        } else {
		        	padded_img[i*(h+2)*(w+2) + (w+2)*j + k] = img[i*h*w + w*(j-1) + k-1];
		        }
		      	}
		    }
		}

		int padded_h = h + 2;
		int padded_w = w + 2;
		  
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

	    for (int i = 0; i < h; i++) {
	        for (int j = 0; j < w; j++) {
	            the_m[w * i + j] = 0;
	        }
	    }

	    //calculate the total 3d energy 
	      for(int j = 1 ; j < padded_h-1 ; j ++){
	        for(int k = 1 ; k < padded_w-1 ; k++){
	          for(int i = 0 ; i < 3 ; i ++){
	            //add elementwise along the z axis 
	          *(the_m+ w*(j-1)+k-1) += *(partial_x + i*padded_h*padded_w + j*padded_w + k) + *(partial_y + i*padded_h*padded_w + j*padded_w + k);
	        }
	      } 
	    }

	    free(partial_x);
	    free(partial_y);

		// contains index of the value from the prev row/column from where we came here
		int *backtrack = (int *) malloc(h * w * sizeof(int)); //different from what we return

		// find horizontal min seam
		for (int i = 1; i < w; i++) { //start from second col
		
			for (int j = 0; j < h; j++) {

				int where = j * w + i;
				int where_before = where - 1;
				int min_idx;
				int min_val;

				// first row
				if (j == 0) {
					MIN2(the_m[where_before], 
						 the_m[where_before + w], 
						 min_val, 
						 min_idx)

					backtrack[where] = min_idx;

				// last row
				} else if (j == h - 1) {
					MIN2(the_m[where_before - w],
						 the_m[where_before],
						 min_val,
						 min_idx)

					min_idx--;

					backtrack[where] = j + min_idx;
				} else {
					MIN3(the_m[where_before - w], 
						 the_m[where_before], 
						 the_m[where_before + w], 
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
		int last_col = w - 1;
		
		for (int cnt = 0; cnt < h; cnt++) {
			int current = last_col + (cnt * w);
			if (the_m[current] < ret) {
				ret = the_m[current];
				direction = cnt;
			}
		}

		optimal_cost = ret;

		//return the 1D backtrack (only the min seam)
		// direction -= last_start;

		for (int i = w - 1; i >= 0; i--) {
			backtrack_indexes[i] = direction;
			direction = backtrack[last_col + (direction * w)];
			last_col -= 1;
		}
		free(the_m);
		free(backtrack);
		free(padded_img);
		





		T[T_index].optimal_cost = T[T_index_up].optimal_cost + optimal_cost;
		T[T_index].i = (int *)malloc(3 * (image_height - 1) * width * sizeof(int));

		int k, l;
		for (k = 0; k < width; ++k) { // construct each column at a time
			int crr_row = 0;
			for (l = 0; l < image_height; ++l)
				if (l != backtrack_indexes[k]) { // check for elem to remove
					T[T_index].i[crr_row*width+k] = T[T_index_up].i[l*width + k];
					T[T_index].i[(image_height-1)*width+crr_row*width+k] = T[T_index_up].i[image_height*width+l*width + k];
					T[T_index].i[2*(image_height-1)*width+crr_row*width+k] = T[T_index_up].i[2*image_height*width+l*width + k];
					crr_row++;
				}
		}
		free(backtrack_indexes);

		// continue from the 2nd column
		for (j = 1; j < width_diff; ++j) {
			// height -> horizontal seam -> row
			// width -> vertical seam -> column
			int T_index = i * width_diff + j;
			int T_index_left = i * width_diff + j - 1;
			int *image_left = T[T_index_left].i;
			int image_left_width = width - j + 1;
			int image_left_height = height - i;

			int *backtrack_left = (int *)malloc(image_left_height * sizeof(int));
			//int optimal_cost_left = min_seam(image_left_height, image_left_width, image_left, 1, backtrack_left);


			//------------------
		
			int optimal_cost_left;
			int h = image_left_height;
			int w = image_left_width;
			int *img = image_left;

			


			int *the_m = (int *) malloc(h * w * sizeof(int));

			int size = 3*(h+2)*(w+2);
			int* padded_img = (int*) malloc(size*sizeof(int));

			//int padded_img[3][n+2][m+2];
			for(int i = 0 ; i < 3 ; i++) {
			   	for(int j = 0 ; j < h+2 ; j++) {
			     	for(int k = 0 ; k < w+2 ; k++) {
			        //if the column is 0 or m+1 or the row is 0 or n+1 we set 0 otherwise copy the value 
			        if(j == 0 || k == 0 || j == h+1 || k == w+1) {
			        	padded_img[i*(h+2)*(w+2) + (w+2)*j + k] = 0;
			        } else {
			        	padded_img[i*(h+2)*(w+2) + (w+2)*j + k] = img[i*h*w + w*(j-1) + k-1];
			        }
			      	}
			    }
			}

			int padded_h = h + 2;
			int padded_w = w + 2;
			  
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

		    for (int i = 0; i < h; i++) {
		        for (int j = 0; j < w; j++) {
		            the_m[w * i + j] = 0;
		        }
		    }

		    //calculate the total 3d energy 
		      for(int j = 1 ; j < padded_h-1 ; j ++){
		        for(int k = 1 ; k < padded_w-1 ; k++){
		          for(int i = 0 ; i < 3 ; i ++){
		            //add elementwise along the z axis 
		          *(the_m+ w*(j-1)+k-1) += *(partial_x + i*padded_h*padded_w + j*padded_w + k) + *(partial_y + i*padded_h*padded_w + j*padded_w + k);
		        }
		      } 
		    }

		    free(partial_x);
		    free(partial_y);

			// contains index of the value from the prev row/column from where we came here
			int *backtrack = (int *) malloc(h * w * sizeof(int)); //different from what we return

			// find vertical min seam
			for (int i = 1; i < h; i++) { //start from second row	
				for (int j = 0; j < w; j++) {

					int where = i * w + j;
					int where_before = where - w;
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
					} else if (j == w - 1) {
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
			int last_row = (h - 1)  * w;
			for (int cnt = 0; cnt < w; cnt++) {
				int current = last_row + cnt;
				if (the_m[current] < ret) {
					ret = the_m[current];
					direction = cnt;
				}
			}

			optimal_cost_left = ret;

			//return the 1D backtrack (only the min seam)
			// direction -= last_start;

			for (int i = h - 1; i >= 0; i--) {
				backtrack_left[i] = direction;
				direction = backtrack[last_row + direction];
				last_row -= w;
			}

			free(the_m);
			free(backtrack);
			free(padded_img);

			//------------------

			int T_index_up = (i - 1) * width_diff + j;
			int *image_up = T[T_index_up].i;
			int image_up_width = width - j;
			int image_up_height = height - i + 1;

			int *backtrack_up = (int *)malloc(image_up_width * sizeof(int));

			// ----------------

		{
			int optimal_cost_up;
			int h = image_up_height;
			int w = image_up_width;
			int *img = image_up;

			int *the_m = (int *) malloc(h * w * sizeof(int));

			int size = 3*(h+2)*(w+2);
			int* padded_img = (int*) malloc(size*sizeof(int));

			//int padded_img[3][n+2][m+2];
			for(int i = 0 ; i < 3 ; i++) {
			   	for(int j = 0 ; j < h+2 ; j++) {
			     	for(int k = 0 ; k < w+2 ; k++) {
			        //if the column is 0 or m+1 or the row is 0 or n+1 we set 0 otherwise copy the value 
			        if(j == 0 || k == 0 || j == h+1 || k == w+1) {
			        	padded_img[i*(h+2)*(w+2) + (w+2)*j + k] = 0;
			        } else {
			        	padded_img[i*(h+2)*(w+2) + (w+2)*j + k] = img[i*h*w + w*(j-1) + k-1];
			        }
			      	}
			    }
			}

			int padded_h = h + 2;
			int padded_w = w + 2;
			  
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

		    for (int i = 0; i < h; i++) {
		        for (int j = 0; j < w; j++) {
		            the_m[w * i + j] = 0;
		        }
		    }

		    //calculate the total 3d energy 
		      for(int j = 1 ; j < padded_h-1 ; j ++){
		        for(int k = 1 ; k < padded_w-1 ; k++){
		          for(int i = 0 ; i < 3 ; i ++){
		            //add elementwise along the z axis 
		          *(the_m+ w*(j-1)+k-1) += *(partial_x + i*padded_h*padded_w + j*padded_w + k) + *(partial_y + i*padded_h*padded_w + j*padded_w + k);
		        }
		      } 
		    }

		    free(partial_x);
		    free(partial_y);

			// contains index of the value from the prev row/column from where we came here
			int *backtrack = (int *) malloc(h * w * sizeof(int)); //different from what we return

			// find horizontal min seam
			for (int i = 1; i < w; i++) { //start from second col
			
				for (int j = 0; j < h; j++) {

					int where = j * w + i;
					int where_before = where - 1;
					int min_idx;
					int min_val;

					// first row
					if (j == 0) {
						MIN2(the_m[where_before], 
							 the_m[where_before + w], 
							 min_val, 
							 min_idx)

						backtrack[where] = min_idx;

					// last row
					} else if (j == h - 1) {
						MIN2(the_m[where_before - w],
							 the_m[where_before],
							 min_val,
							 min_idx)

						min_idx--;

						backtrack[where] = j + min_idx;
					} else {
						MIN3(the_m[where_before - w], 
							 the_m[where_before], 
							 the_m[where_before + w], 
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
			int last_col = w - 1;
			
			for (int cnt = 0; cnt < h; cnt++) {
				int current = last_col + (cnt * w);
				if (the_m[current] < ret) {
					ret = the_m[current];
					direction = cnt;
				}
			}

			optimal_cost_up = ret;

			//return the 1D backtrack (only the min seam)
			// direction -= last_start;

			for (int i = w - 1; i >= 0; i--) {
				backtrack_up[i] = direction;
				direction = backtrack[last_col + (direction * w)];
				last_col -= 1;
			}
			free(the_m);
			free(backtrack);
			free(padded_img);
		
		





			// remove column
			if (T[T_index_left].optimal_cost + optimal_cost_left <=
				T[T_index_up].optimal_cost + optimal_cost_up) {
				T[T_index].optimal_cost = T[T_index_left].optimal_cost + optimal_cost_left;
				T[T_index].i = (int *)malloc(3 * image_left_height * (image_left_width - 1) * sizeof(int));
				
				int k, l;
				for (k = 0; k < image_left_height; ++k) { // construct each row at a time
					int crr_col = 0;
					for (l = 0; l < image_left_width; ++l)
						if (l != backtrack_left[k]) { // check for elem to remove
							T[T_index].i[k*(image_left_width-1)+crr_col] = T[T_index_left].i[k*image_left_width + l];
							T[T_index].i[image_left_height*(image_left_width-1)+k*(image_left_width-1)+crr_col] = T[T_index_left].i[image_left_height*(image_left_width)+k*image_left_width + l];
							T[T_index].i[2*image_left_height*(image_left_width-1)+k*(image_left_width-1)+crr_col] = T[T_index_left].i[2*image_left_height*(image_left_width)+k*image_left_width + l];
							crr_col++;
						}
				}

			// remove row
			} else {
				T[T_index].optimal_cost = T[T_index_up].optimal_cost + optimal_cost_up;
				T[T_index].i = (int *)malloc(3 * (image_up_height - 1) * image_up_width * sizeof(int));

				int k, l;
				for (k = 0; k < image_up_width; ++k) { // construct each column at a time
					int crr_row = 0;
					for (l = 0; l < image_up_height; ++l)
						if (l != backtrack_up[k]) { // check for elem to remove
							T[T_index].i[crr_row*image_up_width+k] = T[T_index_up].i[l*image_up_width + k];
							T[T_index].i[(image_up_height-1)*image_up_width+crr_row*image_up_width+k] = T[T_index_up].i[image_up_height*image_up_width+l*image_up_width + k];
							T[T_index].i[2*(image_up_height-1)*image_up_width+crr_row*image_up_width+k] = T[T_index_up].i[2*image_up_height*image_up_width+l*image_up_width + k];
							crr_row++;
						}
				}
			}
			free(backtrack_left);
			free(backtrack_up);
		
		}
	}
	}
		
	// copy 
	int *res = malloc(3 * (width - width_diff + 1) * (height - height_diff + 1) * sizeof(int));
	memcpy(res, T[width_diff * height_diff - 1].i, 3 * (width - width_diff + 1) * (height - height_diff + 1) * sizeof(int));
	for (int i = 1; i < width_diff * height_diff; ++i) {
		free(T[i].i);
	}

	free(T);

	save_image(output_file_name, width - width_diff + 1, height - height_diff + 1, res);
	free(res);
	return 0;
}

/* 
benchmarking function using rdtsc instruction
*/

double rdtsc(int width, int height, int *output, const char *output_file_name, int width_diff, int height_diff) {
    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;


#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i) {
            run(width, height, output, output_file_name, width_diff, height_diff);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
		run(width, height, output, output_file_name, width_diff, height_diff);	
	}
    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}

/* main function , has two options, first is the image path second is -t/<anything except nothing> */ 

int main(int argc, char const *argv[]) {

	if (argc < 6) {
		printf("Usage: %s <image_path> <output_file_name> <width_diff> <height_diff> <timing boolean>\n", argv[0]);
		return 1;
	}
	// normal run

	int width_diff = atoi(argv[3]);
	int height_diff = atoi(argv[4]);
	int width, height;
	int *output;

	if (!load_image(argv[1], &width, &height, &output)) {
		return 1;
	}

	if(strcmp(argv[5],"1")){
		int out = run(width, height, output, argv[2], width_diff, height_diff);
		free(output);
		return out;
	}
	//time it
	double r = rdtsc(width, height, output, argv[2], width_diff, height_diff);

	free(output);
    printf("RDTSC instruction: %.0lf cycles measured\n", r);
    return 0;
}
