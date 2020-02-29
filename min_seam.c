#include <stdio.h>
#include <stdlib.h>

//calculate minimum of an array (returning both value and index)
double min(double *in, int size, int *idx) {
	double min = in[0]; //could potentially cause errors if size = 0
	(*idx) = 0;

	for (int i = 1; i < size; i++) {

		if (in[i] < min) {
			min = in[i];
			(*idx) = i;
		}

	}

	return min;
}

double **create_2d_array(int rsize, int csize) {
	double **ret = (double **) malloc(rsize * sizeof(double *));

	for (int i = 0; i < rsize; i++) {
		ret[i] = (double *) malloc(csize * sizeof(double));
	}

	return ret;
}

void delete_2d_array(int rsize, int csize, double **mat) {

	for (int i = 0; i < rsize; i++) {
		free(mat[i]);
	}

	free(mat);
}

void min_seam(int rsize, int csize, double **img) {
	double M[100][100];
	int backtrack[100][100];
	double e1[100][100];

	for (int i = 1; i < rsize; i++) {
		
		for (int j = 0; j < csize; j++) {
			int minIdx;
			double minVal;
			double minEnergy;

			if (j == 0) {
				minVal = min(&img[i - 1][j], 2, &minIdx);
				backtrack[i][j] = minIdx;
				minEnergy = img[i - 1][minIdx];
			} else {
				minVal = min(&img[i - 1][j - 1], 3, &minIdx);
				backtrack[i][j] = minIdx + j - 1;
				minEnergy = img[i - 1][minIdx + j - 1];
			}

			M[i][j] = e1[i][j] + minEnergy;
		}

	}

}

int main(int argc, char const *argv[]) {
	return 0;
}