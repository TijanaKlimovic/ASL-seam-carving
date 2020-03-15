#include <stdio.h>
#include <stdlib.h>

//calculate minimum of an array (returning both value and index)
double min(const double *in, int size, int *idx) {
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

void min_seam(int rsize, int csize, const double **img, const double **e1, double ***retM, int ***retBacktrack) {
	double **theM = (*retM);
	int **backtrack = (*retBacktrack);

	for (int i = 1; i < rsize; i++) { //start from second row
		
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

			theM[i][j] = e1[i][j] + minEnergy;
		}

	}

}

void test() {
	double a[] = {4.1, 2.3, 5.5, 4.5, 8.2, 0.9, 10, 1.2, 3.3, 3, 3, 2, 9.0};
	int idx;
	double m = min(a, 13, &idx);
	printf("[%d] = %lf\n", idx, m);
	double **mat = create_2d_array(5, 3);

	for (int i = 0; i < 5; i++) {

		for (int j = 0; j < 3; j++) {
			mat[i][j] = (i * 3 + j) / 2.0;
		}

	}

	for (int i = 0; i < 5; i++) {

		for (int j = 0; j < 3; j++) {
			printf("%lf ", mat[i][j]);
		}

		printf("\n");
	}

	double am[5][3] = {4.1, 2.3, 5.5, 4.5, 8.2, 0.9, 10, 1.2, 3.3, 3, 3, 2, 9.0, 4, -2.4};

	for (int i = 0; i < 5; i++) {
		printf("%p: \n", am[i]);
		printf("%p: \n", am[i][0]);

		for (int j = 0; j < 3; j++) {
			printf("%lf ", am[i][j]);
		}

		printf("\n");
	}
	
}

int main(int argc, char const *argv[]) {

	return 0;
}