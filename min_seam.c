#include <stdio.h>
#include <stdlib.h>

#define MIN2(X, Y, M, IDX) if (X < Y) {M = X; IDX = 0;} else {M = Y; IDX = 1;}

#define MIN3(X, Y, Z, M, IDX) if ((Z < X) && (Z < Y)) {M = Z; IDX = 2;} else {MIN2(X, Y, M, IDX)}


void min_seam(int rsize, int csize, const double *img, const double *e1, int isVer, int **retBacktrack) {
	double *theM = (double *) malloc(rsize * csize * sizeof(double));
	int *backtrack = (*retBacktrack);

	for (int i = 1; i < rsize; i++) { //start from second row
		
		for (int j = 0; j < csize; j++) {
			int minIdx;
			double minVal;
			double minEnergy;

			if (j == 0) {
				MIN2(img[(i - 1) * csize + j], 
					img[(i - 1) * csize + j + 1], 
					minVal, minIdx)
				minVal = img[(i - 1) * csize + j + minIdx];
				backtrack[i * csize + j] = minIdx;
				minEnergy = img[(i - 1) * csize + minIdx];
			} else {
				MIN3(img[(i - 1) * csize + j - 1], 
					img[(i - 1) * csize + j], 
					img[(i - 1) * csize + j + 1], 
					minVal, minIdx)
				minVal = img[(i - 1) * csize + j - 1 + minIdx];
				backtrack[i * csize + j] = minIdx + j - 1;
				minEnergy = img[(i - 1) * csize + minIdx + j - 1];
			}

			theM[i * csize + j] = e1[i * csize + j] + minEnergy;
		}

	}


	free(theM);
}

void test() {
	int idx, idx2;
	double val, val2;
	MIN2(1, 10, val, idx)
	MIN3(-4, 9, -5, val2, idx2)
	printf("[%d] = %lf\n[%d] = %lf\n", idx, val, idx2, val2);
	double *mat = (double *) malloc(5 * 3 * sizeof(double));

	for (int i = 0; i < 5; i++) {

		for (int j = 0; j < 3; j++) {
			mat[i * 3 + j] = (i * 3 + j) / 2.0;
		}

	}

	for (int i = 0; i < 5; i++) {

		for (int j = 0; j < 3; j++) {
			printf("%lf ", mat[i * 3 + j]);
		}

		printf("\n");
	}

	free(mat);
}

int main(int argc, char const *argv[]) {
	test();
	return 0;
}