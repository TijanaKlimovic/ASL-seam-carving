#include <stdio.h>
#include <stdlib.h>

#define MIN2_IDX(X, Y) ((X < Y) ? (0) : (1))
#define MIN3_IDX(X, Y, Z) ( ((Z < X) && (Z < Y)) ? (2) : MIN2_IDX(X, Y) )


void min_seam(int rsize, int csize, const double *img, const double *e1, double **retM, int **retBacktrack) {
	double *theM = (*retM);
	int *backtrack = (*retBacktrack);

	for (int i = 1; i < rsize; i++) { //start from second row
		
		for (int j = 0; j < csize; j++) {
			int minIdx;
			double minVal;
			double minEnergy;

			if (j == 0) {
				minIdx = MIN2_IDX(img[(i - 1) * csize + j], img[(i - 1) * csize + j + 1]);
				minVal = img[(i - 1) * csize + j + minIdx];
				backtrack[i * csize + j] = minIdx;
				minEnergy = img[(i - 1) * csize + minIdx];
			} else {
				minIdx = MIN3_IDX(img[(i - 1) * csize + j - 1], img[(i - 1) * csize + j], img[(i - 1) * csize + j + 1]);
				minVal = img[(i - 1) * csize + j - 1 + minIdx];
				backtrack[i * csize + j] = minIdx + j - 1;
				minEnergy = img[(i - 1) * csize + minIdx + j - 1];
			}

			theM[i * csize + j] = e1[i * csize + j] + minEnergy;
		}

	}

}

void test() {
	int idx = MIN2_IDX(1, 10);
	int idx2 = MIN3_IDX(-4, 9, -5);
	printf("0 =?= %d\n2 =?= %d\n", idx, idx2);
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

}

int main(int argc, char const *argv[]) {
	test();
	return 0;
}