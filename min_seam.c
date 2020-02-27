#include <stdio.h>

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

void min_seam() {
	//placeholders for variables
	double img[100][100];
	int rsize = 100;
	int csize = 100;
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