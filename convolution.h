#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "count.h"

void calc_RGB_energy(int n, int m, double* channels, double* result);
double* padd0_image(int n, int m, double* channels);
