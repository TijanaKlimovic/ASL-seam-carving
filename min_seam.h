#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "count.h"

double min_seam(int rsize, int csize, double *img, int is_ver, int *ret_backtrack);