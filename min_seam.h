#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "count.h"

int min_seam(int rsize, int csize, unsigned char *img, int is_ver, int *ret_backtrack);