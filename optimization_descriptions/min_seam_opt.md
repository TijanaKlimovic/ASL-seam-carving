# Optimalizations done on min_seam.c

min_seam.c identifies one minimum seam (horizontal / vertical) based on the energy map of an image.

## Track 1 optimizations

### 1. Scalar replacement

### 2. Unrolling

Done for the last row of the dynamic programming matrix, when the index of the min seam start is found. Similar problem to the unrolling for rotation, one iteration already contains multiple operations, more exactly 3. there are 4 execution units available. Unrolling of 2 can sometimes be better, but not guaranteed. 

## Track 2 optimizations

### 1. Rotate energy map when we remove horizontal seam
When we remove a horizontal seam then we need to traverse the whole energy map in column major order for the DP for selecting the seam with the minimal cost. Instead, at horizontal seam removal at the beginning we rotate tho whole energy map, so with this we converted it to the vertical seam removal case, where we need to traverse the energy map in row major order which is cache friendly.
Even though this extra rotation adds an overhead, at bigger images we clearly gain improvement using this. It reduced the cycles needed compared to the previous version to 98-90% of it.


### 2. SIMD (AVX2)
We checked and saw that compiling with -O3 -march=native the compiler vectorizess (AVX2) min_seam.c. We also tried vectorizing by hand so maybe we can achieve better result but our results were the same or maybe slightly worse.

We vectorized:

* the dynamic programming part that computes the minimum seams for the entire energy map (considers minimum from 3 possible choices on row above, basic SIMD, move row by row at a stride of 8)
* the find of the minimum seam (done on the last row, meaning on the total seam value, basic SIMD, move onfinal row at a stride of 8)

### 3. Scalar replacement