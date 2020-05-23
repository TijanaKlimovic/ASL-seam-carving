# Optimalizations done on min_seam.c

## 1. Rotate energy map when we remove horizontal seam
When we remove a horizontal seam then we need to traverse the whole energy map in column major order for the DP for selecting the seam with the minimal cost. Instead, at horizontal seam removal at the beginning we rotate tho whole energy map, so with this we converted it to the vertical seam removal case, where we need to traverse the energy map in row major order which is cache friendly.
Even though this extra rotation adds an overhead, at bigger images we clearly gain improvement using this. It reduced the cycles needed compared to the previous version to 98-90% of it.


## 2. SIMD (AVX2)
We checked and saw that compiling with -O3 -march=native the compiler vectorizess (AVX2) min_seam.c. We also tried vectorizing by hand so maybe we can achieve better result but our results were the same or maybe slightly worse.

## 3. Scalar replacement

