# Optimalizations done on convoluton.c

## 1. Optimized sobel filter
Describe

## 2. SIMD (AVX2)
### Potential
We leverage the fact that calculating the convolution for different are independent. For example, if we want to calculate the convolution for pixel `(a, b)` we don't need to know the convolution at point `(a, b - 1)` or any other pixel. Thus, there is a very good potential to use SIMD to parallelize convolution calculations. 
### Method
For the convolution of each pixel we need the pixel value at the 8 surrouding pixels. We want to avoid doing different calculations for the subset of values in each register because that would require us extracting those data from the register, performing the calculations and then putting them back which is costly. Thus, for each register, we load the data that needs the same computation. Here's a figure for the pixels: (we are calculating the convolution for the pixel at \*)
```
NW  N   NE
W   *   E 
SW  S   SE
```
If we are calculating the convolution for pixels A and B, the best candidates who need the same computations are the cells that are located the same relative to the pixel A, B; i.e. one register gets loaded with NW pixels of A and then B, another register gets loaded with N pixels of A and then B, and ...
Since each register has a capacity of 16 shorts, and we're doing all the channels together (so NW actually is 3 shorts, one for R, one for G, and another for B), we can fit at most 5 pixel RGB values (`3 x 5 = 15` so we actually waste 1 short per register). As a result, with 8 registers, we can do the convolution computation for 5 pixels together as follows:
```
Reg[0] would contain RGB values of NW cells for the 5 pixels
Reg[1] would contain RGB values of N  cells for the 5 pixels
Reg[2] would contain RGB values of NE cells for the 5 pixels
Reg[3] would contain RGB values of W  cells for the 5 pixels
Reg[4] would contain RGB values of E  cells for the 5 pixels
Reg[5] would contain RGB values of SW cells for the 5 pixels
Reg[6] would contain RGB values of S  cells for the 5 pixels
Reg[7] would contain RGB values of SE cells for the 5 pixels
```
In other words, if the pixels of a part of the image look like this: (each letter is a pixel and has its own R, G, B value)
```
A B C D E F G
H I J K L M N
O P Q R S T U
```
we would calculate the convolution for pixels `I`, `J`, `K`, `L`, and `M` like this:
```
Reg[0]: (NW of I i.e. A) and (NW of J i.e. B) and (NW of K i.e. C) and (NW of L i.e. D) and (NW of M i.e. E): A B C D E
Reg[1]: (N  of I i.e. B) and (N  of J i.e. C) and (N  of K i.e. D) and (N  of L i.e. E) and (N  of M i.e. F): B C D E F
Reg[2]: (NE of I i.e. C) and (NE of J i.e. D) and (NE of K i.e. E) and (NE of L i.e. F) and (NE of M i.e. G): C D E F G
Reg[3]: (W  of I i.e. H) and (W  of J i.e. I) and (W  of K i.e. J) and (W  of L i.e. K) and (W  of M i.e. L): H I J K L
Reg[4]: (E  of I i.e. J) and (E  of J i.e. K) and (E  of K i.e. L) and (E  of L i.e. M) and (E  of M i.e. N): J K L M N
Reg[5]: (SW of I i.e. O) and (SW of J i.e. P) and (SW of K i.e. Q) and (SW of L i.e. R) and (SW of M i.e. S): O P Q R S
Reg[6]: (S  of I i.e. P) and (S  of J i.e. Q) and (S  of K i.e. R) and (S  of L i.e. S) and (S  of M i.e. T): P Q R S T
Reg[7]: (SE of I i.e. Q) and (SE of J i.e. R) and (SE of K i.e. S) and (SE of L i.e. T) and (SE of M i.e. U): Q R S T U
```
then, e.g. for calculting the horizontal sobel filter:
```
-1  0   1
-2  0   2
-1  0   1
```
we need `NE - NW + 2 * (E - W) + SE - SW`, which translated into our registers would look like: `Reg[2] - Reg[0] + 2 * (Reg[4] - Reg[3]) + Reg[7] - Reg[5]` (side note: we implement \*2 operation by adding the number to itself)

## 3. Blocking
Describe (+ stride)

## 4. Scalar replacement

## 5. Inlining