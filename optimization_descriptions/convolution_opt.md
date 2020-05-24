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

### Extension
So we've seen that with 8 registers, we can calculate the convolution for 5 pixels. Thus, we can extend this and use all the 16 registers to calculate the convolution for 10 pixels.

## 3. Blocking
Describe (+ stride)

### Blocking

The blocking was only done for L1 cache and registers because L2/L3 were of size
256KB and 2MB, which was larger than any pitcure that could be given as input. 

The intuitive way in which we could exploit blocking is by having the convolution 
be computed over smaller partitions of the image of the same size, and combinginig 
these small results into the result of the convolution of the original image. 

The problem with this is that the convolutions done on these smaller parts are not
independent w.r.t to convolution computation of the neighbouring blocks like in the 
matrix multiplication case showed throughout the course. Namely we have that the 
very last row and last column of the smaller blocks serve as padding, and hence
if the block is of size M X N we would be computing the convolution over M-2 X N-2
pixels. What this implies, is that if the blocks are disjoint, there would be pixels
of our image whose convolution wouldn't be computed as they are using entries from
two distinct blocks.

Hence the intuitive idea needs to be modified to the idea of blocks that overlap with
their bottom and top block neighbour in 2 rows, and the left and right neighbour
blocks overlap by 2 columns. In this way the neighbour blocks would compute the 
convolution of the pixels that serve as padding to the current block. 

Now that we know how the blocks would subdivide the image, we turn to the analysis of 
the working set, in order to find out what the M and N should be for our image.
Deciding M and N for the image of course immidiatelly decides the size of the block 
in the partial matrix since, as we said before, the convolution is performed over
M-2 x N-2 amount of pixels, where these results correspond to the M-2 x N-2 pixels
in the partial matrix (forming a way of blocking for the partial matrix). 

When we perform the convolution of a pixel we make use of the 8 neighbour pixels 
around it. If we decide to move horizontally, and compute the pixel at position 
x+1,y where x,y is the position of the just computed pixel, we would reuse the pixels
in the central and right column of the first pixel. Similarily if we continue to 
calculate the x+2,y pixel convolution we would reuse the right column of the x,y
pixel and central right column of the x+1,y pixel. Therefore we of course wish to
keep the central and right column of pixel x,y until they are not used anymore for
computing values in the same row. (therefore considering only the same row we would 
want to keep these values until we reach x+3,y) 

But is it only useful to keep the entries of x,y only until they are not used in the same
row anymore? 

Well if the unpadded image consists of only one row, yes, otherwise no!
We wish to keep it until its never reused in our entire block, and if we blocked well,
until its never used  again in our entire image.

Turns out that, when we move on to compute the convolution of pixels in the second row,
we reuse all the central and bottom rows used to compute convolution of the first row pixels.
Similarily, when we move to compute the convolution of the pixels in the 3rd row we use the 
bottom most row of the first row convolution and bottom row of the 2nd row convolution.
This would make our working set be 3 consecutive rows in the image matrix and one pixel
in the partial matrix. 

However, this again needs to be modified because this size of the working set doesn't
allow for a clean replacement of the first row in the working set used to compute the 
convolution of the first row, with a new row representing the bottom most row of the 
convolution of the 2nd row. Therefore we need a working set of 4(M) + M-2 where M is 
the width of the block. Since this working set is not bound in the height, we have 
blocks of size M x height_of_image and M-2 x height_of_partial for the padded image
and partial respectivelly. 


Now when we want to do blocking for registers, we need to think of how we will 
traverse between between the register blocks. For the L1 blocks it was trivial since
the blocks are of the same height as the original matrices. 

If we let the same working set be equal for the registers, then since the 
working set is of size 4M + M-2 = 16 (no of registers is 16) then we would get 3 
as the width of our block. What this would imply is that then the blocks are of the 
form 3 x height_of_image for the padded image and 1 x height_of_partial for the partial. 

This would imply that we would calculate the convolution of all the pixels that belong
to the first column and only then move to the 2nd column and do convolution of these.

However this creates a problem since the L1 cache contains pixels of 4 consecutive rows 
of the image and one row of the partial, NOT COLUMNS. Therefore all values that can be accessed by the
registers from L1 cache aren't used. 

To fix this problem we realise that for the register blocks we should restrain the height to 3 rows
and that we should go between the blocks in a horizontal way to make use of the saved values in L1.

Since we restrain the height to 3 rows, that means that the convolution computed over the pixels that belong 
to that block will always be in the same row, and hence we dont need to account for the 4th row in the 
working set anymore. Hence we have that the working set is 3M + 2 = 16 => M=4. 

Hence we have now blocks of size 4 in width and 3 in height for the image and 2 in width and 1 in height
for the partial. 


## 4. Scalar replacement

## 5. Inlining