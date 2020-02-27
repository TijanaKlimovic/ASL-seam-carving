# asl_seam_carving

Optimizing the hell out of seam carving.

## 1. Image parsing  
Used lib: stb_image.h  
Compile: gcc -Wall parse_img.c  -o parse_img -lm  
Usage: ./parse_img <input file path>  
Expects a PNG image and dissects it to its R, G, B channels (allocating a separate width x height matrix for each).  

Aydin added.
