// compile with : gcc -Wall parse_img.c  -o parse_img -lm
// usage: ./parse_img <filename>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


int load_image(const char *filename);

void save_image(char *filename, int width, int height);