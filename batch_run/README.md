# Batch Run Script

`br.sh` can be used to run the executable on multiple image sizes or multiple images to get statistics.


## Syntax

```
./br.sh <img_file/img_folder> <src_dir> <output_file.txt>
```


### Arguments

*  `<img_file/img_folder>`: Either the address of directory of images or the address of a single image file
	*  `<img_file>`: If a single image is specified, the executable is run on the same image with different sizes as defined by [parameters](#Parameters).
	*  `<img_folder>`: If a directory is specified, the executable is run on all the images in the directory and they are resized to the same size as defined by [parameters](#Parameters).
*  `<src_dir>`: the address of the source code to compile in order to generate the executable.
*  `<output_file.txt>`: Results will be **appended** to this file in table format.


## Dependencies

`br.sh` uses `imagemagick` package.


## Parameters

You can change the script based on your needs by tuning the parameters at the beginning section of the code.
*  `w_red_p` and `h_red_p` (width/height reduction percentage): The number of columns/rows to remove as a percentage of the input image size.
*  `w_step` and `h_step`: Specifies the stride of the iteration over width/height. **Only used in `<img_file>` mode.** For example if input image is 500x400, `w_step=100`, and `h_step=150`, the executable would be run on the same image with sizes 500x400, 400x250, and 300x100 respectively.
*  `w_const` and `h_const`: Specifies the width/height to resize all the images to. **Only used in `<img_folder>` mode.** In this mode we would not have iteration over different image sizes.
*  `gcc flags`: GCC Flags! :)
*  `timing`: Whether to run the executable in timing mode or not.


## Example

Let's say we are in the current directory.
```
./br.sh ../validation/no_alpha_images/ ../ output.txt
```
Would compile the source code in `../` and run it on all the images in `../validation/no_alpha_images/` tuned by the [parameters](#Parameters) and the result can be seen in `output.txt`. If `output.txt` already exists, the results are appended to it, otherwise a new file would be created. Moreover, `inputs` and `outputs` directory would be created in current directory which contain inputs and outputs of all the runs respectively.

Moreover, in `stdout`, you can see the command being run and the instruction counts of different parts.

