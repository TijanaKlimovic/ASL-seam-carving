 
`./time_test.sh` can be used to make measurements for creating performance plots.

`imgs/` contains images grouped by aspect ratio.

## Usage

```
./time_test.sh <img_file> <lower_bound_img_width> <upper_bound_img_width> <step> <seams_to_remove>"
```

## Example

```
./time_test.sh imgs/1_5/9.png 200 600 100 10
```

Takes imgs/1_5/9.png and resized it from width of 200 to 600 pixels (keeping the aspect ratio) by steps of 100 pixel. At each resized image it removes 10 horizontal and 10 vertical seams. For each resized image it outputs the measured cycles, flops (adds + mults) and from these the calculated performance [flops/cycle].
