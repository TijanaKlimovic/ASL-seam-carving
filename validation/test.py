import sys

from tqdm import trange
import numpy as np
from imageio import imread, imwrite
from scipy.ndimage.filters import convolve

def crop_c(img, diff):
    for i in trange(diff):
        img = carve_column(img)

    return img

def crop_r(img, diff):
    img = np.rot90(img, 1, (0, 1))
    img = crop_c(img, diff)
    img = np.rot90(img, 3, (0, 1))
    return img

def carve_column(img):
    r, c, _ = img.shape

    M, backtrack = minimum_seam(img)
    mask = np.ones((r, c), dtype=np.bool)

    j = np.argmin(M[-1]) # min cost
    for i in reversed(range(r)):
        mask[i, j] = False
        j = backtrack[i, j]

    mask = np.stack([mask] * 3, axis=2)
    img = img[mask].reshape((r, c - 1, 3))
    return img

def minimum_seam(img):
    r, c, _ = img.shape
    energy_map = calc_energy(img)

    M = energy_map.copy()
    backtrack = np.zeros_like(M, dtype=np.int)

    for i in range(1, r):
        for j in range(0, c):
            # Handle the left edge of the image, to ensure we don't index a -1
            if j == 0:
                idx = np.argmin(M[i-1, j:j + 2])
                backtrack[i, j] = idx + j
                min_energy = M[i-1, idx + j]
            else:
                idx = np.argmin(M[i - 1, j - 1:j + 2])
                backtrack[i, j] = idx + j - 1
                min_energy = M[i - 1, idx + j - 1]

            M[i, j] += min_energy

    return M, backtrack

def test_pad_img(img, fname):
    r = img[:,:,0]
    g = img[:,:,1]
    b = img[:,:,2]
    padded_r = np.pad(r, 1, mode='constant')
    padded_g = np.pad(g, 1, mode='constant')
    padded_b = np.pad(b, 1, mode='constant')

    with open(fname, "w") as f:
        for img in [padded_r, padded_g, padded_b]:
            h, w = img.shape
            for i in range(h):
                for j in range(w):
                    f.write(str(int(img[i, j])) + " ")
                f.write("\n")
            f.write("\n")

def calc_energy(img):
    """
    Modified energy map calculation.

    Calculates energy with Sobel filter on each channel then sums them.
    """
    img = img.astype('float32')

    R = img[:,:,0]
    G = img[:,:,1]
    B = img[:,:,2]

    filter_du = np.array([
        [1.0, 2.0, 1.0],
        [0.0, 0.0, 0.0],
        [-1.0, -2.0, -1.0],
    ])
    
    filter_dv = np.array([
        [1.0, 0.0, -1.0],
        [2.0, 0.0, -2.0],
        [1.0, 0.0, -1.0],
    ])

    # Sobel filter
    Gx = np.array([[1, 0, -1], [2, 0, -2], [1, 0, -1]])
    Gy = np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]])
    convolved_r = np.absolute(convolve(R, filter_dv, mode='constant', cval=0)) + np.absolute(convolve(R, filter_du, mode='constant', cval=0))
    convolved_g = np.absolute(convolve(G, filter_dv, mode='constant', cval=0)) + np.absolute(convolve(G, filter_du, mode='constant', cval=0))
    convolved_b = np.absolute(convolve(B, filter_dv, mode='constant', cval=0)) + np.absolute(convolve(B, filter_du, mode='constant', cval=0))

    energy_map = convolved_r + convolved_g + convolved_b

    return energy_map

def test_calc_energy(img, fname):
    """
    Modified energy map calculation.

    Calculates energy with Sobel filter on each channel then sums them.
    """
    img = img.astype('float32')

    R = img[:,:,0]
    G = img[:,:,1]
    B = img[:,:,2]

    filter_du = np.array([
        [1.0, 2.0, 1.0],
        [0.0, 0.0, 0.0],
        [-1.0, -2.0, -1.0],
    ])
    
    filter_dv = np.array([
        [1.0, 0.0, -1.0],
        [2.0, 0.0, -2.0],
        [1.0, 0.0, -1.0],
    ])

    # Sobel filter
    Gx = np.array([[1, 0, -1], [2, 0, -2], [1, 0, -1]])
    Gy = np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]])
    convolved_r = np.absolute(convolve(R, filter_dv, mode='constant', cval=0)) + np.absolute(convolve(R, filter_du, mode='constant', cval=0))
    convolved_g = np.absolute(convolve(G, filter_dv, mode='constant', cval=0)) + np.absolute(convolve(G, filter_du, mode='constant', cval=0))
    convolved_b = np.absolute(convolve(B, filter_dv, mode='constant', cval=0)) + np.absolute(convolve(B, filter_du, mode='constant', cval=0))

    energy_map = convolved_r + convolved_g + convolved_b

    h, w = energy_map.shape
    with open(fname, "w") as f:
        for i in range(h):
            for j in range(w):
                f.write(str(int(energy_map[i, j])) + " ")
            f.write("\n")
        f.write("\n")


    return energy_map

def test_min_seam(img, is_vertical, fname):
    if is_vertical:
        M, backtrack = minimum_seam(img)
        min_cost = M[-1, np.argmin(M[-1])]
    else:
        rotated = np.rot90(img, 1, (0, 1))
        M, backtrack = minimum_seam(rotated)
        min_cost = M[-1, np.argmin(M[-1])]
    with open(fname, "w") as f:
        f.write(str(int(min_cost)) + "\n")

def test_seam_carving(img, is_vertical, diff, fname):
    if is_vertical:
        out = crop_c(img, diff)
    else:
        out = crop_r(img, diff)
    r = out[:,:,0]
    g = out[:,:,1]
    b = out[:,:,2]

    with open(fname, "w") as f:
        for img in [r, g, b]:
            h, w = img.shape
            for i in range(h):
                for j in range(w):
                    f.write(str(int(img[i, j])) + " ")
                f.write("\n")
            f.write("\n")

def main():
    if len(sys.argv) != 2:
        print('usage: test.py <image_in>', file=sys.stderr)
        sys.exit(1)

    in_filename = sys.argv[1]

    img = imread(in_filename)

    # test pad image
    test_pad_img(img, "unit_tests/out_reference/padding.txt")

    # test energy map
    test_calc_energy(img, "unit_tests/out_reference/energy_map.txt")

    # test min seam - vertical
    test_min_seam(img, True, "unit_tests/out_reference/min_seam_cost_v.txt")

    # test min seam - horizontal
    test_min_seam(img, False, "unit_tests/out_reference/min_seam_cost_h.txt")

    # test seam carving - remove 1 vertical seam
    test_seam_carving(img, True, 1, "unit_tests/out_reference/seam_carved_1.txt")

    # test seam carving - remove 1 horizontal seam
    img = imread(in_filename)
    test_seam_carving(img, False, 1, "unit_tests/out_reference/seam_carved_2.txt")

    # test seam carving - remove 2 vertical seams
    img = imread(in_filename)
    test_seam_carving(img, True, 2, "unit_tests/out_reference/seam_carved_3.txt")

    # test seam carving - remove 2 horizontal seams
    img = imread(in_filename)
    test_seam_carving(img, False, 2, "unit_tests/out_reference/seam_carved_4.txt")

if __name__ == '__main__':
    main()