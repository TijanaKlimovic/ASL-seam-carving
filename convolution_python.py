import sys

import numpy as np
from imageio import imread, imwrite
from scipy.ndimage.filters import convolve

# tqdm is not strictly necessary, but it gives us a pretty progress bar
# to visualize progress.
from tqdm import trange

def calc_energy(img):
    filter_du = np.array([
        [1.0, 2.0, 1.0],
        [0.0, 0.0, 0.0],
        [-1.0, -2.0, -1.0],
    ])
    # This converts it from a 2D filter to a 3D filter, replicating the same
    # filter for each channel: R, G, B
    # filter_du = np.stack([filter_du] * 3, axis=2)

    filter_dv = np.array([
        [1.0, 0.0, -1.0],
        [2.0, 0.0, -2.0],
        [1.0, 0.0, -1.0],
    ])
    # This converts it from a 2D filter to a 3D filter, replicating the same
    # filter for each channel: R, G, B
    # filter_dv = np.stack([filter_dv] * 3, axis=2)

    img = img.astype('float32')
    R = img[:, :, 0]
    G = img[:, :, 1]
    B = img[:, :, 2]
    # convolved = np.absolute(convolve(img, filter_du)) + np.absolute(convolve(img, filter_dv))
    convolvedR = np.absolute(convolve(R, filter_du, mode = 'constant', cval = 0.0)) + np.absolute(convolve(R, filter_dv, mode = 'constant', cval = 0.0))
    convolvedG = np.absolute(convolve(G, filter_du, mode = 'constant', cval = 0.0)) + np.absolute(convolve(G, filter_dv, mode = 'constant', cval = 0.0))
    convolvedB = np.absolute(convolve(B, filter_du, mode = 'constant', cval = 0.0)) + np.absolute(convolve(B, filter_dv, mode = 'constant', cval = 0.0))

    # We sum the energies in the red, green, and blue channels
    # energy_map = convolved.sum(axis=2)
    # for row in energy_map:
    # 	print(row)
    #print("the type is: " , type(energy_map));
    energy_map = convolvedR + convolvedG + convolvedB
    
    return energy_map


def main():
    #scale = float(sys.argv[1])
    in_filename = sys.argv[1] #file path
    #out_filename = sys.argv[3]

    img = imread(in_filename)
    calc_energy(img)
    #out = crop_c(img, scale)
    #imwrite(out_filename, out)

if __name__ == '__main__':
    main()