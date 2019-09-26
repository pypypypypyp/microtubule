#!/usr/bin/python
#coding: utf-8

import numpy as np
from EMAN2 import *

CutOffResolution = 10.
APIX = 1.32
file = "particles4.mrcs"
N = 58999
SIZE = 600

center_x = 0
center_y = SIZE/2

x, y = np.meshgrid(np.arange(SIZE/2+1), np.arange(SIZE))
y = np.where(y > SIZE/2, y-SIZE, y)
mask = np.sqrt(x*x + y*y) / (SIZE * APIX)

mask_random   = mask > 1.0/CutOffResolution
mask_nochange = mask < 1.0/CutOffResolution

for i in range(N):
        e = EMData(file, i)
        data = EMNumPy.em2numpy(e)
        fft = np.fft.rfft2(data)
        random = np.random.rand(SIZE, SIZE/2+1)*np.pi*2
        random_phi = np.exp(random*1j)
        mask = mask_random*random_phi + mask_nochange
        fftmask = fft*mask
        EMNumPy.numpy2em(np.fft.irfft2(fftmask).astype(np.float32)).write_image("random.mrcs", -1)
        if i%1000 == 0: print i
