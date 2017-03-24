# coding=utf-8
#MEASURE = measure
__author__ = 'Anastasia Bazhutina'

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import sys
#sys.path.append('/usr/lib/python2.7')
print sys.path
import skimage
import skimage.io as io
from skimage import measure
import scipy.io
from skimage import filters
from skimage import util
from skimage import color
from skimage.morphology import square
from skimage.morphology import dilation
from skimage.morphology import erosion
from skimage.filters.rank import mean
from skimage.morphology import disk
from skimage.morphology import remove_small_holes
from skimage.morphology import remove_small_objects



from mpl_toolkits.mplot3d.art3d import Poly3DCollection
print skimage.__version__