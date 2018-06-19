import pylab
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import scipy
from qutip import *
from scipy import linalg
from math import *

def nndist(x, n):
    return np.exp(-4*pi*x*x*x*n/3.)*4.*pi*x*x*n
