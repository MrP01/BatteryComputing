# Python code for solving chronoamperometre problem

from scipy.sparse import dia_matrix, dok_matrix
from scipy.sparse.linalg import spsolve
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Set up for solving the heat equation for a
# Boundary values for space and time
a = 0 ; b = 10
T0 = 0
T = 1


dx = 0.1
dt = 0.01
BC2 = 1
IC = 1
BC1 = 0