# %%
import sys
sys.path.append('/Users/125715/python')
import math as m
import numpy as np

from fluidlib import thermal_solver as ts
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d



# %%
'''
Stability Criterion: Fo*(1 + Bi) <= .5
Bi = h*dx/k
Fo = alpha*t/dx**2
'''
n = 2
t = 400
step = 1
dt = int(t/step)
T = np.zeros(shape=(n,dt))
T[:,0] = np.random.random_integers(90,700, n)
Tinf = 298
dx = 10
h = 50
k = 50
print(ts.ThermalSolve(h, dx, k, 0.001, dt, T, Tinf))
t = np.linspace(0,t,dt)
for i in range(len(T[:,0])):
    plt.plot(t, T[i, :])


# %%
