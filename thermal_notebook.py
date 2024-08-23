# %%
import sys
sys.path.append('/Users/125715/python')
import math as m
import numpy as np
from fluidlib import thermal_solver as ts
import matplotlib.pyplot as plt




# %%
'''
Stability Criterion: Fo*(1 + Bi) <= .5
Bi = h*dx/k
Fo = alpha*t/dx**2
'''
n = 2
t = 3600
step = 1
dt = int(t/step)
T = np.zeros(shape=(n,dt))
T[:,0] = [200,150]
Tinf = 298
dx = 1
h = 50
k = 45
A=.01
print(ts.ThermalSolve(h, dx, k, 0.001, dt, T, Tinf))
t = np.linspace(0,t,dt)
for i in range(len(T[:,0])):
    plt.plot(t, T[i, :])


# %%
rho = 7980
T1 = ts.node(300, 7980, .1)
T2 = ts.node(100, rho, 0.1)
P1 = ts.path(.01, k = 25, dx=1)

T_array = ts.nodeArray()
T_array.build_array(T1,T2)
print(T_array.get_array())
solution = ts.buildNetwork(T_array, P1)
# %%
