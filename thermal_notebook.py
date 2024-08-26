# %%
import sys
sys.path.append('/Users/125715/python')
import math as m
import numpy as np
from fluidlib import thermal_solver as ts
import matplotlib.pyplot as plt

# %%



# %%

# %%
'''Define Nodes'''
T1 = ts.Node(340, Ac = 0.1,rho=800, c=500, h = 0, e = 0, V =0.1)
T2 = ts.Node(1000, Ac = 0.1, c=500, Eg = 0,h = 20)
T3 = ts.Node(300, Ac = 0.1, c= 500, h = 20, e= 0.2)
# T4 = ts.Node(130, Ac = 0.1,m=1, c=500)
# T5 = ts.Node(700, Ac = 0.4, m = 2, c = 200, q = 10, Eg = 0, h = 20, e = .6)
'''Define Paths'''
P1 = ts.Path(.01,k= 45,dx = .01)
P2 = ts.Path(.01, k = 45, dx = .01)
# P3 = ts.Path(.01, k = 45, dx = .01)
# P4 = ts.Path(.5, k = 25, dx = .02)
t = [0,10000]
t_eval = np.linspace(t[0], t[1], 3600)
nodeMatrix = ts.TMatrix((T1,T2,T3))
pathMatrix = ts.pMatrix((P1,P2),n=2)
sol = ts.T_vs_t(t, t_eval = t_eval,nodeMatrix=nodeMatrix, pathMatrix=pathMatrix)
n = len(nodeMatrix[0])
for i in range (n):
    plt.figure(1)
    plt.plot(sol.t, sol.y[i],label=str(sol.y[i,0]))
    plt.legend()

# %%
ts.ThermalSolve()
# %%
