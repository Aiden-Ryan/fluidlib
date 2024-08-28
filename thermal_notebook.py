# %%
import sys
sys.path.append('/Users/125715/python')
import math as m
import numpy as np
from fluidlib import thermal_solver as ts
import matplotlib.pyplot as plt

# %%
'''Define Nodes'''
T1 = ts.Node(340, kcond=45, Ac = 1,rho=800, c=500, h = 10, e = 0, V =0.1)
T2 = ts.Node(1000, kcond = 45,Ac = 1, c=500, Eg = 0,h = 10)
T3 = ts.Node(300, kcond= 45,Ac = 1, c= 500, h = 10, e= 0.2)
'''Define Paths'''
P1 = ts.Path(.01,k= 45,dx = .01)
P2 = ts.Path(.01, k = 45, dx = .01)
'''Define Time Window'''
t = [0,10000]
t_eval = np.linspace(t[0], t[1], 3600)
'''Create Node and Path Matrices to pass into T_vs_t solver'''
nodeMatrix = ts.TMatrix((T1,T2,T3))
pathMatrix = ts.pMatrix((P1,P2),n=2)
sol = ts.T_vs_t(t, t_eval = t_eval,nodeMatrix=nodeMatrix, pathMatrix=pathMatrix)
n = len(nodeMatrix[0])
for i in range (n):
    plt.figure(1)
    plt.plot(sol.t, sol.y[i],label=str(sol.y[i,0]))
    plt.legend()
