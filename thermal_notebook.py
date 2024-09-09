# %%y
import sys
sys.path.append('/Users/125715/python')
import math as m
import numpy as np
from fluidlib import thermal_solver as ts
import matplotlib.pyplot as plt
import mat_lib as p
import conversion_library as c

# # %%
'''Define Nodes'''
n = 1
nodes= np.empty(n,dtype=object)
paths = [0,0,0]
nodes[0] = ts.Node(293, mat= "SS316", Ac = c.si2sm(76.969),Ar = c.si2sm(76.969), h = 50, V=c.cintocm(53.26), q = 2520)

'''Define Paths'''
P1 = ts.Path(.01, k= 45,dx = .01)
P2 = ts.Path(.01, k = 45, dx = .01)
'''Define Time Window'''
t = [0,6]
t_eval = np.linspace(t[0], t[1], 400)
'''Create Node and Path Matrices to pass into T_vs_t solver'''
nodeMatrix = ts.TMatrix(nodes)
pathMatrix = paths
sol = ts.T_vs_t(t, t_eval = t_eval,nodeMatrix=nodeMatrix, pathMatrix=pathMatrix)
n = len(nodeMatrix[0])
print(sol.y[0,388])
print(sol.t)
for i in range (n):
    plt.figure(1)
    plt.plot(sol.t, sol.y[i],label="T0 = " + str(sol.y[i,0]) + "K")
    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (K)")
    y = max(sol.y[0,:])
    plt.hlines(float(max(sol.y[0,:])), 0, 6, colors='red', linestyles='--', label='T = ' +str(float(max(sol.y[0,:]))) + ' K')
    plt.legend()

# %%
# %%
