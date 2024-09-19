# %%
import sys
sys.path.append('/Users/125715/python')
import math as m
import numpy as np
from fluidlib import thermal_solver as ts
import matplotlib.pyplot as plt
import mat_lib as p
import conversion_library as c

#
n = 1
nodes= np.empty(n,dtype=object)
convPath = np.empty(n, dtype=object)
radPath = np.empty(n, dtype=object)

'''Define Nodes'''
nodes[0]= ts.Node(300, "SS316", "SOLID",V = 0.01, Eg = 500, isothermal=0)
# nodes[1] = ts.Node(500, "SS316", "SOLID", V = 0.01,isothermal=1)
# nodes[2] = ts.Node(300, "WATER", "FLUID", V = 0.01,Eg=200, P =101325)
# nodes[3] = ts.Node(400, 'HYDROGEN', "FLUID", V = 0.1, P = 202322)

'''Define Paths'''
if n > 1:
    '''Conduction Path(s)'''
    condPath = np.empty(n-1, dtype=object)
    condPath[0] = ts.Path(Acond= 0.01, k = 16.3, dx = 0.25)
    condPath[1] = ts.Path(Acond = 0.1, k = 45, dx = 0.25)
    # condPath[2] = ts.Path(Acond = 0.1, k = .01, dx = 0.13)
'''Convection Path(s)'''
convPath[0] = ts.Path(Aconv = 0.2, h = 0, Tinf = 298.15)
# convPath[1] = ts.Path(Aconv = 0.2, h = 10, Tinf = 298.15)
# convPath[2] = ts.Path(Aconv = 0.1, h = 5, Tinf = 298.15)
# convPath[3] = ts.Path(Aconv = 0.1, h = 10)
'''Radiation Path(s)'''
radPath[0] = ts.Path(Arad  = 0, e = 0.0, Tsurr = 298.15)
# radPath[1] = ts.Path(Arad = 0.1, e = 0.0 , Tsurr=298.15)
# radPath[2] = ts.Path(Arad = 0.1, e = 0.0 , Tsurr=298.15)
# radPath[3] = ts.Path()

'''Define Time Window'''
t = [0,900]
t_eval = np.linspace(t[0], t[1], 100)
'''Create Node and Path Matrices to pass into T_vs_t solver'''
nodeMatrix = ts.nodeMatrix(nodes)
if n > 1:
    sol = ts.T_vs_t(t, t_eval = t_eval, nodeMatrix=nodeMatrix,condPath=ts.pathMatrix(condPath),convPath=ts.pathMatrix(convPath), radPath = ts.pathMatrix(radPath))
else:
     sol = ts.T_vs_t(t, t_eval=t_eval, nodeMatrix=nodeMatrix, convPath=ts.pathMatrix(convPath))
n = len(nodeMatrix[0])
print(sol.y-273.15)
for i in range(n):
    plt.figure(1)
    plt.plot(sol.t, sol.y[i]-273.15,label="T0 = " + str(round(sol.y[i,0]-273.15,2)) + " C")
    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (C)")
    # y = max(sol.y[0,:])
    # plt.hlines(float(max(sol.y[0,:])), 0, 6, colors='red', linestyles='--', label='T = ' +str(float(max(sol.y[0,:]))) + ' K')
    plt.legend()

# %%
# %%
