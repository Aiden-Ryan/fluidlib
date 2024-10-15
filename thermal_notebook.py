

# %%
# import sys
# sys.path.append('/Users/125715/python/fluidlib')
import math as m
import matplotlib.pyplot as plt
import numpy as np
import conversion_library as c; 
import volume_calc as v
import mat_lib as mat
import thermal_solver as t
from scipy.integrate import solve_ivp
# from fluidlib import flow_prop_calc as f
import gc

node1 = t.Node(
    T = 300,
    medium="SS316",
    medium_type="SOLID",
    V = 0.001
)
node2 = t.Node(
    T = 400,
    medium="SS316",
    medium_type="SOLID",
    V = 0.001 
)
node3 = t.Node(
    T = 350, 
    medium="SS316", 
    medium_type="SOLID",
    V = 0.001
)


a = t.Path(
    nodeA = node1,
    nodeB = node2,
    Acond=0.1,
    dx=.02
)
b = t.Path(
    nodeA = node1,
    nodeB= node3,
    Acond=0.1,
    dx = .3
)

n = 6
paths = (a,b)
nodes = (node1, node2,node3)
T = np.array([node1.T, node2.T, node3.T])
t_span=[0,1000]
t_eval = np.linspace(t_span[0], t_span[1], 1000)
sol = t.T_vs_t(t_span, t_eval, paths,nodes)
print(sol.y)
# for i in range(n):
#     plt.figure(1)
#     plt.plot(sol.t, sol.y[i]-273.15,label="T" + str(i+1) + " = " + str(round(sol.y[i,0]-273.15,2)) + " C")
#     plt.xlabel("Time (s)") 
#     plt.ylabel("Temperature (C)")
# plt.legend()