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
T = ts.node(198)
T2 = ts.node(400)
p = ts.path(1, h=2,L=2,k=3)
args = p.getArgs() + T.getArgs()
T_array = [T.getT(), T2.getT()]
T_array = [165,400]
t = [1,2]
T = solve_ivp(ts.ODE, t, T_array, args=args)
print(T.y)
plt.plot(T.t, T.y[0,:])
plt.plot(T.t, T.y[1,:])







# %%
