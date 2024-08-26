# %%
import sys
import math as m
import numpy as np
import thermal_solver as ts
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


# %%
T = ts.node(198, 8000, .05, 500)
T2 = ts.node(400, 8000, .05, 500)
p = ts.path(1, h=15)
args = p.getArgs() + T.getArgs()
T_array = [T.getT(), T2.getT()]
t = [1, 10000]
T = solve_ivp(ts.ODE, t, T_array, args=args)
print(T.y)
plt.plot(T.t, T.y[1,:])







# %%
