
'''
Assumptions: 
ignore radiation 

'''
import numpy as np
from scipy.integrate import solve_ivp
import mat_lib as p
import conversion_library as c

material_mapping = {
    "SS316": "p.SS316",
    "Al6061": "p.Al6061",
}

class Node:
    '''
    Creates Temperature node with HT attributes
    '''
    def __init__(self, T, medium, medium_type, V= 0.01, Eg = 0.0, Pressure = 0,isothermal = 0,connectedPaths = None, Tsurr = 298.15, e = 0.0,Arad = 0.0): 
        '''
        T = type::float; initial temperature, K
        
        medium_type = "FLUID" or "SOLID"
        
        '''
        if medium_type == "FLUID":
            self.density = c.getFluidProperty('D', 'TP', T, Pressure, medium) #kg/m3 
            self.k = c.getFluidProperty('CONDUCTIVITY', 'TP', T, Pressure, medium) #W/m-K
            self.c = c.getFluidProperty('CPMASS', 'TP', T, Pressure, medium) #J/kgK
        elif medium_type == "SOLID": 
            self.medium = eval(material_mapping.get(medium)) # type: ignore
            self.density = p.getMatProp(self.medium, "density") #kg/m3 
            self.k= p.getMatProp(self.medium, "k") #W/m-K
            self.c = p.getMatProp(self.medium, "specificHeatCapacity") #J/kgK
        self.isothermal = isothermal #isothermal condition is true if isothermal == 1 
        self.V = V #m3
        self.T = T #K 
        self.Eg = Eg #W, heat generated
        self.Tsurr = Tsurr #K
        self.e = e #emissivity
        self.Arad = Arad #m^2
        if connectedPaths is None:
            self.connectedPaths = [] #specifies which paths are connected to node
    
class Path:
    '''
    Creates Conduction, Convection, or Radiation Pathway between Nodes
    Assumptions: Path is of Constant Area
    '''
    def __init__(self, nodeA, nodeB, Acond=0.0, Aconv = 0.0, Arad = 0.0, dx=1.0, h=0, e=0.0, Tinf = 298.15, Tsurr = 298.15):
        self.Acond = Acond #Conduction Area [m^2]
        self.Aconv = Aconv #Convection Area [m^2]
        self.Arad = Arad #Radiation Area [m^2]
        self.dx = dx #length [m]
        self.h = h #Overall heat transfer coefficient [W/m2K]
        self.e = e # dimensionless, emissivity of body's surface
        self.nodeA = nodeA 
        self.nodeB = nodeB
        self.k = (nodeA.k+nodeB.k)/2  #conduction coefficient [W/mK]

def T_vs_t(t, t_eval, path, nodes):
    n = len(nodes)
    T = []
    for node in nodes: 
        T.append(node.T)
    if n > 2:
        for i in range(n):
            for j in range(len(path)):  
                if path[j].nodeA == nodes[i] or path[j].nodeB == nodes[i]: #if either of path j's nodes are equal to current node, append path to connected attribute path
                    nodes[i].connectedPaths.append(j)
    def func(t,T):
        dTdt = np.zeros(n, dtype=float)
        if n == 2: #Only two nodes
            k = path.k
            h = path.h
            dx = path.dx

            Acond = path.Acond
            Aconv = path.Aconv
            
            Arad0 = path.nodeA.Arad
            Arad1 = path.nodeB.Arad

            e0 = path.nodeA.e
            e1 = path.nodeB.e
            
            Tsurr0 = path.nodeA.Tsurr
            Tsurr1 = path.nodeB.Tsurr
            sig = 5.678e-8

            dTdt[0] = (- k * Acond*(T[0] - T[1]) / dx - h * Aconv *(T[0] - T[[1]]) - (sig * e0 * Arad0 * (T[0]**4 - Tsurr0**4)))/(path.nodeA.density*path.nodeA.V*path.nodeA.c)

            dTdt[1] =  (- k * Acond*(T[1] - T[0]) / dx - h * Aconv *(T[1] - T[0]) - (sig * e1 * Arad1 * (T[1]**4 - Tsurr1**4)))/(path.nodeB.density*path.nodeB.V*path.nodeB.c)

            return dTdt
        else:
            for i in range(n):
                # for node i get connectedPaths
                P = nodes[i].connectedPaths
                for j in range(len(P)):
                    #Iterate thru each connected to node i
                    '''Conduction Terms'''
                    k = path[P[j]].k
                    dx = path[P[j]].dx
                    Acond = path[P[j]].Acond
                    '''Convection Terms'''
                    h = path[P[j]].h
                    Aconv = path[P[j]].Aconv
                    '''Radiation Terms'''
                    e = path[P[j]].nodeA.e
                    Tsurr = path[P[j]].nodeA.Tsurr
                    Arad = path[P[j]].nodeA.Arad
                    sig = 5.678e-8
                    #get temperature associated with path[p[j]] order must be maintained. 
                    T1 = T[nodes.index(path[P[j]].nodeA)]
                    T2 = T[nodes.index(path[P[j]].nodeB)]

                    if nodes.index(path[P[j]].nodeA) == i: #If T being accessed is from iteration node; use this equation
                        dTdt[i] = dTdt[i] - (k*Acond*(T1-T2) / dx + h * Aconv * (T1 - T2) + e*sig*Arad*(T1**4 - Tsurr**4))/(path[P[j]].nodeA.density * path[P[j]].nodeA.V * path[P[j]].nodeA.c) if nodes[i].isothermal == 0 else 0
                    elif nodes.index(path[P[j]].nodeA) != i: #If T being accessed is not from iteration node; use this equation
                        dTdt[i] =  dTdt[i] + (k*Acond*(T1-T2) / dx + h * Aconv * (T1 - T2) + e*sig*Arad*(T1**4 - Tsurr**4))/(path[P[j]].nodeA.density * path[P[j]].nodeA.V * path[P[j]].nodeA.c) if nodes[i].isothermal == 0 else 0
            return dTdt 
        
    return solve_ivp(func, t_span=t, t_eval=t_eval, y0=T,method='LSODA',atol=1e-7)
