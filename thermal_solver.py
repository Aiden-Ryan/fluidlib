
'''
Assumptions: 
ignore radiation 

'''
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
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
    def __init__(self, order, T, medium, medium_type, V= 0.01, Eg = 0.0, Pressure = 0,isothermal = 0,connectedPaths = None, Tsurr = 298.15, e = 0.0,Arad = 0.0): 
        '''
        T = type::float; initial temperature, K
        
        medium_type = "FLUID" or "SOLID"
        
        order attribute assigns sequence position
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
        self.order = order # type- int. dictates iteration order  
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

def T_vs_t2(t, t_eval, T, path, nodes):
    '''PATHS MUST BE SEQUENCED ACCORDING TO NODE ORDER. EX: path1(node1, node2) NOT path1(node2,node1)'''
    n = len(nodes)
    if n > 2:
        for i in range(n):
            for j in range(len(path)):  
                if path[j].nodeA == nodes[i] or path[j].nodeB == nodes[i]: #if either of path j's nodes are equal to current node, append path to connected attribute path
                    nodes[i].connectedPaths.append(j)
    def func2(t,T):
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

            dTdt[path.nodeA.order] =  - k * Acond*(T[path.nodeA.order] - T[path.nodeB.order]) / (dx*path.nodeA.density*path.nodeA.V*path.nodeA.c) - h * Aconv *(T[path.nodeA.order] - T[[path.nodeB.order]]) / (path.nodeA.density*path.nodeA.V*path.nodeA.c) - (sig * e0 * Arad0 * (T[path.nodeA.order]**4 - Tsurr0**4)) /(path.nodeA.density*path.nodeA.V*path.nodeA.c)

            dTdt[path.nodeB.order] =  - k * Acond*(T[path.nodeB.order] - T[path.nodeA.order]) / (dx*path.nodeB.density*path.nodeA.V*path.nodeA.c) - h * Aconv *(T[path.nodeB.order] - T[path.nodeA.order]) / (path.nodeB.density*path.nodeB.V*path.nodeB.c) - (sig * e1 * Arad1 * (T[path.nodeB.order]**4 - Tsurr1**4)) /(path.nodeB.density*path.nodeB.V*path.nodeB.c)

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
                    T1 = T[path[P[j]].nodeA.order] 
                    T2 = T[path[P[j]].nodeB.order] 
                    if path[j].nodeA.order == i: #If T being accessed is from iteration node; use this equation
                        dTdt[i] = dTdt[i] - k*Acond*(T1-T2) / (dx*path[P[j]].nodeA.density*path[P[j]].nodeA.V*path[P[j]].nodeA.c) - h * Aconv * (T1 - T2) / (path[P[j]].nodeA.density*path[P[j]].nodeA.V * path[P[j]].nodeA.c) - e*sig*Arad*(T1**4 - Tsurr**4)/(path[P[j]].nodeA.density * path[P[j]].nodeA.V * path[P[j]].nodeA.c) if nodes[i].isothermal == 0 else 0
                    elif path[j].nodeA.order != i: #If T being accessed is not from iteration node; use this equation
                        dTdt[i] =  dTdt[i] + k*Acond*(T1-T2) / (dx*path[P[j]].nodeA.density*path[P[j]].nodeA.V*path[P[j]].nodeA.c) + h * Aconv * (T1 - T2) / (path[P[j]].nodeA.density*path[P[j]].nodeA.V * path[P[j]].nodeA.c) + e*sig*Arad*(T1**4 - Tsurr**4)/(path[P[j]].nodeA.density * path[P[j]].nodeA.V * path[P[j]].nodeA.c) if nodes[i].isothermal == 0 else 0
            return dTdt 
        
    return solve_ivp(func2, t_span=t, t_eval=t_eval, y0=T,method='LSODA',atol=1e-7)

# def nodeMatrix(nodes):
    '''
    Creates Node matrix to be passed into T_vs_t function
    '''
    attr = ['T', 'density', 'c', 'q_flux', 'Eg', 'V', 'k','isothermal']
    attribute_arrays=[]
    for attr in attr:
        attribute_arrays.append([getattr(node, attr) for node in nodes])
    return np.vstack(attribute_arrays)

# def pathMatrix(paths): 
    '''
    Creates path attribute matrix to be passed into T_vs_t function
    '''
    attr = ['Acond', 'Aconv', 'Arad', 'k', 'dx', 'h', 'e', 'Tinf', 'Tsurr']
    attribute_arrays=[]
    for attr in attr:
        attribute_arrays.append([getattr(path, attr) for path in paths])
    return np.vstack(attribute_arrays)
 
# def T_vs_t(t_span, t_eval, nodeMatrix, condPath=np.zeros(10), convPath=np.zeros(10), radPath=np.zeros(10)):
    '''
    Iteratively solves Thermal Network. Returns Temperature vs time plot.
    '''
    n = len(nodeMatrix[0, :])
    T = nodeMatrix[0,:]
    rho = nodeMatrix[1,:]
    c = nodeMatrix[2,:]
    q_flux = nodeMatrix[3,:]
    Eg = nodeMatrix[4,:]
    V = nodeMatrix[5,:]
    kcheck = nodeMatrix[6,:]
    isoT = nodeMatrix[7,:]

    if n == 1:
        '''Conduction Path'''
        Acond = condPath[0]
        k = condPath[3]
        dx = condPath[4]
        '''Convection Path'''
        Aconv = convPath[1]
        h = convPath[5]
        Tinf = convPath[7]
        '''Radiation Path'''
        Arad = radPath[2]
        e = radPath[6]
        Tsurr = radPath[8]
    elif n == 2: 
        '''Conduction Path'''
        Acond = condPath[0]
        k = condPath[3]
        dx = condPath[4]
        '''Convection Path'''
        Aconv = convPath[1,:]
        h = convPath[5,:]
        Tinf = convPath[7,:]
        '''Radiation Path'''
        Arad = radPath[2,:]
        e = radPath[6,:]
        Tsurr = radPath[8,:]
    else:
        '''
        if more than one path exists
        '''
        '''Conduction Path'''
        Acond = condPath[0,:]
        k = condPath[3,:]
        dx = condPath[4,:]
        '''Convection Path'''
        Aconv = convPath[1,:]
        h = convPath[5,:]
        Tinf = convPath[7,:]
        '''Radiation Path'''
        Arad = radPath[2,:]
        e = radPath[6,:]
        Tsurr = radPath[8,:]

    '''BIOT CHECK'''
    for i in range(n-1):
        a = h[i]*V[i]/(Aconv[i]*kcheck[i])
        if a > 0.1:
            raise Exception("BIOT NUMBER GREATER THAN 0.1, DOES NOT SATISFY LUMPED CAPACITANCE CRITERIA")
    dTdt = np.zeros(n, dtype=float)

    def func(t, T,rho,c,q_flux,Eg,Tinf,h,e,k,dx,Aconv,V,Acond=np.zeros(n),Arad = np.zeros(n),isoT=np.zeros(n)):
        sig = 5.67e-8 #Stefan-Boltzmann Constant [W/m^2/K^4]
        if n == 2:
            dTdt[0] = (q_flux[0] +Eg[0] - k[0] * Acond[0] * (T[0] - T[1])/(dx) 
            - h[0]*Aconv[0]*(T[0] - Tinf[0]) - sig*e[0]*Arad[0]*(T[0]**4 - Tsurr[0]**4))/(rho[0]*c[0]*V[0]) if isoT[0] == 0 else 0
            dTdt[1] = (q_flux[1] +Eg[1] + k[0] * Acond[0] * (T[0] - T[1])/(dx) 
            - h[1]*Aconv[1]*(T[1] - Tinf[1]) - sig*e[1]*Arad[1]*(T[1]**4 - Tsurr[1]**4))/(rho[1]*c[1]*V[1]) if isoT[1] == 0 else 0
            return dTdt       
        elif n == 1:
            dTdt[0] = (q_flux+Eg[0]- h*Aconv*(T - Tinf) - sig*e*Arad*(T**4 - Tsurr**4))/(rho[0]*c[0]*V[0]) if isoT == 0 else 0
            return dTdt
        else:
            for i in range(n):
                if i == 0:
                    dTdt[i] = (q_flux[i] + Eg[i] -k[i]*Acond[i]*(T[i] - T[i+1])/dx[i] 
                            -h[i]*Aconv[i]*(T[i] - Tinf[i])- sig*e[i]*Arad[i]*(T[i]**4 - Tsurr[i]**4))/(rho[i] * c[i]*V[i]) if isoT[i] == 0 else 0
                elif i == n-1:
                    dTdt[i] = (q_flux[i] + Eg[i] -k[i-1]*Acond[i-1]*(T[i] - T[i-1])/dx[i-1] 
                            -h[i]*Aconv[i]*(T[i] - Tinf[i])- sig*e[i]*Arad[i]*(T[i]**4 - Tsurr[i]**4))/(rho[i] * c[i]*V[i]) if isoT[i] == 0 else 0
                else:
                    dTdt[i] = q_flux[i] + Eg[i] + (-k[i-1]*Acond[i-1]*(T[i] - T[i-1]) / (rho[i] * c[i] * dx[i-1] * V[i])
                            - k[i]*Acond[i]*(T[i] - T[i+1]) / (rho[i] * c[i] * dx[i] * V[i]))+(-h[i]*Aconv[i]*(T[i] - Tinf[i])- sig*e[i]*Arad[i]*(T[i]**4 - Tsurr[i]**4))/(rho[i] * c[i]*V[i]) if isoT[i] == 0 else 0
            return dTdt   
    return solve_ivp(func, t_span=t_span, t_eval=t_eval, y0=T, args=(rho,c,q_flux,Eg,Tinf,h,e,k,dx,Aconv,V,Acond,Arad,isoT))


# return solve_ivp(func2, t_span = t_span, t_eval=t_eval, y0=T, args=()  )

# def ThermalSolve():
    print("--------------------------\n1-D TRANSIENT THERMAL SOLVER\nBy: Aiden Ryan\nDescription: \nSolver uses lumped capacitance model and scipy's solve_ivp package to output Acond Temperature vs Time Plot \nof N number of thermal nodes. \n---\nREQUIRED PACKAGES: \n-numpy\n-scipy\n-matplotlib ")
    print("---\nOTHER REQUIREMENTS: \nEnsure Biot Number for each node is < 0.1 \nBi = h*L_c/k \n")
    # x = input("Proceed with user input data? Y/N: \n")
    x = 'Y'
    if x == 'N' or x=='n':
        n = 2
        Acond = 0.1
        dx = 3
        kcond = 45
        nodes = np.empty(n, dtype=object)
        paths = np.empty(n-1, dtype=object)
        paths = Path(Acond, 45,dx)
        t = [0,360]
        t_eval = np.linspace(t[0], t[1], t[1] - t[0])
        nMatrix = nodeMatrix(nodes)
        pMatrix = pathMatrix(paths)
        Biot = nMatrix[6,:]*pMatrix[2]/pMatrix[1]
        for i in range(len(Biot)):
            if Biot[i] > 1:
                print("Conditions do not satisfy Lumped Capacitance Criterion.")
                return 0
    else:
        n = int(input("Define Number of Nodes: "))
        nodes= np.empty(n,dtype=object)
        paths = np.empty(n-1, dtype=object)
        print("ENTER NODE ATTRIBUTES. \n")
        if n == 1:
                T = float(input("ENTER TEMPERATURE OF LUMPED NODE: "))
                mat = input("ENTER MATERIAL: ")
                h = float(input("ENTER HEAT TRANSFER COEFFICIENT: " ))
                V = float(input("ENTER VOLUME OF LUMPED NODE: "))
                Aconv = float(input("ENTER CONVECTIVE AREA OF LUMPED NODE: "))
                Arad = float(input("ENTER RADIATIVE AREA OF LUMPED NODE: "))
                q_flux = float(input("ENTER HEATFLUX THROUGH LUMPED NODE: "))
                Eg = float(input("ENTER HEAT GENERATED FROM LUMPED NODE: "))
                Tinf = float(input("ENTER TEMPERATURE OF FREESTREAM: "))
                Tsurr = float(input("ENTER TEMPERATURE OF SURROUNDINGS: "))
                e = float(input("ENTER EMISSIVITY OF LUMPED NODE: "))
                nodes[0] = Node(T,mat, Aconv, Arad, V, q_flux, Eg, e, h , Tinf, Tsurr)
                paths = [0,0,0]
        else:
            for i in range(n):
                T = float(input("ENTER TEMPERATURE OF LUMPED NODE " +str(i+1)+": "))
                mat = input("ENTER MATERIAL OF LUMPED NODE " +str(i+1)+": ")
                h = float(input("ENTER HEAT TRANSFER COEFFICIENT " +str(i+1)+": "))
                Aconv = float(input("ENTER CONVECTIVE AREA OF LUMPED NODE " +str(i+1)+": "))
                V = float(input("ENTER VOLUME OF LUMPED NODE " +str(i+1)+": "))
                Arad = float(input("ENTER RADIATIVE AREA OF LUMPED NODE " +str(i+1)+": "))
                q_flux = float(input("ENTER HEATFLUX THROUGH LUMPED NODE " +str(i+1)+": "))
                Eg = float(input("ENTER HEAT GENERATED FROM LUMPED NODE " +str(i+1)+": "))
                Tinf = float(input("ENTER TEMPERATURE OF FREESTREAM " +str(i+1)+": "))
                Tsurr = float(input("ENTER TEMPERATURE OF SURROUNDINGS: "))
                e = float(input("ENTER EMISSIVITY OF LUMPED NODE " +str(i+1)+": "))
                nodes[i] = Node(T,mat, Aconv, Arad, V, q_flux, Eg, e, h , Tinf, Tsurr)
        
        if n == 2:
            '''SINGLE PATH'''
            print("ENTER PATH ATTRIBUTES. ENTER 0 TO IGNORE PARAMETER. ")
            Acond = float(input("ENTER CONDUCTION PATH AREA: "))
            k = float(input("ENTER CONDUCTION COEFFICIENT: "))
            dx = float(input("ENTER DISTANCE BETWEEN NODES: "))
            print('\n')
            paths= Path(Acond, k, dx)
        elif n > 2:
            '''MORE THAN ONE PATH'''
            print("ENTER PATH ATTRIBUTES. ENTER 0 TO IGNORE PARAMETER. ")
            paths = np.empty(n-1, dtype=object)
            for i in range(n-1):
                Acond = float(input("ENTER CONDUCTION PATH AREA: "))
                k = float(input("ENTER CONDUCTION COEFFICIENT: "))
                dx = float(input("ENTER DISTANCE BETWEEN NODES: "))
                paths[i] = Path(Acond, k, dx)
                print('\n')

        #sets up nodes for Tvst function
        if n >1:
            nodeMatrix = nodeMatrix(nodes)
            pathMatrix = pathMatrix(paths, n-1)
        else:
            nodeMatrix = nodeMatrix(nodes)
            pathMatrix = paths
        '''
        Check if nodes satisfy Lumped Capacitance and Convergence Criterion
        Biot = hL/k
        Fo = k*dt/(rho*cp*dx**2)
        ''' 
        #Set time condition
        t0 = int(input("ENTER TIME START: "))
        tf = int(input("ENTER TIME STOP: "))
        t = [t0,tf]
        t_eval = np.linspace(t[0], t[1], (tf-t0)*1000)
    
    #Solve
    sol = T_vs_t(t, t_eval = t_eval,nodeMatrix=nodeMatrix, pathMatrix=pathMatrix)
    n = len(nodeMatrix[0])
    plt.figure(1)
    for i in range(n):
        plt.plot(sol.t, sol.y[i],label=str(sol.y[i,0]))
        plt.xlabel('Time (s)')
        plt.ylabel('Temperature (K)')
        plt.title('Temperature Node vs Time Plot')
        plt.legend()
    plt.show()
    y = input("Would you like to change your time window? Y/N: ")
    if y == "Y" or y == "y":
        while y != "N": 
            t0 = int(input("ENTER TIME START: "))
            tf = int(input("ENTER TIME STOP: "))
            t = [t0,tf]
            t_eval = np.linspace(t[0], t[1], int((tf-t0)/.001))
            sol = T_vs_t(t, t_eval = t_eval,nodeMatrix=nodeMatrix, pathMatrix=pathMatrix)
            for i in range(n):
                plt.plot(sol.t, sol.y[i],label=str(sol.y[i,0]))
                plt.xlabel('Time (s)')
                plt.ylabel('Temperature (K)')
                plt.title('Temperature Node vs Time Plot')
                plt.legend()
            plt.show()
            y = input("Would you like to change your time window? Y/N: ")
    print("goodbye.")
# ThermalSolve()
