
'''
Assumptions: 
ignore radiation 

'''
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import mat_lib as p


material_mapping = {
    'SS316': "p.SS316",
    "SS316": "p.SS316",

}
class Node:
    '''
    Creates Temperature node with HT attributes
    '''
    # def __init__(self, T, kcond, Ac=0, Ar = 0,density = 803, V=.01, c=1, q=0, Eg=0, h=0,Tinf=298,Tsurr=298,e = 0, mat = 'SS316'):
    def __init__(self, T, mat=p.SS316,  Ac=0, Ar=0, V= 0.01, q=0, Eg = 0, e = 0,h= 0,Tinf=298,Tsurr=298): 


        #316 is default material
        self.mat = eval(material_mapping.get(mat))
        self.V = V #m3
        self.T = T #K 
        self.q = q #W
        self.Eg = Eg #W
        self.h = h #W/m^2K
        self.Tinf = Tinf #K
        self.Tsurr = Tsurr #K
        self.e = e # dimensionless, emissivity of body
        self.Ac = Ac #Convection Area
        self.Ar = Ar #Radiation Area
        self.density = p.getMatProp(self.mat, "density")#kg/m3 
        self.k= p.getMatProp(self.mat, "k") #W/m-K
        self.c = p.getMatProp(self.mat, "specificHeatCapacity") #J/kgK
        if (self.h*self.V/(Ac*self.k) > 0.1 ):
                raise Exception("Biot Number is greater than 0.1, does not satisfy Lumped Capacitance Criterion.")
        
    def get_args(self):
        return [self.density, self.c, self.q, self.Eg, self.h, self.Tinf, self.Tsurr, self.e,self.Ac,self.Ar,self.V]
    
class Path:
    '''
    Creates Conduction Pathway between Nodes
    Assumptions: Path is of Constant Area
    '''
    def __init__(self, A1, k=0, dx=1):
        self.A1 = A1 #Conduction Area [m^2]
        self.k = k #conduction coefficient [W/mK]
        self.dx = dx #length [m]
    def get_args(self):
        return [self.A1, self.k, self.dx]

def TMatrix(nodes):
    '''
    Creates Temperature attribute matrix to be passed into T_vs_t function
    '''
    array1 = np.array([node.T for node in nodes], dtype= float)
    array2 = np.array([node.density for node in nodes],dtype= float)
    array3 = np.array([node.c for node in nodes],dtype= float)
    array4 = np.array([node.q for node in nodes],dtype= float)
    array5 = np.array([node.Eg for node in nodes],dtype= float)
    array6 = np.array([node.Tinf for node in nodes],dtype= float)
    array7 = np.array([node.h for node in nodes],dtype= float)
    array8 = np.array([node.e for node in nodes],dtype= float)
    array9 = np.array([node.Tsurr for node in nodes],dtype= float)
    array10 = np.array([node.Ac for node in nodes],dtype= float)
    array11 = np.array([node.Ar for node in nodes],dtype= float)
    array12 = np.array([node.V for node in nodes],dtype= float)
    array13 = np.array([node.k for node in nodes], dtype=float)
    return np.vstack((array1, array2, array3, array4, array5,array6, array7, array8, array9,array10,array11, array12,array13))

def pMatrix(paths, n):
    '''
    Creates path attribute matrix to be passed into T_vs_t function
    '''
    if n > 1:
        array1 = [path.A1 for path in paths]
        array2 = [path.k for path in paths]
        array3 = [path.dx for path in paths]
        return np.vstack((array1, array2, array3))
    else:
        return paths.get_args()

def T_vs_t(t_span, t_eval, nodeMatrix, pathMatrix):
    '''
    Iteratively solves Thermal Network. Returns Temperature vs time plot.
    '''
    n = len(nodeMatrix[0, :])
    T_array = nodeMatrix[0,:]
    rho = nodeMatrix[1,:]
    c = nodeMatrix[2,:]
    q = nodeMatrix[3,:]
    Eg = nodeMatrix[4,:]
    Tinf = nodeMatrix[5,:]
    h = nodeMatrix[6,:]
    e = nodeMatrix[7,:]
    Tsurr = nodeMatrix[8,:]
    Ac = nodeMatrix[9,:]
    Ar = nodeMatrix[10,:]
    V = nodeMatrix[11,:]
    kcond = nodeMatrix[12,:]

    if isinstance(pathMatrix, list) == True:
        '''
        If only one path, assign path attributes to respective list indices 
        '''
        A = pathMatrix[0]
        k = pathMatrix[1]
        dx = pathMatrix[2]
    else:
        '''
        if more than one path exists
        '''
        A = pathMatrix[0,:]
        k = pathMatrix[1,:]
        dx = pathMatrix[2,:]

    dTdt = np.zeros(n,dtype=float)
    def func(t, T_array,rho,c,q,Eg,Tinf,h,e,A,k,dx,Ac,Ar,V):
        sig = 5.67e-8 #Stefan-Boltzmann Constant
        if n == 2:
            # Two nodes and one path
            dTdt[0] = (q[0] +Eg[0] - k * A * (T_array[0] - T_array[1]) /(dx) 
            - h[0]*Ac[0]*(T_array[0] - Tinf[0]) - sig*e[0]*Ar[0]*(T_array[0]**4 - Tsurr[0]**4))/(rho[0] *c[0]*V[0])
            dTdt[1] = (q[1] +Eg[1] + k * A * (T_array[0] - T_array[1]) / (dx) 
            - h[1]*Ac[1]*(T_array[1] - Tinf[1]) - sig*e[1]*Ar[1]*(T_array[1]**4 - Tsurr[1]**4))/(rho[1] *c[1]*V[1])
            return dTdt       
        elif n == 1:
            return q+Eg- h*Ac*(T_array - Tinf) - sig*e*Ar*(T_array**4 - Tsurr**4)/(rho *c*V)
        else:
            for i in range(n):
                if i == 0:
                    dTdt[i] = (q[i] + Eg[i] -k[i]*A[i]*(T_array[i] - T_array[i+1])/dx[i] 
                            -h[i]*Ac[i]*(T_array[i] - Tinf[i])- sig*e[i]*Ar[i]*(T_array[i]**4 - Tsurr[i]**4))/(rho[i] * c[i]*V[i])
                elif i == n-1:
                    dTdt[i] = (q[i] + Eg[i] -k[i-1]*A[i-1]*(T_array[i] - T_array[i-1])/dx[i-1] 
                            -h[i]*Ac[i]*(T_array[i] - Tinf[i])- sig*e[i]*Ar[i]*(T_array[i]**4 - Tsurr[i]**4))/(rho[i] * c[i]*V[i])
                else:
                    dTdt[i] = q[i] + Eg[i] + (-k[i-1]*A[i-1]*(T_array[i] - T_array[i-1]) / (rho[i] * c[i] * dx[i-1] * V[i])
                            - k[i]*A[i]*(T_array[i] - T_array[i+1]) / (rho[i] * c[i] * dx[i] * V[i]))+(-h[i]*Ac[i]*(T_array[i] - Tinf[i])- sig*e[i]*Ar[i]*(T_array[i]**4 - Tsurr[i]**4))/(rho[i] * c[i]*V[i])
            return dTdt
    
    return solve_ivp(func, t_span=t_span, t_eval=t_eval, y0=T_array, args=(rho,c,q,Eg,Tinf,h,e,A,k,dx,Ac,Ar,V))
    

def ThermalSolve():
    print("--------------------------\n1-D TRANSIENT THERMAL SOLVER\nBy: Aiden Ryan\nDescription: \nSolver uses lumped capacitance model and scipy's solve_ivp package to output a Temperature vs Time Plot \nof N number of thermal nodes. \n---\nREQUIRED PACKAGES: \n-numpy\n-scipy\n-matplotlib ")
    print("---\nOTHER REQUIREMENTS: \nEnsure Biot Number for each node is < 0.1 \nBi = h*L_c/k \n")
    x = input("Proceed with user input data? Y/N: \n")
    if x == 'N' or x=='n':
        n = 2
        A = 0.1
        dx = 3
        kcond = 45
        nodes = np.empty(n, dtype=object)
        paths = np.empty(n-1, dtype=object)
        nodes[0] = Node(500, kcond,  A, A, c = 500, h = 10)
        nodes[1] = Node(200, kcond, A, A, c= 500, h = 10)
        paths = Path(A, 45,dx)
        t = [0,360]
        t_eval = np.linspace(t[0], t[1], t[1] - t[0])
        nodeMatrix = TMatrix(nodes)
        pathMatrix = pMatrix(paths, n-1)
        Biot = nodeMatrix[6,:]*pathMatrix[2]/pathMatrix[1]
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
                Ac = float(input("ENTER CONVECTIVE AREA OF LUMPED NODE: "))
                Ar = float(input("ENTER RADIATIVE AREA OF LUMPED NODE: "))
                q = float(input("ENTER HEATFLUX THROUGH LUMPED NODE: "))
                Eg = float(input("ENTER HEAT GENERATED FROM LUMPED NODE: "))
                Tinf = float(input("ENTER TEMPERATURE OF FREESTREAM: "))
                Tsurr = float(input("ENTER TEMPERATURE OF SURROUNDINGS: "))
                e = float(input("ENTER EMISSIVITY OF LUMPED NODE: "))
                nodes[0] = Node(T,mat, Ac, Ar, V, q, Eg, e, h , Tinf, Tsurr)
                paths = [0,0,0]
        else:
            for i in range(n):
                T = float(input("ENTER TEMPERATURE OF LUMPED NODE " +str(i+1)+": "))
                kcond = float(input("ENTER THERMAL CONDUCTIVITY OF NODE " +str(i+1)+": "))
                h = float(input("ENTER HEAT TRANSFER COEFFICIENT " +str(i+1)+": "))
                Ac = float(input("ENTER CONVECTIVE AREA OF LUMPED NODE " +str(i+1)+": "))
                V = float(input("ENTER VOLUME OF LUMPED NODE " +str(i+1)+": "))
                if (h*V/(Ac*kcond) > 0.1 ):
                    raise Exception("Biot Number is greater than 0.1, does not satisfy Lumped Capacitance Criterion.")
                Ar = float(input("ENTER RADIATIVE AREA OF LUMPED NODE " +str(i+1)+": "))
                rho = float(input("ENTER DENSITY OF LUMPED NODE " +str(i+1)+": "))
                c = float(input("ENTER SPECIFIC HEAT CAPACITY OF LUMPED NODE " +str(i+1)+": "))
                q = float(input("ENTER HEATFLUX THROUGH LUMPED NODE " +str(i+1)+": "))
                Eg = float(input("ENTER HEAT GENERATED FROM LUMPED NODE " +str(i+1)+": "))
                Tinf = float(input("ENTER TEMPERATURE OF FREESTREAM " +str(i+1)+": "))
                Tsurr = float(input("ENTER TEMPERATURE OF SURROUNDINGS: "))
                e = float(input("ENTER EMISSIVITY OF LUMPED NODE " +str(i+1)+": "))
                mat = input("ENTER MATERIAL OF LUMPED NODE " +str(i+1)+": ")
                print('\n')
                nodes[i] = Node(T,kcond, Ac, Ar, rho, V, c, q, Eg, h, Tinf, Tsurr, e, mat)
        
        if n == 2:
                print("ENTER PATH ATTRIBUTES. ENTER 0 TO IGNORE PARAMETER. ")
                A = float(input("ENTER CONDUCTION PATH AREA: "))
                k = float(input("ENTER CONDUCTION COEFFICIENT: "))
                dx = float(input("ENTER DISTANCE BETWEEN NODES: "))
                print('\n')
                paths= Path(A, k, dx)
        elif n > 2:
            print("ENTER PATH ATTRIBUTES. ENTER 0 TO IGNORE PARAMETER. ")
            for i in range(n-1):
                A = float(input("ENTER CONDUCTION PATH AREA: "))
                k = float(input("ENTER CONDUCTION COEFFICIENT: "))
                dx = float(input("ENTER DISTANCE BETWEEN NODES: "))
                paths[i] = Path(A, k, dx)
                print('\n')

        #sets up nodes for Tvst function
        if n >1:
            nodeMatrix = TMatrix(nodes)
            pathMatrix = pMatrix(paths, n-1)
        else:
            nodeMatrix = TMatrix(nodes)
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
