'''
Assumptions: 
ignore radiation 
constant heat transfer coefficient

'''
import numpy as np
import CoolProp.CoolProp as CP

class node:
    def __init__(n, T, rho=0, V=0, c=0, q=0, E_g=0):
        #Build a Thermal node object that has attributes T, rho, V, C, q, and E_g
        n.rho = rho #kg/m3
        n.T = T #K
        n.V = V #m3
        n.c = c #J/kgK
        n.q = q #W
        n.E_g = E_g #W
    def getT(n):
        return n.T
    def getArgs(n):
        return [n.rho, n.V, n.c, n.q, n.E_g]
    
class nodeArray:
    def __init__(self) -> None:
        self.array = None
    def build_array(self, *nodes):
        params_list = [node.T for node in nodes]
        self.array = np.array(params_list)
    def get_array(self):
        return self.array

class path: 
    def __init__(p, A1, Tinf = 298, h=0, e=0, k=0, dx=0):
        #Builds a thermal path object that has attributes A, Tinf, h, e, k, and dx
        p.A1 = A1
        p.h = h
        p.e = e
        p.k =k
        p.dx = dx
        p.Tinf = Tinf
    def getArgs(p): 
        return [p.Tinf, p.h, p.A1, p.e, p.k, p.dx]
        
class buildNetwork(node, path):
    
    def sol(nodeArray, path):
        pathArgs = path.getArgs(node)
        return[
            '''-k*A*dT/dx'''
            -pathArgs[4]*pathArgs[2]*nodeArray[0]-nodeArray[1]/pathArgs[6]
        ]




def ThermalSolve(h,dx,k,alpha,t,T,Tinf=298):
    Bi = h*dx/k
    Fo = alpha*t/dx**2
    if Fo*(1 + Bi) > .5: 
        print(Fo*(1 + Bi)) 
        return("Does not meet stability criterion") 

    jsize = len(T[:,0]-1)
    for i in range(t-1):
        for j in range(jsize):
            if j == 0: 
                T[j,i+1] = 2*Fo*(T[j+1, i] + Bi*Tinf) + (1-2*Fo - 2*Bi*Fo)*T[j,i]
            elif j == (jsize-1):
                T[j,i+1] = 2*Fo*(T[j-1, i] + Bi*Tinf) + (1-2*Fo - 2*Bi*Fo)*T[j,i]
            else:
                T[j, i+1] = Fo*(T[j+1, i] + T[j-1, i]) + (1-2*Fo)*T[j,i]
            
    return T





    


    