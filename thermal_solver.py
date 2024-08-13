'''
Assumptions: 
ignore radiation 
constant heat transfer coefficient

'''
import numpy as np
   

class node:
    def __init__(n, T, rho=0, V=0, c=0, q=0, E_g=0):
        n.rho = rho
        n.T = T
        n.V = V
        n.c = c 
        n.q = q
        n.E_g = E_g
    def getT(n):
        return n.T
    def getArgs(n):
        return [n.rho, n.V, n.c, n.q, n.E_g]


class path: 
    def __init__(p, A1, A2=0, Tinf = 298, h=0, e=0, k=0, L=0):
        p.A1 = A1
        if A2 == 0:
            p.A2 = A1
        p.h = h
        p.e = e
        p.k =k
        p.L = L
        p.Tinf = Tinf
    def getArgs(p): 
        return [p.Tinf, p.h, p.A1, p.A2, p.e, p.k, p.L]
        

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





    


    