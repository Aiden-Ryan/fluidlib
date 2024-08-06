'''
Assumptions: 
ignore radiation 
constant heat transfer coefficient

'''

   

class node:
    def __init__(n, T, rho, V, c, q=0, E_g=0):
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
        

def ODE(t, T, Tinf, h, A1, A2, e, k ,L, rho, V, c, q, E_g):

    return q*A1 + E_g - (h*A1)/(rho*V*c)*(T - Tinf)

# def ODE(t, T, path.getArgs()):
#     return T.q*T.A + T.E_g - (path.h*T.A)/(T.rho*T.V*T.c)*(T.T- path.Tinf)


    


    