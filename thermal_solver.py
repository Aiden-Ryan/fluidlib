'''
Assumptions: 
ignore radiation 
constant heat transfer coefficient

'''
def transientODE(t, T, h, A_surface, rho, V, c, q, E_g, Tinf):
    return q*A_surface + E_g - (h*A_surface)/(rho*V*c)*(T - Tinf) 

    