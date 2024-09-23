import conversion_library as c 
import numpy as np
'''NUSSELT NUMBER CALCULATOR: SOLVE FOR NUSSELT NUMBER USING SIMPLIFIED GEOMETRIES'''

def getRe(u, fluid, D, P2prop, P2, Lc):
    mu = c.getFluidProperty(hOut = 'V', hIn='D' + P2prop, val1 = D, val2=P2, hFLd=fluid)
    return D*u*Lc/mu

    

    

def flatPlateNu(Re, Pr,type='AVERAGE'):
    '''TYPE - LOCAL OR AVERAGE'''
    if Pr < 0.6: 
        raise Exception("Pr must be greater than or equal to 0.6")
    
    if Re < 2300:
        '''LAMINAR FLOW'''
        if type == 'LOCAL':
            return 0.332*Re**0.5*Pr**(1/3)
        elif type == 'AVERAGE':
            return 0.664*Re**0.5*Pr**(1/3)
    elif Re > 2300 and Re < 4000:
        '''TRANSIENT'''
    elif Re > 4000:
        '''TURBULENT'''
        if type == 'LOCAL':
            return 0.0296*Re**(4/5)*Pr**(1/3)
        
