import conversion_library as c 
'''NUSSELT NUMBER CALCULATOR: SOLVE FOR NUSSELT NUMBER USING SIMPLIFIED GEOMETRIES'''

def getRe(u, fluid, D, Pressure, Lc):
    '''Returns Reynolds number, requires pressure and density input
    Lc == 4*Ac/P'''
    Viscosity = c.getFluidProperty(hOut = 'V', hIn='DP', val1 = D, val2=Pressure, hFLd=fluid)
    return D*u*Lc/Viscosity
    
def flatPlateNu(Re, Pr,type='AVERAGE'):
    '''Returns Nusselt Number for a flat plate for a provided Reynolds number and Prandtl number'''
    '''TYPE - LOCAL OR AVERAGE'''
    if Pr < 0.6: 
        raise Exception("Pr must be greater than or equal to 0.6")
    if Re < 5e5:
        '''LAMINAR FLOW'''
        if type == 'LOCAL':
            return 0.332*Re**0.5*Pr**(1/3)
        elif type == 'AVERAGE':
            return 0.664*Re**0.5*Pr**(1/3)
    elif Re == 5e5:
        '''TRANSIENT'''
    elif Re > 5e5 and Re < 1e7:
        '''TURBULENT'''
        if type == 'LOCAL' and (Pr > 0.6 and Pr < 60):
            return 0.0296*Re**(4/5)*Pr**(1/3)
        
def internalFlowNu(Re, Pr=0.6, BoundaryCondition = 'constantQ', f = 1):
    '''Returns Nusselt Number for circular tube: 
    Boundary Conditions:
    constantQ;
    constantTs'''
    if Re < 2300 and BoundaryCondition == 'constantQ':
        '''LAMINAR FLOW fully developed and constant surface heat flux'''
        return 4.36
    elif Re < 2300 and BoundaryCondition == 'constantTs': 
        return 3.66
    elif Re > 10^4 and Re < 5e6 and Pr > 0.5 and Pr < 2000:
        '''Applies to constantQ and constantTs'''
        return (f/8)*(Re*Pr)/(1.07 + 12.7*(f/8)**0.5*(Pr**(2/3)-1))
