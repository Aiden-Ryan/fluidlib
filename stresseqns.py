def shearstress(P, A):
    return P/A

def maxDeflectionIntLoad(P, a, L, E, I):
    return (P*a**2)*(3*L - a)/(6*E*I) 
    
def areaMomentofInertiaRect(b, h):
    return b*h**3/12