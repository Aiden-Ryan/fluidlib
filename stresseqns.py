def shearstress(P, A):
    return P/A

def maxDeflectionIntLoad(P, a, L, E, I):
    return (P*a**2)*(3*L - a)/(6*E*I) 