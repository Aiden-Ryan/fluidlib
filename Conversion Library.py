import numpy as np
from CoolProp.CoolProp import PropsSI 
import math as m

# Globals 
g = 9.81 #m/s^2
g_EE = 32.2 #ft/s^2

# Getters
def getR(cp, cv):
    return

# Cv Conversions 
def Cv2Do(Cv):
    Do = 0.236 * Cv**0.5
    return Do

def Cv2F(Cv): 
    F = 0.556*Cv
    return F

def Cv2CdA(Cv):
    Do = Cv2Do(Cv)
    CdA = 0.6*m.pi()*Do**2/4
    return CdA

#Do conversions
def Do2Cv(Do):
    Cv = 18.0*pow(Do,2) #power function
    return Cv

def Do2CdA(Do, Cd):
    Cv = Do2Cv(Do)
    Do = Cv2CdA(Cv)
    return Do

def Do2F(Do):
    F = 10*pow(Do,2)
    return F


#Flow Factor Conversions
def F2Do(F): 
    Do = 0.316*pow(F, 0.5)
    return Do

def F2Cv(F): 
    Cv = 1.8*F
    return Cv


#Equivalent Orifice Diameter
def Do(Q, w, dP):
    Do = 1.445*pow(Q,0.5)*pow(w,0.25)/pow(dP,0.25)
    return Do

# Flow Equations
def mdotIG(gamma, P1, P2, T1, CdA, M):
    R = getR(cp, cv)
    if M == 1:
        mdot = CdA*P1*m.sqrt(g*gamma/(R*T1)*(2/(gamma + 1))**((gamma+1)/gamma-1))
    else:
        mdot = CdA*(P2/P1)**(1/gamma)*m.sqrt(2*g*R*T1*gamma/(gamma-1)*(1-(P2/P1))^(gamma-1)/gamma)
    return mdot