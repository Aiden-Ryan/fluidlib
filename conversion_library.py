import sys

import numpy as np
import CoolProp.CoolProp as CP
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import math as m
import os

check = 0
try:
    RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
    RP.SETPATHdll(os.environ['RPPREFIX'])
    MOLAR_BASE_SI = RP.GETENUMdll(0,"MOLAR BASE SI").iEnum
except:
    print("REFPROP unavailable, defaulting to CoolProp")
    check = 1

# Getters
def getR(fluid):
    ''' 
    Returns Specific Gas Constant for provided fluid in [J/(kg*K)]

    ''' 
    M = CP.PropsSI('M', fluid)
    R = 8.314/M
    return R


def getFluidProperty(hOut, hIn, val1, val2, hFLd): 
    '''
    Returns thermodynamic property for two provided properties and specified fluid
    hFLd - fluid
    hIn ~ "P1P2"
    hOut ~ "P3"
    '''
    if check == 0:
        prop3 = RP.REFPROPdll(hFLd,hIn,hOut, MOLAR_BASE_SI,0,0,val1, val2, [1.0]) # Don't know what the two zeros are for exactly but 
    else:
        prop1 = hIn[0]
        prop2 = hIn[1]
        prop3 = CP.PropsSI(hOut, prop1, val1, prop2, val2, hFLd)
    return prop3


# Cv Conversions 
def Cv2Do(Cv):
    '''
    Returns equivalent ESEOD (Cd = 0.6) for a provided Cv
    Cv: Flow coefficient defined as GPM H2O given 1psid across component
    '''
    return 0.236*Cv**0.5

def Cv2F(Cv): 
    '''
    Returns flow factor for a provided Cv
    Cv: Flow coefficient defined as GPM H2O given 1psid across component at 60F
    ''' 
    return 0.556*Cv 

def Cv2CdA(Cv):
    '''
    Returns effective flow Area for a provided Cv.
    Assume Cd = 0.6
    Cv: Flow coefficient defined as GPM H2O given 1psid across component
    '''
    Do = Cv2Do(Cv)
    return 0.6*m.pi*Do**2/4

#Do conversions
def Do2Cv(Do):
    '''
    Returns CV for a provided ESEOD
    Do: equivalent sharp-edged orifice diameter (ESEOD)

    '''
    return (Do/.236)**2 


def Do2CdA(Do):
    '''
    Returns effective flow area for a provided Do
    Do: equivalent sharp-edged orifice diameter (ESEOD)
    '''
    Cv = Do2Cv(Do)
    return Cv2CdA(Cv)
    

def Do2F(Do):
    '''
    Returns Flow Factor for a provided ESOD (Cd = 0.6)
    
    '''
    return 10.01*Do**2
    


#Flow Factor Conversions
def F2Do(F): 
    '''
    Returns ESOD for a provided flow factor

    '''
    return 0.316*F**.5


def F2Cv(F): 
    '''
    Returns Cv for a provided flow Factor

    '''
    return 1.8*F


#Equivalent Orifice Diameter
def Do(Q, w, dP):
    '''
    Returns ESEOD for a provided flow rate Q, specific weight w, and Pressure drop dP
    '''
    return 1.445*Q**2*w**.25/dP**.25
    

# Flow Equations

def mdotIG(P1, P2, T1, CdA, fluid):
    '''
    Returns mass flow rate for an ideal gas for a provided pressure differential, Temperature, CdA, and Fluid

    '''
    g = 9.81 #m/s^2
    R = getR(fluid)
    cp = getFluidProperty(fluid, 'PT', 'Cp0mass', P1, T1) #Need to be able to distinguish between Coolprop and REFPROP
    gamma = cp/(cp - R)
    PcOverP1 = (2/(gamma+1))**(gamma/(gamma-1))

    if P2/P1 <= PcOverP1:
        '''CHOKED'''
        return CdA*P1*((g*gamma/(R*T1))*(2/(gamma+1))**((gamma+1)/(gamma-1)))**.5
    else:
        return CdA*(P2/P1)**(1/gamma)*(2*g*R*T1*(gamma/(gamma-1))*(1-(P2/P1)**((gamma-1)/gamma)))**.5

# Mass Conversions
def kg2lb(kg):
    '''
    kilogram to pound
    '''
    return 2.204*kg
    

def lb2kg(lb):
    '''
    Pound to Kilogram
    '''
    return lb/2.204
 
# Length Conversions
def m2ft(m):
    '''
    meter to foot
    '''
    return m/.3048
   

def ft2m(ft):
    '''
    Foot to Meter
    '''
    return ft*.3048
  

def inch2m(inch):
    '''
    Inch to Meter 
    '''
    return inch*.0254
    

def m2inch(m):
    ''' 
    Meter to Inch
    '''
    return m/.0254
    
# Area Conversions

def sm2sf(sm):  
    '''
    Meter squared to foot squared
    '''
    return sm*3.281**2

def sf2sm(sf): 
    '''
    Foot squared to meter squared
    '''
    return sf/(3.281**2)

def si2sm(si):
    '''
    inches squared to meters squared
    '''
    return si*.0254**2
    

def sm2si(sm): 
    '''
    Meters squared to inches squared 
    '''
    return sm/(.0254**2)
    

# Volume Conversions

def cintocm(cin): 
    '''
    cubic inch to cubic meter
    '''
    return cin/61023.7441
    

def cmtocin(cm):
    '''
    cubic meter to cubic inch
    '''
    return cm*61023.7441
   

# Density Conversions

def kgcmtolbcf(kgcm):
    '''
    
    kilogram per cubic meter to pound per cubic foot
    '''
    return kgcm/16.018


def lbcftokgcm(lbcf):
    '''
    pound per cubic foot to kilogram per cubic meter
    '''
    return lbcf*16.018
    

# Temperature Coversions
def F2K(F):
    '''
    Fahrenheit to Kelvin 
    '''
    return (F-32)/1.8+273.15

def K2F(K): 
    '''
    Kelvin to Fahrenheit
    '''
    return ((K-273.15)*1.8)+32


def C2F(C): 
    '''
    Celsius to Fahrenheit
    '''
    return C*1.8+32

def F2C(F):
    '''
    Fahrenheit to Celsius
    '''
    return (F-32)/1.8

# Pressure Conversions

def psi2pa(psi):
    '''
    Psi to Pascal
    '''
    return psi*6894.76

def pa2psi(pa): 
    '''
    Pascal to Psi
    '''
    return pa/6894.76

def cmm2cin(cmm):
    '''
    Cubic millimeters to Cubic inches
    '''
    return cmm/16387.064

