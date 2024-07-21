import numpy as np
import CoolProp.CoolProp as CP
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import math as m
import os

g = 9.81 #m/s^2_
g_EE = 32.2 #ft/s^2


# Getters
def getR():
    fluid = input("input desired fluid: ")
    M = CP.PropsSI("M", fluid)
    R = 8.314/M
    return R

def getFluidProperty(hFLd, hIn, hOut, val1, val2): # hIn~"P1P2"; hOut~"P3"
    try:
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        RP.SETPATHdll(os.environ['RPPREFIX'])
        MOLAR_BASE_SI = RP.GETENUMdll(0,"MOLAR BASE SI").iEnum
        prop3 = RP.REFPROPdll(hFLd,hIn,hOut, MOLAR_BASE_SI,0,0,val1, val2, [1.0]) # Don't know what the two zeros are for exactly but 
    except:
        print("REFPROP unavailable, defaulting to CoolProp")
        prop1 = hIn[0]
        prop2 = hIn[1]
        prop3 = CP.PropsSI(hOut, prop1, val1, prop2, val2, hFLd)
    return prop3


# Cv Conversions 
def Cv2Do():
    Cv = float(input("input Cv: "))
    Do = 0.236*m.sqrt(Cv)
    return Do

def Cv2F(): 
    Cv = float(input("input Cv: "))
    F = 0.556*Cv
    return F

def Cv2CdA():
    Cv = float(input("input Cv: "))
    Do = Cv2Do(Cv)
    CdA = 0.6*m.pi*Do**2/4
    return CdA

#Do conversions
def Do2Cv():
    Do = float(input("input Do: "))
    Cv = (Do/.236)**2 #power function
    return Cv

def Do2CdA():
    Do = float(input("input Do: "))
    Cv = float(input("input Cv: "))
    Cv = Do2Cv(Do)
    CdA = Cv2CdA(Cv)
    return CdA

def Do2F():
    Do = float(input("input Do: "))
    F = 10.01*Do**2
    return F


#Flow Factor Conversions
def F2Do(): 
    F = float(input("input F: "))
    Do = 0.316*m.sqrt(F)
    return Do

def F2Cv(): 
    F = float(input("input F: "))
    Cv = 1.8*F 
    return Cv 
 
 
#Equivalent Orifice Diameter 
def Do(Q, w, dP):
    Do = 1.445*m.sqrt(Q)*pow(w,0.25)/pow(dP,0.25)
    return Do

# Flow Equations
def mdotIG(gamma, P1, P2, T1, CdA, M, fluid):
    R = getR(fluid)
    if M == 1:
        mdot = CdA*P1*m.sqrt((g*gamma/(R*T1))*(2/(gamma+1))**((gamma+1)/(gamma-1)))
    else:
        mdot = CdA*(P2/P1)**(1/gamma)*m.sqrt(2*g*R*T1*(gamma/(gamma-1))*(1-(P2/P1)**((gamma-1)/gamma)))
    return mdot

# Mass Conversions
def kg2lb(kg):
    lb = 2.204*kg
    return lb

def lb2kg(lb):
    kg = lb/2.204
    return kg

# Length Conversions
def m2ft(m):
    ft = m/.3048
    return ft

def ft2m(ft):
    m = ft*.3048
    return m

def inch2m(inch):
    m = inch*.0254
    return m

def m2inch(m):
    inch = m/.0254
    return inch
    
# Area Conversions
def sm2sf(sm): 
    sf = sm*3.281**2
    return sf

def sf2sm(sf): 
    sm = sf/(3.281**2)
    return sm

def si2sm(si):
    sm = si*.0254**2
    return sm

def sm2si(sm): 
    si = sm/(.0254**2)
    return si

# Volume Conversions
def in3tom3(in3): 
    m3 = in3/61023.7441
    return m3

def m3toin3(m3):
    in3 = m3*61023.7441
    return in3

# Density Conversions
def kgm3tolbf3(kgm3):
    lbf3 = kgm3/16.018
    return lbf3

def lbf3tokgm3(lbf3):
    kgm3 = lbf3*16.018
    return kgm3

# Temperature Coversions
def F2K(F):
    K = (F-32)/1.8+273.15
    return K

def K2F(K): 
    F = ((K-273.15)*1.8)+32
    return F

def C2F(C): 
    F = C*1.8+32
    return F

def F2C(F):
    C = (F-32)/1.8
    return C

def unknown_function():
    print("Unknown function. Please try again.")

# Dictionary
functions = {
    'getR': getR,
    'Cv2Do': Cv2Do, 
    'Cv2F': Cv2F, 
    'Cv2CdA': Cv2CdA, 
    'Do2Cv':Do2Cv,
    'Do2CdA':Do2CdA,
    'Do2F':Do2F,
    'F2Do':F2Do,
    'F2Cv':F2Cv,


}

while True:
    user_input = input("Enter function name or exit: ").strip()
    if user_input == "exit":
        break
    result = functions.get(user_input, unknown_function)()
    if result is not None:
        print("Result:", result)

#Test Cases
# Cv = .6
# D = Cv2Do(Cv)
# F = Cv2F(Cv)
# CdA = Cv2CdA(Cv)
# Cv=Do2Cv(D)
# CdA= Do2CdA(D,Cv)
# F = Do2F(D)
# D = F2Do(F)
# print(D)
# met = 1
# ft = 1
# m = ft2m(ft)
# print(m)
# ft = m2ft(m)
# print(ft)
# sm = 1
# sf = sm2sf(sm)
# print(sf)
# sm = sf2sm(sf)
# print(sm)
# F = 32
# K= F2K(F)
# print(K)
# F = K2F(K)
# print(F)
# C = 0
# F = C2F(C)
# print(F)
# C = F2C(F)
# print(C)
# si = 1
# sm = si2sm(si)
# print(sm)
# si = sm2si(sm)
# print(si)
# in3 = 1
# m3 = in3tom3(in3)
# print(m3)
# in3 = m3toin3(m3)
# print(in3)
# lbf3 = kgm3tolbf3(1)
# print(lbf3)
# kgm3 = lbf3tokgm3(lbf3)
# print(kgm3)
