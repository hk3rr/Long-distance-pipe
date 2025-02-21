# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 16:48:07 2021

@author: hk3rr
"""
import numpy as np
import math as math

#Physical properties

rho = 1000
#Density of the fluid in kg/m^3
mu = 1*10**-3
#Dynamic viscosity of the fluid in Pas
g = 9.81
#Acceleration experienced due to gravity in ms^-2

# The system properties     

D1 = 0.05
D2 = 0.1
#Diameters in m
A1 = math.pi*(D1/2)**2
A2 = math.pi*(D2/2)**2
#Areas of pipes in m 
L1 = 100
L2 = 2300
#Lengths in m
H1 = 55
H2 = 49
HLtot = H1-H2
#Heads in m
Kabs = 1E-5
#Absolute roughness in m
Krel1 = Kabs/D1
Krel2 = Kabs/D2
#Relative roughness (unitless)

#Function definitions
def swamja(Re,Kabs,D):
    return (1/(-2*(math.log(((Kabs/(3.7*D))+(5.7/(Re**0.9))),10))))**2

def darcy(HL,L,lamb,g,D):
    return ((HL*2*g*D)/L*lamb)**0.5

def darcyandbalancegivingV(g,HLtot,lamb1,lamb2,L1,L2,A1,A2):
    return ((HLtot*2*g)/(((lamb1*L1)/D1)+((lamb2*L2*(A1**2/A2**2))/D2)))**0.5

def Reeq(rho,V,D,mu):
    return (rho*V*D)/mu

#Initial guess of the Reynolds Number
Re1 = 1E5
Re2 = 1E5

for i in range (10):
    lamb1 = swamja(Re1,Kabs,D1)
    lamb2 = swamja(Re2,Kabs,D1)
    V1 = darcyandbalancegivingV(g,HLtot,lamb1,lamb2,L1,L2,A1,A2)
    V2 = V1 * (D1**2 / D2**2)
    Re1 = Reeq(rho,V1,D1,mu)
    Re2 = Reeq(rho,V2,D2,mu)

#Check that the Reynlds number for the last few iterations hasn't changed
#The final flowrate

Q1 = V1*(np.pi*((D1/2)**2))
#Flow rate through pipe 1 in m^3s^-1
Q2 = V2*(np.pi*((D2/2)**2))
#Flow rate through pipe 2 in m^3s^-1

L1 = Q1*1000
L2 = Q2*1000
print('V1 =',V1)
print('V2 =',V2)

print('lamb1 =',lamb1)
print('lamb2 =',lamb2)

print('Q1 =',Q1)
print('Q2 =',Q2)

print('Flow in Litres/s',L1)
