# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 18:06:30 2021

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

def horrorfunc(lamb1,lamb2,L1,L2,A1,A2,D1,D2,QD,HLtot):
    a = lamb1*L1*A2**2*D2 
    a1 = D1*lamb2*L2*A1**2
    
    b = D1*lamb2*L2*A1*QD
    
    a2 = D1*A2**2*D2
    b2 = 2*HLtot*lamb1*L1*g*A2**2*D2 
    c2 = 2*HLtot*D1*lamb2*g*L2*A1**2 
    d2 = lamb1*L1*lamb2*L2*QD**2
    
    c = a2*(b2+c2-d2)
    
    return (b+c**0.5)/(a+a1)

LQD = 1.1
#Flow out of the system inbetween the two pipes in ls^-1

QD = LQD / 1000

#Initil guess of the Reynolds Number
Re1 = 1E5
Re2 = 1E5


for i in range (10):
    lamb1 = swamja(Re1,Kabs,D1)
    lamb2 = swamja(Re2,Kabs,D2)
    
    V1 = horrorfunc(lamb1,lamb2,L1,L2,A1,A2,D1,D2,QD,HLtot)
    
    V2 = V1 * D1**2 / D2**2
    Re1 = Reeq(rho,V1,D1,mu)
    Re2 = Reeq(rho,V2,D2,mu)
    
    #print(Re1)
#Check that the Reynlds number f o r the last few iterations hasn't changed
#The final flowrate

Q1 = V1*(np.pi*((D1/2)**2))
#Flow rate through pipe 1 in m^3s^-1
Q2 = Q1 - QD
#Flow rate through pipe 2 in m^3s^-1

print('V1 =',V1)
print('V2 =',V2)

print('lamb1 =',lamb1)
print('lamb2 =',lamb2)

print('Q1 =',Q1)
print('Q2 =',Q2)

HL1 = HLtot - ((lamb2*L2*V2**2)/(2*g*D2))

Hmid = H1 - HL1
print('Total head at the junction = ',Hmid)

L1 = Q1*1000
L2 = Q2*1000
print('Flow in pipe 1 in Litres/s',L1)
print('Flow in pipe 2 in Litres/s',L2)
