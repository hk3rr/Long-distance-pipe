# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 19:28:14 2021

@author: hk3rr
"""
Dst = [0, 32,59,1093,1727,2116,2170,2214,2574,4668,5262,5644,6017,6062,7493,8171,11287,13042,13119,13172,14207,14437,14461,14680,15412,15926,16145,17055,17670,18379,18553,19467,19915,20332,21350,22216,22265,22665,22748,22902,23838,25852,26221,26734,27770,29125,30343,30846,31433,31567,32320,34168,35364,35987,36985,38750,39077,39330,39511,40605,41172,41702,42131,42591,43230,44884,47027,47459,47642,48350,48613,49187,49463,51365,52647,53310,53870,54668,55000,55075,55554,56019,56637,57360,57416,57476,57926,57926,57970,59425,61220,61448,61633,61976,62772,63360,63551,63705,63914,64151,64273,64695,65218,65526,66230,66655,67053,67255,67379,67853,67961,68258,68763,69442,69700,70031,70741,71030,71257,71399,71605,72028,73195,73395,73522,73700,75070,75341,76096,76288,76664,76939,77139,77370,77545,77791,78050,78296,78545,78755,78900,79133,79275,79718,80059,80345,80539,80686,80968,81704,82491,83353,84000]

Elv = [139,135,135,132,128,128,128,128,127,116,112,107,106,106,119,122,95,87,87,87,91,92,88,88,80,76,65,45,62,52,42,37,38,48,73,80,80,82,83,83,82,83,79,75,80,77,79,79,78,74,71,75,76,76,75,76,74,74,74,71,69,67,65,71,74,78,75,71,71,69,68,73,73,67,65,66,67,59,60,61,63,60,67,67,67,67,65,65,65,61,55,50,50,48,51,50,49,49,49,46,46,46,14,17,2,5,10,10,10,29,29,35,33,32,32,35,34,33,33,36,37,38,38,37,37,40,48,53,59,62,63,66,68,67,67,66,67,37,37,25,25,19,18,12,7,10,7,7,8,6,8,6,8]

import numpy as np
import math as math
import matplotlib.pyplot as plt

#Physical properties

rho = 1000
#Density of the fluid in kg/m^3
mu = 1*10**-3
#Dynamic viscosity of the fluid in Pas
g = 9.81
#Acceleration experienced due to gravity in ms^-2

Kabs = 1E-5
#Absolute roughness in m
#Krel1 = Kabs/D1
#Krel2 = Kabs/D2
#Relative roughness (unitless)

#Function definitions
def swamja(Re,Kabs,D):
    return (1/(-2*(math.log(((Kabs/(3.7*D))+(5.7/(Re**0.9))),10))))**2

def darcy(HL,L,lamb,g,D):
    return ((HL*2*g*D)/L*lamb)**0.5

def swamjaforD(Kabs,lamb,Re):
    return (Kabs*(10**(1/(2*np.sqrt(lamb)))-((Re**0.9)/5.74)))/3.7

def Reeq(rho,V,D,mu):
    return (rho*V*D)/mu

def horrorfuncforV(lamb1,lamb2,L1,L2,A1,A2,D1,D2,QD,HLtot):
    a = lamb1*L1*A2**2*D2 
    a1 = D1*lamb2*L2*A1**2
    
    b = D1*lamb2*L2*A1*QD
    
    a2 = D1*A2**2*D2
    b2 = 2*HLtot*lamb1*L1*g*A2**2*D2 
    c2 = 2*HLtot*D1*lamb2*g*L2*A1**2 
    d2 = lamb1*L1*lamb2*L2*QD**2
    
    c = a2*(b2+c2-d2)
    
    return (b+c**0.5)/(a+a1)

def TotalEnergyLine(lamb,L,D,V,g):
    return (lamb*L*V**2)/(D*2*g)

def HydraulicGradeLine(lamb,L,D,V,g):
    return TotalEnergyLine(lamb,L,D,V,g) - (V**2)/2*g

def HEGradient(lamb,D,V,g):
    return (lamb*V**2)/(D*2*g)

def trenchCSA(D):
    if D<0.08:
        TW = 0.3
        TD = D+1
    else:
        if D>0.6:
            TW = D+0.6
            TD = D+1
        else:
            TW = D+0.3
            TD = D+1
    return TW*TD

HeadChange = []
for i in range(1,len(Elv)):
    HeadChange.append(Elv[i-1]-Elv[i])

HydDist = []
for i in range(1,len(Elv)):
    HydDist.append((((Dst[i-1]-Dst[i])**2) + (HeadChange[i-1]**2))**0.5)

shortDst = Dst[1:]
HGLgrad = np.divide(HeadChange,shortDst)

#Time to draw HGLs for every pipe

HGLDst = [Dst[0],Dst[15],Dst[79],Dst[138],Dst[152]]
HGLElv = [Elv[0],Elv[15],Elv[79]+20,Elv[136]+3,Elv[152]]

plt.figure(0)
plt.plot(Dst,Elv)
plt.plot(HGLDst,HGLElv)
plt.xlabel("Distance (m)")
plt.ylabel("Elevation (m)")

HGLLength = []
HGLLoss = []
for i in range(1,len(HGLDst)):
    HGLLength.append(HGLDst[i]-HGLDst[i-1])
    HGLLoss.append(HGLElv[i-1]-HGLElv[i])

Dratios = []
for i in range(1,len(HGLLength)):
    Dratios.append(round(((HGLLoss[i-1]/HGLLength[i-1])/(HGLLoss[i]/HGLLength[i]))**0.2,1))

for i in range(0,len(Dratios)):
    print('D',[i+2],'=','D',[i+1],'*',Dratios[i])

RAFLs = 10
RAF = RAFLs/1000

Deez = []
Queues = []
for D1 in np.arange(0.01,1,0.01):
  Deez.append(D1)
  D2 = D1*Dratios[0]
  D3 = D2*Dratios[1]
  D4 = D3*Dratios[2]

  Diam = [D1,D2,D3,D4]

  Area = []
  for i in range(0,len(Diam)):
     Area.append(math.pi*(Diam[i]/2)**2)

  HLtot = []
  for i  in range(1,len(HGLLoss)):
     HLtot.append(HGLLoss[i-1]+HGLLoss[i])

#Initil guess of the Reynolds Number
  Re1 = 1E5
  Re2 = 1E5

  Q = []
  Qb = []
  for j in range(1,len(Diam)):
     if j == 2:
        QD = RAF
     else:
        QD = 0
     for i in range (20):
        lamb1 = swamja(Re1,Kabs,Diam[j-1])
        lamb2 = swamja(Re2,Kabs,Diam[j])
   
        V1 = horrorfuncforV(lamb1,lamb2,HGLLength[j-1],HGLLength[j],Area[j-1],Area[j],Diam[j-1],Diam[j],QD,HLtot[j-1])
    
        V2 = V1 * Diam[j-1]**2 / Diam[j]**2
        Re1 = Reeq(rho,V1,Diam[j-1],mu)
        Re2 = Reeq(rho,V2,Diam[j],mu)
 
        Q1 = V1*(np.pi*((Diam[j-1]/2)**2))
        Q2 = V2*(np.pi*((Diam[j]/2)**2))
    
     Q.append(Q1)
     Qb.append(Q2)
  Queues.append(np.mean(Q))

plt.figure(1)
plt.plot(Deez,Queues)
plt.grid(True)
plt.xlabel("Diameter of pipe 1 (m)")
plt.ylabel("Flow rate for pipes (m^3s^-1)")

minElvs = [12,31,104,144]

for i in range(0,len(Diam)):
    grad = (HGLElv[i+1]-HGLElv[i])/(HGLDst[i+1]-HGLDst[i])
    c = HGLElv[i] - grad*HGLDst[i]
 
    y = grad*Dst[minElvs[i]] + c
 
    Pmaxhead = y - Elv[minElvs[i]]

    Pmaxbar = Pmaxhead*0.0981
    print('Max pressure in pipe',i,'=',Pmaxbar)

#Now to use stupid monkey brain 
#Be on the conservative side and give us some pipes that easily 
#fill ratio requirement and also withstand the pressure
#And you'll have to do the kg/m because THERE'S NO GODDAMN EQ
#EQUATION WHICH GIVES WEIGHT FROM DIAMETER AND RATING GODDAMN

D1 = 0.450
PRating1 = 2.5
W1 = 15.9

D2 = 0.56
PRating2 = 10
W2 = 84

D3 = 0.630
PRating3 = 10
W3 =  106

D4 = 0.315
PRating4 = 6
W4 = 17

Pipe1elms = HydDist[:16]
Pipe2elms = HydDist[15:80]
Pipe3elms = HydDist[79:139]
Pipe4elms = HydDist[138:]
kg = W1 * sum(Pipe1elms)
CO21 = 2.52*kg
kg = W2 * sum(Pipe2elms)
CO22 = 2.52*kg
kg = W3 * sum(Pipe3elms)
CO23 = 2.52*kg
kg = W4 * sum(Pipe4elms)
CO24 = 2.52*kg

TA1 = trenchCSA(D1)*sum(Pipe1elms) 
TCO21 = TA1*4.18 
TA2 = trenchCSA(D2)*sum(Pipe2elms) 
TCO22 = TA2*4.18 
TA3 = trenchCSA(D3)*sum(Pipe3elms) 
TCO23 = TA3*4.18 
TA4 = trenchCSA(D4)*sum(Pipe4elms) 
TCO24 = TA4*4.18 

TotalCO2 = CO21+CO22+CO23+CO24+TCO21+TCO22+TCO23+TCO24
print('Total CO2 produced in kg =',TotalCO2)
print('')
 








