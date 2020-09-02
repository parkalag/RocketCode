# RocketEquations
# Written by Brendan

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ipywidgets as widgets
import time
import math
from IPython.display import display

def PropType(Ox,Fuel):

    # O/Fm, O/Fv, SG, T, c*, MW, Isp, k
    
    # Oxygen Based
    if Ox=='Oxygen' and Fuel=='Methane':
        C=[3.20,1.19,0.81,3526,1835,20.3,296,1.2]
    elif Ox=='Oxygen' and Fuel=='Hydrazine':
        C=[0.74,0.66,1.06,3285,1871,18.3,301,1.25]
    elif Ox=='Oxygen' and Fuel=='Hydrogen':
        C=[3.40,0.21,0.26,2959,2428,8.9,386,1.26]
    elif Ox=='Oxygen' and Fuel=='RP-1':
        C=[2.24,1.59,1.01,3571,1774,21.9,285.4,1.24]
    elif Ox=='Oxygen' and Fuel=='UDMH':
        C=[1.39,0.96,0.96,3542,1835,19.8,295,1.25]

    # Flourine Based
    elif Ox=='Flourine' and Fuel=='Hydrazine':
        C=[2.30,1.54,1.31,4713,2208,19.4,365,1.33]
    elif Ox=='Flourine' and Fuel=='Hydrogen':
        C=[4.54,0.21,0.33,3080,2534,8.9,389,1.33]

    # Nitrogen Tetroxide Based    
    elif Ox=='Nitrogen Tetroxide' and Fuel=='Hydrazine':
        C=[1.08,0.75,1.20,3258,1765,19.5,283,1.26]
    elif Ox=='Nitrogen Tetroxide' and Fuel=='UDMH/Hydrazine':
        C=[1.62,1.01,1.18,3242,1652,21.0,278,1.24]
    elif Ox=='Nitrogen Tetroxide' and Fuel=='RP-1':
        C=[3.4,1.05,1.23,3290,'Unknown',24.1,297,1.23]
    elif Ox=='Nitrogen Tetroxide' and Fuel=='MMH':
        C=[1.65,1.00,1.16,3200,1591,21.7,278,1.23]

    # Nitric Acid Based
    elif Ox=='Nitric Acid' and Fuel=='RP-1':
        C=[4.1,2.12,1.35,3175,1594,24.6,258,1.22]
    elif Ox=='Nitric Acid' and Fuel=='UDMH/Hydrazine':
        C=[1.73,1.00,1.23,2997,1682,20.6,272,1.22]

    # Hydrogen Peroxide Based
    elif Ox=='Hydrogen Peroxide' and Fuel=='RP-1':
        C=[7.0,4.01,1.29,2760,'Unknown',21.7,297,1.19]

    # Nitrous Oxide Based
    elif Ox=='Nitrous Oxide' and Fuel=='ABS':
        C=[5.0,10.0,'SG',3100.0,1560.0,25.6,190.0,1.15]

    # Monoprop, Solid, and Cold Gas
    elif Ox=='N/A' and Fuel=='Air':
        C=[1.0,1.0,1.3,273.15,408.904,28.9647,55.0,1.4]
    elif Ox=='N/A' and Fuel=='Sorbitol':
        C=[1.857,1.0,'SG',1600,909,39.86,164,1.1362]

    # Error
    else:
        raise ValueError('Incompatible or Incorrect Propellant Combination Chosen')

    return C

def LARA(L15,AR):       # Length and Arera Ratio Give Angles

    if L15>=0.9:
        if 1<AR<=4:
            C=[18,10]
        elif 4<AR<=5:
            C=[19,9]
        elif 5<AR<=6:
            C=[20,8]
        elif 6<AR<=8:
            C=[21,7]
        elif 8<AR<=10:
            C=[22,6]
        elif 10<AR<=20:
            C=[24,6]
        elif 20<AR<=30:
            C=[26,5]
        elif 30<AR<=40:
            C=[28,5]
        elif 40<AR<=50:
            C=[29,5]
        elif 50<AR<=100:
            C=[31,5]
        else:
            C=[33,4]
    elif L15>=0.8:
        if 1<AR<=4:
            C=[21,14]
        elif 4<AR<=5:
            C=[22,13]
        elif 5<AR<=6:
            C=[24,13]
        elif 6<AR<=8:
            C=[25,11]
        elif 8<AR<=10:
            C=[26,11]
        elif 10<AR<=20:
            C=[28,10]
        elif 20<AR<=30:
            C=[29,9]
        elif 30<AR<=40:
            C=[30,8]
        elif 40<AR<=50:
            C=[31,8]
        elif 50<AR<=100:
            C=[33,7]
        else:
            C=[34,7]
    elif L15>=0.7:
        if 1<AR<=4:
            C=[23,18]
        elif 4<AR<=5:
            C=[24,17]
        elif 5<AR<=6:
            C=[26,16]
        elif 6<AR<=8:
            C=[27,14]
        elif 8<AR<=10:
            C=[28,14]
        elif 10<AR<=20:
            C=[30,13]
        elif 20<AR<=30:
            C=[32,12]
        elif 30<AR<=40:
            C=[33,11]
        elif 40<AR<=50:
            C=[34,11]
        elif 50<AR<=100:
            C=[35,10]
        else:
            C=[37,9]
    else:
        raise ValueError('Data for this range is not available.')
        
    return C


def ChamberP(pambient,k,M):
    po=pambient*(1+0.5*(k-1)*M**2)**(k/(k-1))
    return po

def ExpansionRatio(k,M):
    ER=(1/M)*((1+((k-1)/2)*M**2)/((k+1)/2))**((k+1)/(2*(k-1)))
    return ER

def AreaRatio(k,p2p1):
    AR=((k+1)/2)**(1/(k-1))*(p2p1)**(1/k)*(((k+1)/(k-1))*(1-p2p1**((k-1)/k)))**0.5
    AR=AR**-1
    return AR

def ThroatArea(mdot,p1,R,T1,k):
    At=(mdot/p1)*((R*T1)/(k*(2/(k+1))**((k+1)/(k-1))))**0.5
    return At

def PressureRatio(k,M):
    PPo=((1+0.5*(k-1)*M**2)**(k/(k-1)))**(-1)
    return PPo

def TempRatio(k,M):
    TTo=((1+0.5*(k-1)*M**2))**(-1)
    return TTo

def PropMass(InitialM , FinalM):
    m=InitialM-FinalM
    return m

def MRVehicle(InitialM , FinalM):
    mr=FinalM/InitialM
    return mr

def MRSubSys(InitialM , FinalM , NonSysMass):
    mr=(FinalM-NonSysMass)/(InitialM-NonSysMass)
    return mr

def PMFVehicle(PropellantMass , InitialM):
    PMF=PropellantMass/InitialM
    return PMF

def PMFSubSys(PropellantMass , InitialM , NonSysMass):
    PMF=PropellantMass/(InitialM-NonSysMass)
    return PMF

def MassFlow(ThrustOrPropM , EVelOrBurnTime):
    mflow=ThrustOrPropM/EVelOrBurnTime
    return mflow

def WeightFlow(MFR , gravitySL):  
    WFR=MFR*gravitySL  
    return WFR

def EEV_ISP(Isp , Gravity):
    c=Isp*Gravity
    return c

def EEV_Force(Thrust , MFR):    
    c=Thrust/MFR    
    return c

def IdealThrust(C , MFR):
    F=MFR*C
    return F

def MomThrust(MFR , V2):    
    mt=MFR*V2    
    return mt

def RealThrust(MomentumThrust , PressureThrust):    
    T=MomentumThrust+PressureThrust    
    return T

def TotalImpulse(Force , Time):    
    I=Force*Time   
    return I

def Acceleration(Thrust , Mass):   
    a=Thrust/Mass    
    return a

def CStar(ChamberP , ThroatArea , MFR):    
    cstar=(ChamberP*ThroatArea)/MFR    
    return cstar

def SpecificImpulseT(Thrust , WeightFlowRate):    
    Isp=Thrust/WeightFlowRate    
    return Isp

def SpecificImpulseV(V , SpecGrav):    
    Isp=V/SpecGrav    
    return Isp

def PressThrust(p2 , p3 , A2):    
    pt=(p2-p3)*A2
    return pt

def ExhaustVmass(TotalThrust , PThrust , MFR):
    v=(TotalThrust-PThrust)/MFR
    return v

def ExhaustVMet(k,R,T1,p2p1):
    v=(((2*k)/(k-1))*R*T1*(1-(p2p1)**((k-1)/k)))**(1/2)
    return v

def ExhaustVUS(k,R,t1,MW,PR):
    v=(((2*k*R*t1*32.2)/((k-1)*MW))*(1-(PR**((k-1)/k))))**0.5
    return v

def ExhaustVMol(k,R,t1,MW,PR):
    v=(((2*k*R*t1)/((k-1)*MW))*(1-(PR**((k-1)/k))))**0.5
    return v

def ExhaustVMAX(R,k,t1):
    v=((2*R*t1*k)/(k-1))**0.5
    return v

def ThroatV(k,R,ToM):
    v=(2*k*R*ToM/(k+1))**0.5
    return v

def Sonic(k,R,t):
    a=(k*R*t)**0.5
    return a

def CfIdeal(k,p2p1):
    Cf=(((2*k**2)/(k-1))*(2/(k+1))**((k+1)/(k-1))*(1-(p2p1)**((k-1)/k)))**0.5
    return Cf

def T2Ideal(t1,p1,p2,k):
    t2=t1*(p1/p2)**((1-k)/k)
    return t2

def ExitP(p1,k,Me):
    p=p1*(1+0.5*(k-1)*Me**2)**(-k/(k-1))
    return p

def FindEpsln(numbers,eps):   
     numbers = np.asarray(numbers) 
     i=(np.abs(numbers - eps)).argmin() 
     return numbers[i]

def TotalWdot(F,Isp):
    w=F/Isp
    return w

def OxWdot(Wdot,r):
    wo=Wdot*r/(r+1)
    return wo

def FuelWdot(Wdot,r):
    wf=Wdot/(r+1)
    return wf

