import numpy as np
import math

def etaToR(eta, Z):
    return Z*math.tan(2*math.atan(pow(math.e, -1*eta)))

particle=[22]  # pion
nparticles=1
eMin=100.
eMax=100.
etaMin=1.6
etaMax=1.6
#phiMin=3.1415926
#phiMax=3.1415926
phiMin=0
phiMax=0
#phiMin=0
#phiMax=0
zMin=410
zMax=410
#rMin=etaToR(etaMax, zMin)
#Max=etaToR(etaMin, zMax)
rMin=160.
rMax=160.
#print rMin, rMax
#print etaToR(etaMin, 500), etaToR(etaMax, 500)
#delta=15.0
#delta=12.5
#delta=10.0
delta=7.5
isOverlapping=False

if nparticles>1:
    if (abs(phiMin)==3.1415926 and abs(phiMax)==3.1415926) or (phiMin==0 and phiMax==0):
        filename="step1_r"+str(int(rMin))+"_r"+str(int(rMax))+"_e"+str(int(eMin))+"_p"+str(particle[0])+"_delta"+str(delta).replace(".","p")+"_phiPi.root"
        nevt=1
    else:
        filename="step1_r"+str(int(rMin))+"_r"+str(int(rMax))+"_e"+str(int(eMin))+"_p"+str(particle[0])+"_delta"+str(delta).replace(".","p")+".root"
        nevt=1000
else:
    if (abs(phiMin)==3.1415926 and abs(phiMax)==3.1415926) or (phiMin==0 and phiMax==0):
        filename="step1_r"+str(int(rMin))+"_r"+str(int(rMax))+"_e"+str(int(eMin))+"_p"+str(particle[0])+"_phiPi.root"
        nevt=1
    else: 
        filename="step1_r"+str(int(rMin))+"_r"+str(int(rMax))+"_e"+str(int(eMin))+"_p"+str(particle[0])+".root"
        nevt=1000
