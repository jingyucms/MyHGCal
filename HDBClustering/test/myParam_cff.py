import numpy as np
import math

def etaToR(eta, Z):
    return Z*math.tan(2*math.atan(pow(math.e, -1*eta)))

# To get the variables used in CloseByParticleGunProducer see: https://github.com/cms-sw/cmssw/blob/master/IOMC/ParticleGuns/src/CloseByParticleGunProducer.cc

#fromVtx=True
fromVtx=False

isPU200=True
#isPU200=False

#isOverlapping=True
isOverlapping=False

nevt=100

#particle=[211, 22]  # pion
particle=[211]
nparticles=1
#ptMin=5
#ptMax=10
ptMin=10   #used in FlatRandomPtGunProducer
ptMax=20
#etaMin=1.6    #used in FlatRandomPtGunProducer
#etaMax=1.6
etaMin=2.5    #used in FlatRandomPtGunProducer
etaMax=2.5
#phiMin=-3.1415926
#phiMax=3.1415926
phiMin=0     
phiMax=0
#eMin=eMax=300    # used in CloseByParticleGunProducer
eMin=eMax=100
#eMin=eMax=50
#eMin=eMax=20
#zMin=410     # front face of HGC scint
#zMax=410    
#rMin=160
#rMax=160
zMin=zMax=320     # ~front of EE 
#rMin=40
#rMax=60
rMin=rMax=50
#rMin=rMax=100
#delta=15.0
#delta=12.5
delta=100.0
#delta=7.5

print(etaToR(1.8, 320), etaToR(2.2, 320), etaToR(2.55, 320), etaToR(2.7, 320))

if not fromVtx:
    if nparticles<2:
        filename="step1_r"+str(rMin).replace(".","p")+"_z"+str(zMin).replace(".","p")+"_e"+str(eMin).replace(".","p")+"_p"+str(particle[0])+".root"
        if isPU200:
            filename = filename.replace(".root","_PU200_"+str(nevt)+".root")
    else:
        if len(set(particle))>1:
            filename="step1_r"+str(rMin).replace(".","p")+"_z"+str(zMin).replace(".","p")+"_e"+str(eMin).replace(".","p")+"_pMix"+"_delta"+str(delta).replace(".","p")+".root"
        else:
            filename="step1_r"+str(rMin).replace(".","p")+"_z"+str(zMin).replace(".","p")+"_e"+str(eMin).replace(".","p")+"_p"+str(particle[0])+"_delta"+str(delta).replace(".","p")+".root"
else:
    filename="step1_eta"+str(etaMin).replace(".","p")+"_eta"+str(etaMin).replace(".","p")+"_pt"+str(ptMin).replace(".","p")+"_pt"+str(ptMax).replace(".","p")+"_p"+str(particle[0])+".root"
    if isPU200:
        filename = filename.replace(".root","_PU200_"+str(nevt)+".root")
