import numpy as np
import math

def etaToR(eta, Z):
    return Z*math.tan(2*math.atan(pow(math.e, -1*eta)))

# To get the variables used in CloseByParticleGunProducer see: https://github.com/cms-sw/cmssw/blob/master/IOMC/ParticleGuns/src/CloseByParticleGunProducer.cc

#fromVtx=True
fromVtx=False

#isD41=True
isD41=False

particle=[211]  # pion
nparticles=1
ptMin=20   #used in FlatRandomPtGunProducer
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
eMin=eMax=50
#zMin=410     # front face of HGC scint
#zMax=410    
#rMin=160         
#rMax=160
zMin=300     # ~front of HGC
zMax=300
rMin=40
rMax=60
#rMin=rMax=50
#delta=15.0
#delta=12.5
delta=10.0
#delta=7.5
isOverlapping=False

nevt=1

#print etaToR(2.5, 300), etaToR(2.2, 300), etaToR(2.7, 300)

if not fromVtx:
    if nparticles<2:
        filename="step1_r"+str(rMin).replace(".","p")+"_r"+str(rMax).replace(".","p")+"_e"+str(eMin).replace(".","p")+"_p"+str(particle[0])+".root"
        if isD41:
            filename = filename.replace(".root","_D41PU200.root")
    else:
        filename="step1_r"+str(rMin).replace(".","p")+"_r"+str(rMax).replace(".","p")+"_e"+str(eMin).replace(".","p")+"_p"+str(particle[0])+"_delta"+str(delta).replace(".","p")+".root"
else:
    filename="step1_eta"+str(etaMin).replace(".","p")+"_eta"+str(etaMin).replace(".","p")+"_pt"+str(ptMin).replace(".","p")+"_p"+str(particle[0])+".root"
