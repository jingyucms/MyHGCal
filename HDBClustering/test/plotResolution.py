import math, sys, ROOT
import numpy as np

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
import matplotlib.axes as axs
#from mpl_toolkits.mplot3d import Axes3D
import networkx as nx

events = np.linspace(1,100,100)

inputFileDir='./'

inputFileName="h_step3_r50_z320_e100_p211_delta100p0.root"

inputFile=inputFileDir+inputFileName

fil=ROOT.TFile(inputFile)

def distance(eta1, phi1, eta2, phi2):
    deta=eta1-eta2
    dphi=phi1-phi2
    if dphi > math.pi:
        dphi -= 2*math.pi
    if dphi < -math.pi:
        dphi += 2*math.pi
    dist_angle=math.sqrt(deta**2+dphi**2)

    return dist_angle

#def calc_mean_phi(x):
#    if max(x)-min(x)>6 and max(x)*min(x)<0:
#        mask = x>0
#        xpos=x[mask]
#        xneg=x[~mask]
#        xneg+=2*math.pi
#        xnew=np.concatenate((xpos, xneg))
#        mean_phi=np.mean(xnew)
#        if mean_phi>math.pi:
#            mean_phi-=2*math.pi
#    else:
#        mean_phi=np.mean(x)
#
#    return mean_phi

def calc_mean_phi(x):
    #mask = (x < 3.14) & (x > 0)
    #x = x[mask]
    return np.mean(x)

cpTree=fil.Get("hgcalHitNtuple/CPTree")
cps=cpTree.GetBranch("cpInfo")
CPs = np.asarray([[cps.GetLeaf("event").GetValue(),
                   cps.GetLeaf("idx").GetValue(),
                   cps.GetLeaf("energy").GetValue(),
                   cps.GetLeaf("energy_rec").GetValue(),
                   cps.GetLeaf("eta").GetValue(),
                   cps.GetLeaf("phi").GetValue(),
                   cps.GetLeaf("id").GetValue(),
                   cps.GetLeaf("pt").GetValue()] for cps in cpTree])


simClusterTree=fil.Get("hgcalHitNtuple/SimHitsInSimClusterTree")
simclusters=simClusterTree.GetBranch("simhitsInSimClusterInfo")
SCs = np.asarray([[simclusters.GetLeaf("event").GetValue(),
                   simclusters.GetLeaf("x").GetValue(), 
                   simclusters.GetLeaf("y").GetValue(), 
                   simclusters.GetLeaf("z").GetValue(), 
                   simclusters.GetLeaf("energy").GetValue(), 
                   simclusters.GetLeaf("eta").GetValue(), 
                   simclusters.GetLeaf("phi").GetValue(),
                   simclusters.GetLeaf("cpidx").GetValue(),
                   simclusters.GetLeaf("scidx").GetValue()] for simcluster in simClusterTree])

recoClusterTree=fil.Get("hgcalClusterNtuple/RecHitsInClusterTree")
recclusters=recoClusterTree.GetBranch("RecHitsInClusterInfo")
RCs = np.asarray([[recclusters.GetLeaf("event").GetValue(),
                   recclusters.GetLeaf("x").GetValue(), 
                   recclusters.GetLeaf("y").GetValue(), 
                   recclusters.GetLeaf("z").GetValue(), 
                   recclusters.GetLeaf("energy").GetValue(), 
                   recclusters.GetLeaf("eta").GetValue(),
                   recclusters.GetLeaf("phi").GetValue(),
                   recclusters.GetLeaf("layer").GetValue(),
                   recclusters.GetLeaf("label").GetValue()] for reccluster in recoClusterTree])


responses=[]
responses_rec=[]
for event in events:
    print("Start to process event:",event)

    sample = RCs[RCs[:,0]==event]
    labels = sample[:, 8]
        
    cluster_labels=set(labels[labels[:]>0]) 

    if not len(sample)==len(labels):
        print("------ERROR------!!!")

    SC = SCs[SCs[:,0]==event]

    cpList = set(SC[:, 7])

    #print(len(cpList))
    #print(CP[:,1])

    #sys.exit()

    distance_cut = 0.1
    for icp in cpList:
        CP = CPs[CPs[:,0]==event]
        #print(CP[:,1])
        cp=CP[CP[:,1]==icp][0]
        cp_energy=cp[2]
        cp_rec_energy=cp[3]
        cp_eta=cp[4]
        cp_phi=cp[5]
        
        distances=[]
        for icl in cluster_labels:
            hits_in_cluster = sample[labels==icl]
            cl_energy = sum(hits_in_cluster[:, 4])
            
            cl_eta = np.mean(hits_in_cluster[:,5])
            cl_phi = calc_mean_phi(hits_in_cluster[:,6])
            #print(hits_in_cluster[:,6])
            #print(cp_eta, cp_phi, cl_eta, cl_phi)
            distances+=[[distance(cp_eta, cp_phi, cl_eta, cl_phi), icl, cl_energy]]

        if len(distances)==0: continue
        distances=np.array(distances)
        #print(distances)
        if not len(distances)==0:
            if min(distances[:,0]) > distance_cut:
                response=0
                response_rec=0
            else:
                matched=distances[distances[:,0]==min(distances[:,0])]
                response=matched[0][2]/cp_energy
                response_rec=matched[0][2]/cp_rec_energy
                cluster_labels.remove(matched[0][1])
                if response > 2:
                    response=2
                    if response_rec > 2:
                        response_rec=2
        else:
            response=0
            response_rec=0

        responses+=[response]
        responses_rec+=[response_rec]

        if response < 0.2 or response > 1.5:
            print(response_rec, "------")
        else:
            print(response_rec)
        

bins=np.linspace(0, 20, 21)*0.1

presponse=plt.figure(1)
#ax=axs.Axes(presponse)
plt.hist(responses, bins, histtype='step')
plt.hist(responses_rec, bins, histtype='step', fill=False)
plt.xlim(0,2)
saveName="response_"+inputFileName.replace(".root", ".pdf")

plt.savefig(saveName)
presponse.clear()
