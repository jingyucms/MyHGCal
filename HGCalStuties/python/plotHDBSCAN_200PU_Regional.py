import math, sys, ROOT
import numpy as np

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
import matplotlib.axes as axs
#from mpl_toolkits.mplot3d import Axes3D

events = np.linspace(1,100,100)

inputFileDir='../test/'

inputFileName="h_step2_eta2p5_eta2p5_pt10_pt20_p211_PU200_100.root"

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

def calc_mean_phi(x):
    if max(x)-min(x)>6 and max(x)*min(x)<0:
        mask = x>0
        xpos=x[mask]
        xneg=x[~mask]
        xneg+=2*math.pi
        xnew=np.concatenate((xpos, xneg))
        mean_phi=np.mean(xnew)
        if mean_phi>math.pi:
            mean_phi-=2*math.pi
    else:
        mean_phi=np.mean(x)

    return mean_phi

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

scTree=fil.Get("hgcalHitNtuple/SimHitsInSimClusterTree")
scs=scTree.GetBranch("simhitsInSimClusterInfo")
SCs = np.asarray([[scs.GetLeaf("event").GetValue(),
                   scs.GetLeaf("x").GetValue(),
                   scs.GetLeaf("y").GetValue(),
                   scs.GetLeaf("z").GetValue(),
                   scs.GetLeaf("energy").GetValue(),
                   scs.GetLeaf("eta").GetValue(),
                   scs.GetLeaf("phi").GetValue(),
                   scs.GetLeaf("cpidx").GetValue(),
                   scs.GetLeaf("scidx").GetValue()] for sc in scTree])


isRecHits=True
isAllLC=True
whichLabel=5
whichSide='positive' ### Don't change this

if isRecHits:
    rechitTree=fil.Get("hgcalHitNtuple/RecHitSignalRegionTree")
    rechits=rechitTree.GetBranch("RecHitsSignalRegionInfo")

    RHs = np.asarray([[rechits.GetLeaf("event").GetValue(),
                       rechits.GetLeaf("x").GetValue(), 
                       rechits.GetLeaf("y").GetValue(), 
                       rechits.GetLeaf("z").GetValue(), 
                       rechits.GetLeaf("energy").GetValue(), 
                       rechits.GetLeaf("eta").GetValue(), 
                       rechits.GetLeaf("phi").GetValue(),
                       rechits.GetLeaf("thickness").GetValue(), 
                       rechits.GetLeaf("layer").GetValue()] for rechit in rechitTree])
else:
    lcTree=fil.Get("hgcalHitNtuple/LCTree")
    lcInfo=lcTree.GetBranch("lcInfo")
    LCs = np.asarray([[lcInfo.GetLeaf("event").GetValue(),
                       lcInfo.GetLeaf("x").GetValue(),
                       lcInfo.GetLeaf("y").GetValue(),
                       lcInfo.GetLeaf("z").GetValue(),
                       lcInfo.GetLeaf("energy").GetValue(),
                       lcInfo.GetLeaf("eta").GetValue(),
                       lcInfo.GetLeaf("phi").GetValue(),
                       lcInfo.GetLeaf("size").GetValue(),
                       lcInfo.GetLeaf("layer").GetValue()] for lc in lcTree])

responses=[]
responses_rec=[]
for event in events:
    print("Start to process event:",event)
    #if not event in [52]: continue

    if isRecHits:
        rechits_in_event=RHs[RHs[:,0]==event]
        mask = (rechits_in_event[:, 3]>0) & (2.2<rechits_in_event[:,5]) & (rechits_in_event[:,5]<2.8) & (-0.3<rechits_in_event[:,6]) & (rechits_in_event[:,6]<0.3)
        RH=rechits_in_event[mask]
        outfilename="out_hdbscan_200PU_pt10_pt20_allRecHits_"+whichSide+"_evt"+str(int(event))+"_v3.npz"
        data = np.load(outfilename, allow_pickle=True)
        Labels=data['Labels'+str(whichLabel)]
        #if len(set(Labels))<10:
        #    Labels=data['Labels'+str(whichLabel+1)]
        labels=Labels
        sample=RH
    else:
        lcs_in_event=LCs[LCs[:,0]==event]
        if isAllLC:
            lcmask = (lcs_in_event[:, 3]>0) & (2.2<lcs_in_event[:,5]) & (lcs_in_event[:,5]<2.8) & (-0.5<lcs_in_event[:,6]) & (lcs_in_event[:,6]<0.5)
            outfilename="out_hdbscan_200PU_pt10_pt20_evt"+str(int(event))+"_"+whichSide+"_v3.npz"
        else:
            lcmask = (lcs_in_event[:, 3]>0) & (lcs_in_event[:, 7]>=2) & (2.2<lcs_in_event[:,5]) & (lcs_in_event[:,5]<2.8) & (-0.5<lcs_in_event[:,6]) & (lcs_in_event[:,6]<0.5)
            outfilename="out_hdbscan_200PU_pt10_pt20_nonMipLCs_evt"+str(int(event))+"_"+whichSide+"_v3.npz"
        LC=lcs_in_event[lcmask]
        data = np.load(outfilename, allow_pickle=True)
        if whichLabel==0:
            Labels=data['Labels']
        else:
            Labels=data['Labels'+str(whichLabel)]
        labels=Labels
        sample=LC
        
    cluster_labels=set(labels[labels[:]>=0])

    if not len(sample)==len(labels):
        print("------ERROR------!!!")

    sc_in_event=SCs[SCs[:,0]==event]

    for icp in [0,1]:
        cpmask = CPs[:,0]==event
        CP = CPs[cpmask]
        cp=CP[CP[:,1]==icp][0]
        cp_energy=cp[2]
        cp_rec_energy=cp[3]
        cp_eta=cp[4]
        cp_phi=cp[5]
        if whichSide=='positive' and cp_eta<0:
            continue

        sc_in_cp=sc_in_event[sc_in_event[:,7]==icp]
        scList=set(sc_in_cp[:,8])
        sc_energy_cut=20
        cl_energy_cut=2
        distance_cut=0.1
        for isc in scList:
            sc_mask = (sc_in_cp[:,8]==isc)
            sc=sc_in_cp[sc_mask]
            sc_energy=sum(sc[:,4])
            sc_rec_energy=sum(sc[:,4])
            #print(isc, sc_energy, sc_eta, sc_phi)
            if sc_energy <= sc_energy_cut: continue
            sc_phi_mask =  (-0.3<sc[:,6]) & (sc[:,6]<0.3)
            sc_phi_masked=sc[sc_phi_mask]
            sc_eta=np.mean(sc_phi_masked[:,5])
            sc_phi=calc_mean_phi(sc_phi_masked[:,6])
            
           
            
            distances=[]
            for cl in cluster_labels:
                hits_in_cluster = sample[labels==cl]
                cl_energy = sum(hits_in_cluster[:, 4])
                if cl_energy <= cl_energy_cut: continue
                cl_eta = np.mean(hits_in_cluster[:,5])
                cl_phi = calc_mean_phi(hits_in_cluster[:,6])
                #print(sc_eta, sc_phi, cl_eta, cl_phi)
                distances+=[[distance(sc_eta, sc_phi, cl_eta, cl_phi), cl, cl_energy]]

            if len(distances)==0: continue
            distances=np.array(distances)
            #print(distances)
            if not len(distances)==0:
                if min(distances[:,0]) > distance_cut:
                    response=0
                    response_rec=0
                else:
                    matched=distances[distances[:,0]==min(distances[:,0])]
                    response=matched[0][2]/sc_energy
                    response_rec=matched[0][2]/sc_rec_energy
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

            print(response)

bins=np.linspace(0, 20, 21)*0.1

presponse=plt.figure(1)
#ax=axs.Axes(presponse)
plt.hist(responses, bins, histtype='step')
plt.hist(responses_rec, bins, histtype='step', fill=False)
plt.xlim(0,2)
if isRecHits:
    saveName="response_recHits_200PU_v3_"+str(whichLabel)+".pdf"
else:
    if isAllLC:
        saveName="response_lcs_200PU_v3_"+str(whichLabel)+".pdf"
    else:
        saveName="response_nonMipLCs_200PU_v3_"+str(whichLabel)+".pdf"
plt.savefig(saveName)
presponse.clear()
