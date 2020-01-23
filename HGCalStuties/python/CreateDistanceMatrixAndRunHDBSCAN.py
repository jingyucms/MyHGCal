import numpy as np
import ROOT, math, sys

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import hdbscan

def distanceMatrices(x, y):
    dl = abs(x[8]-y[8])
    deta=x[5]-y[5]
    dphi=x[6]-y[6]
    if dphi > math.pi:
        dphi -= 2*math.pi
    if dphi < -math.pi:
        dphi += 2*math.pi
    dist_angle=math.sqrt(deta**2+dphi**2)
    
    w=math.sqrt(x[7])

    if dl>=10:
        distance = 9999
    elif x[8] > 28 and x[8] <= 40:
        distance =  dist_angle/(2*w)
    elif x[8] > 40:
        distance =  dist_angle/(4*w)
    else:
        distance = dist_angle/(w)

    return distance

inputFileDir="/uscms_data/d3/jingyu/HGC/3DClustering/CMSSW_11_0_0_pre6/src/MyHGCal/HGCalStuties/test/"

inputFileName="h_step2_r100_z320_e300_p211_PU200_100.root"

inputFile=inputFileDir+inputFileName

fil=ROOT.TFile(inputFile)

lcTree=fil.Get("hgcalHitNtuple/LCTree")
lcInfo=lcTree.GetBranch("lcInfo")

cpTree=fil.Get("hgcalHitNtuple/CPTree")
cps=cpTree.GetBranch("cpInfo")

whichSide='positive'
isAllLC=True

LCs = np.asarray([[lcInfo.GetLeaf("event").GetValue(),
                   lcInfo.GetLeaf("x").GetValue(),
                   lcInfo.GetLeaf("y").GetValue(),
                   lcInfo.GetLeaf("z").GetValue(),
                   lcInfo.GetLeaf("energy").GetValue(),
                   lcInfo.GetLeaf("eta").GetValue(),
                   lcInfo.GetLeaf("phi").GetValue(),
                   lcInfo.GetLeaf("size").GetValue(),
                   lcInfo.GetLeaf("layer").GetValue()] for lc in lcTree])

CPs = np.asarray([[cps.GetLeaf("event").GetValue(),
                   cps.GetLeaf("idx").GetValue(),
                   cps.GetLeaf("energy").GetValue(),
                   cps.GetLeaf("energy_rec").GetValue(),
                   cps.GetLeaf("eta").GetValue(),
                   cps.GetLeaf("phi").GetValue(),
                   cps.GetLeaf("id").GetValue(),
                   cps.GetLeaf("pt").GetValue()] for cps in cpTree])

for event in range(1, 101):

    print("Starting to Process Event:", event)
    
    lcs_in_event=LCs[LCs[:,0]==event]
    cps_in_event=CPs[CPs[:,0]==event]

    signal_cp1 = cps_in_event[0]
    signal_cp2 = cps_in_event[1]

    if signal_cp1[4]>0:
        signal_cp = signal_cp1
    else:
        signal_cp = signal_cp2

    cp_eta = signal_cp[4]
    cp_phi = signal_cp[5]

    #print(cp_eta, cp_phi)

    #continue

    if isAllLC:
        lcmask = (lcs_in_event[:, 3]>0) & (cp_eta-0.3<lcs_in_event[:,5]) & (lcs_in_event[:,5]<cp_eta+0.3) & (cp_phi-0.3<lcs_in_event[:,6]) & (lcs_in_event[:,6]<cp_phi+0.3)
    else:
        lcmask = (lcs_in_event[:, 3]>0) & (lcs_in_event[:, 7]>=2) & (2.2<lcs_in_event[:,5]) & (lcs_in_event[:,5]<2.8) & (-0.5<lcs_in_event[:,6]) & (lcs_in_event[:,6]<0.5)
    LC=lcs_in_event[lcmask]

    n=len(LC)
    matrix=np.zeros((n,n), dtype=np.float32)

    #continue
    
    for i in range(n):
        for j in range(n):
            x=LC[i]
            y=LC[j]
            matrix[i][j]=distanceMatrices(x, y)

    clusterer = hdbscan.HDBSCAN(min_samples=5, min_cluster_size=10, metric='precomputed').fit(matrix.astype(np.float64))
    labels = clusterer.labels_
    probabilities = clusterer.probabilities_

    nw=clusterer.condensed_tree_.to_networkx()
    g=[]
    for e in nw.edges:
        edge=[e[0], e[1], nw.get_edge_data(*e)['weight']]
        edge=np.array(edge)
        g+=[edge]
        
    g=np.array(g)
    
    d=clusterer.condensed_tree_.to_numpy()
    d=[list(i) for i in d]
    d=np.array(d)

    i=0
    clMap=[]
    for l in labels:
        clMap+=[[i, l]]
        i+=1
    clMap=np.array(clMap)

    if isAllLC:
        outfilename=inputFileName.replace('100.root','').replace('h_step2','out_hdbscan_allLCs')+"evt"+str(int(event))+"_v4.npz"
    else:
        outfilename="out_hdbscan_200PU_pt10_pt20_nonMipLCs_evt"+str(int(event))+"_"+whichSide+"_v3.npz"

    np.savez(outfilename,
             Labels=labels,
             Probabilities=probabilities,
             Graph=g,
             Dendrogram=d,
             CLMap=clMap)
