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

inputFileName="h_step2_eta2p5_eta2p5_pt10_pt20_p211_PU200_100.root"

inputFile=inputFileDir+inputFileName

fil=ROOT.TFile(inputFile)

lcTree=fil.Get("hgcalHitNtuple/LCTree")
lcInfo=lcTree.GetBranch("lcInfo")

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

for event in range(1, 101):

    print("Starting to Process Event:", event)
    
    lcs_in_event=LCs[LCs[:,0]==event]

    if isAllLC:
        lcmask = (lcs_in_event[:, 3]>0) & (2.2<lcs_in_event[:,5]) & (lcs_in_event[:,5]<2.8) & (-0.5<lcs_in_event[:,6]) & (lcs_in_event[:,6]<0.5)
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
        outfilename="out_hdbscan_200PU_pt10_pt20_evt"+str(int(event))+"_"+whichSide+"_v3.npz"
    else:
        outfilename="out_hdbscan_200PU_pt10_pt20_nonMipLCs_evt"+str(int(event))+"_"+whichSide+"_v3.npz"

    np.savez(outfilename,
             Labels=labels,
             Probabilities=probabilities,
             Graph=g,
             Dendrogram=d,
             CLMap=clMap)
