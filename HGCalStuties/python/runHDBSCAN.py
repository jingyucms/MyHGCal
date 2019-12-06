import math, sys, os
import numpy as np
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import hdbscan
import networkx as nx

whichSide='positive'
matrixFileDir='root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/HGC/HDBSCAN/'

for event in range(1, 101):
    print("Starting to Process Event:", event)
    #if not event==1: continue

    matrixFile='distance_matrix_hdbscan_allRecHits_'+whichSide+'_evt'+str(int(event))+'_v3.npy'

    os.system("xrdcp "+matrixFileDir+matrixFile+" "+matrixFile)

    print("Use Matrix:", matrixFile)

    matrix=np.load(matrixFile)

    clusterer = hdbscan.HDBSCAN(min_samples=10, min_cluster_size=20, metric='precomputed').fit(matrix.astype(np.float64))
    labels = clusterer.labels_
    probabilities = clusterer.probabilities_

    labels_alt=clusterer.single_linkage_tree_.get_clusters(1/100, min_cluster_size=20)

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

    saveName="out_hdbscan_200PU_pt10_pt20_allRecHits_"+whichSide+"_evt"+str(int(event))+"_v3.npz"


    np.savez(saveName,
             Labels=labels,
             Probabilities=probabilities,
             Labels_Alt=labels_alt,
             Graph=g,
             Dendrogram=d,
             CLMap=clMap)

    os.system("rm "+matrixFile)
