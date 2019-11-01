import math, sys, os
import numpy as np

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

import hdbscan

whichSide='positive'
matrixFileDir='root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/HGC/HDBSCAN/'

for event in range(1, 101):
    print("Starting to Process Event:", event)

    matrixFile='distance_matrix_hdbscan_allRecHits_'+whichSide+'_evt'+str(int(event))+'_v3.npy'

    os.system("xrdcp "+matrixFileDir+matrixFile+" "+matrixFile)

    print("Use Matrix:", matrixFile)

    matrix=np.load(matrixFile)

    db1 = hdbscan.HDBSCAN(min_samples=10, min_cluster_size=10, metric='precomputed').fit(matrix.astype(np.float64))
    labels1 = db1.labels_
    probabilities1 = db1.probabilities_

    #labels1_alt=db1.single_linkage_tree_.get_clusters(min_cluster_size=20)

    db2 = hdbscan.HDBSCAN(min_samples=20, min_cluster_size=20, metric='precomputed').fit(matrix.astype(np.float64))
    labels2 = db2.labels_
    probabilities2 = db2.probabilities_

    db3 = hdbscan.HDBSCAN(cluster_selection_method='leaf', min_samples=10, min_cluster_size=10, metric='precomputed').fit(matrix.astype(np.float64))
    labels3 = db3.labels_
    probabilities3 = db3.probabilities_

    #labels3_alt=db3.single_linkage_tree_.get_clusters(min_cluster_size=20)

    db4 = hdbscan.HDBSCAN(cluster_selection_method='leaf', min_samples=20, min_cluster_size=20, metric='precomputed').fit(matrix.astype(np.float64))
    labels4 = db4.labels_
    probabilities4 = db4.probabilities_

    db5 = hdbscan.HDBSCAN(min_samples=10, min_cluster_size=20, metric='precomputed').fit(matrix.astype(np.float64))
    labels5 = db5.labels_
    probabilities5 = db5.probabilities_
    
    db6 = hdbscan.HDBSCAN(cluster_selection_method='leaf', min_samples=10, min_cluster_size=20, metric='precomputed').fit(matrix.astype(np.float64))
    labels6 = db6.labels_
    probabilities6 = db6.probabilities_

    saveName="out_hdbscan_200PU_pt10_pt20_allRecHits_"+whichSide+"_evt"+str(int(event))+"_v3.npz"
    #np.savez(saveName, Labels=labels, Probabilities=probabilities, Labels_Alt1=labels_alt1, Labels_Alt2=labels_alt2)
    np.savez(saveName,
             Labels1=labels1,
             Probabilities1=probabilities1,
             #Labels1_Alt=labels1_alt,
             Labels2=labels2,
             Probabilities2=probabilities2,
             Labels3=labels3,
             Probabilities3=probabilities3,
             #Labels3_Alt=labels1_alt,
             Labels4=labels4,
             Probabilities4=probabilities4,
             Labels5=labels5,
             Probabilities5=probabilities5,
             Labels6=labels6,
             Probabilities6=probabilities6)

    os.system("rm "+matrixFile)
