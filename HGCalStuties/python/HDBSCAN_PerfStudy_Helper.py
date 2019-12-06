import math, sys, ROOT
import numpy as np

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
import matplotlib.axes as axs
import networkx as nx



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

def do_matching(sample, labels, sc, sc_energy_cut, cl_energy_cut, distance_cut):
    cluster_labels=set(labels[labels[:]>=0])
    sc_energy=sum(sc[:,4])
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
        distances+=[[distance(sc_eta, sc_phi, cl_eta, cl_phi), cl, cl_energy, sc_energy]]

    distances=np.array(distances)
    #print(distances)
    if not len(distances)==0 and not min(distances[:,0]) > distance_cut:
        m = distances[distances[:,0]==min(distances[:,0])]
        matched=m[0]
        cluster_labels.remove(matched[1])
    else:
        matched=[0]
    return matched
