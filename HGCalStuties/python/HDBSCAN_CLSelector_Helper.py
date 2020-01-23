import math, sys, ROOT
import numpy as np

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
import matplotlib.axes as axs
import networkx as nx

def is_end_of_branch(cls, mother_cls):
    endOfBranch=[]
    for cl in cls:
        if not cl in mother_cls: 
            endOfBranch+=[True]
        else:
            endOfBranch+=[False]
    return np.array(endOfBranch)

def find_descendants(node, clusters, mothers_idx):
    cls = node
    while True: 
        for cl in cls:
            if cl in mothers_idx:
                cls=cls[~(cls == cl)]
                daughters_mask = clusters[:,0] == cl
                daughters = clusters[daughters_mask]
                node=np.concatenate((daughters[:,1], node), axis=0)
                cls = np.concatenate((daughters[:,1], cls), axis=0)
                if len(cls)-np.count_nonzero(is_end_of_branch(cls, mothers_idx)) == 0: break
        if len(cls)-np.count_nonzero(is_end_of_branch(cls, mothers_idx)) == 0: break
    node=node.flatten()
    return node

def find_children(node, clusters, mothers_idx):
    if node in mothers_idx:
        daughters_mask = clusters[:,0] == node
        daughters = clusters[daughters_mask]
    else: daughters = []
    return daughters

def calc_stability(cells, clusters, mothers_idx, daughters_idx):
    cl_stability_descendants=[]
    all_nodes = np.concatenate((mothers_idx, daughters_idx), axis=0)
    clusters_list = np.unique(all_nodes)
    for icl in clusters_list:
        if icl == len(cells): continue
        mother_mask = clusters[:,1] == icl
        mother = clusters[mother_mask]
        lambda_birth = mother[0][2]
        cluster_mask = clusters[:,0] == icl
        cluster = clusters[cluster_mask]
        if len(cluster)>0:
            lambda_death = cluster[0][2]
        #else: lambda_death = 9999
        else:
            member_mask = cells[:,0] == icl
            member = cells[member_mask]
            lambda_death = max(member[:,2])
        dlambda = lambda_death - lambda_birth
        cluster = cluster[:, 1]
        if False:    ### use for debugging
            if len(cluster)>0:
                print("-----", icl, cluster[0][3])
                print(find_children(icl, clusters, mothers_idx)[:,1])
            else:
                mask = cells[:,0] == icl
                c = cells[mask]
                print("-----", icl, len(c))
                print(find_children(icl, clusters, mothers_idx))
        
        descendants = find_descendants(cluster, clusters, mothers_idx)
        if len(descendants)>0:
            n=0
            for desc in descendants:
                desc_member_mask = cells[:,0] == desc
                desc_member = cells[desc_member_mask]
                n+=len(desc_member)
            stability = n * (dlambda)
            nhits_desc = n
        else:
            stability = 0
            nhits_desc = 0
        #print(dlambda)
        cluster_member_mask = cells[:,0] == icl
        cluster_member = cells[cluster_member_mask]
        nhits = nhits_desc+len(cluster_member)
        stability+=sum(cluster_member[:,2])-len(cluster_member[:,2])*lambda_birth
        desc_stability = 0
        cl_stability_descendants+=[[icl, stability, descendants, desc_stability, nhits, nhits_desc, lambda_birth, lambda_death]]
    
    cl_stability_descendants = np.array(cl_stability_descendants)

    for i in range(len(cl_stability_descendants)):
        descs = cl_stability_descendants[i][2]
        node = cl_stability_descendants[i][0]
        desc_stability = 0
        if len(find_children(node, clusters, mothers_idx))>0: children = find_children(node, clusters, mothers_idx)[:,1]
        else: children = []
        for d in descs:
        #for d in children:
            cl_mask = cl_stability_descendants[:,0] == d
            cl = cl_stability_descendants[cl_mask]
            if len(cl)>0: desc_stability+=cl[0][1]
        cl_stability_descendants[i][3] = desc_stability

    return cl_stability_descendants

def calc_stability_dlambda(cells, clusters, mothers_idx, daughters_idx):
    cl_stability_descendants=[]
    all_nodes = np.concatenate((mothers_idx, daughters_idx), axis=0)
    clusters_list = np.unique(all_nodes)
    for icl in clusters_list:
        if icl == len(cells): continue
        mother_mask = clusters[:,1] == icl
        mother = clusters[mother_mask]
        lambda_birth = mother[0][2]
        cluster_mask = clusters[:,0] == icl
        cluster = clusters[cluster_mask]
        if len(cluster)>0:
            lambda_death = cluster[0][2]
        else:
            member_mask = cells[:,0] == icl
            member = cells[member_mask]
            lambda_death = max(member[:,2])
        dlambda = lambda_death - lambda_birth
        stability = dlambda
        cluster = cluster[:, 1]
        descendants = find_descendants(cluster, clusters, mothers_idx)
        desc_stability = 0
        if len(descendants)>0:
            n=0
            for desc in descendants:
                desc_member_mask = cells[:,0] == desc
                desc_member = cells[desc_member_mask]
                n+=len(desc_member)
            integral = n * (dlambda)
        else: integral = 0
        cluster_member_mask = cells[:,0] == icl
        cluster_member = cells[cluster_member_mask]
        integral+=sum(cluster_member[:,2])-len(cluster_member[:,2])*lambda_birth
        cl_stability_descendants+=[[icl, stability, descendants, desc_stability, integral]]

    cl_stability_descendants = np.array(cl_stability_descendants)
    for i in range(len(cl_stability_descendants)):
        descs = cl_stability_descendants[i][2]
        node = cl_stability_descendants[i][0]
        desc_stability = 0
        d_stability = []
        for d in descs:
            cl_mask = cl_stability_descendants[:,0] == d
            cl = cl_stability_descendants[cl_mask]
            if len(cl)>0:
                desc_stability+=cl[0][1]
                d_stability = [cl[0][1]]
        if len(descs) > 0:
            cl_stability_descendants[i][3] = desc_stability/len(descs)
        else:
            cl_stability_descendants[i][3] = 0

    return cl_stability_descendants

def select_leaf_clusters(cl_stability_descendants, daughters_idx):
    leaf_clusters = []

    for cl in cl_stability_descendants:
        if len(cl[2])==0: leaf_clusters += [cl]

    return leaf_clusters
    
def select_eom_clusters(cl_stability_descendants, daughters_idx):
    eom_clusters = []
    clusters_ref = np.unique(daughters_idx)
    while len(clusters_ref) > 0:
        for icl in clusters_ref:
            mask = cl_stability_descendants[:,0] == icl
            cl = cl_stability_descendants[mask][0]
            clusters_ref=clusters_ref[~(clusters_ref==cl[0])]
            print("Debug ***",cl[0], cl[1], cl[2], cl[3], cl[4], cl[5], cl[6], cl[7])
            #if (cl[1] > cl[3] and cl[3] < 10000 and cl[6] > 20) or (cl[1] > cl[3] and len(cl[2]) == 0):
            # v0 final result for 200 PU
            if (cl[1] > cl[3] and cl[3] < 10000 and cl[6] > 20 and len(cl[2]) < 5) or (cl[1] > cl[3] and len(cl[2]) == 0):
            # v1 (v2) with dlambda cuts at 10 (3)
            #if (cl[1] > cl[3] and cl[3] < 10000 and cl[6] > 20 and len(cl[2]) < 4) or (cl[1] > cl[3] and len(cl[2]) == 0) or (cl[7] - cl[6] < 10 and cl[3] < 10000 and cl[6] > 20 and len(cl[2]) < 4):
            #### used for 0PU final result (v4)
            #if (cl[1] > cl[3] and cl[3] < 10000 and cl[6] > 20) or (cl[1] > cl[3] and len(cl[2]) == 0) or (cl[7] - cl[6] < 10 and cl[3] < 10000 and cl[6] > 20):
                print("Debug ---",cl[0], cl[1], cl[2], cl[3], cl[4], cl[5], cl[6], cl[7])
                eom_clusters+=[cl]
                for desc in cl[2]:
                    clusters_ref=clusters_ref[~(clusters_ref==desc)]
                break
            
    return eom_clusters

def assign_cluster_labels(selected_clusters, cells):
    i=0
    for cl in selected_clusters:
        cl_idx=cl[0]
        desc_idx = cl[2]
        for c in cells:
            if c[0] == cl_idx or c[0] in desc_idx:
                c[4] = i
        i+=1
        #print(i)
    cells=cells[cells[:, 1].argsort()]
    return cells

def select_new_clusters(cells, clusters, mothers_idx, families_matched_node, mother_matched_node):
    new_clusters = []
    for cl in clusters:
        mother_idx = cl[0]
        daughter_idx = cl[1]
        if (not daughter_idx in mothers_idx and not daughter_idx in families_matched_node) or daughter_idx == mother_matched_node:
            new_clusters+=[cl]
            
    new_clusters=np.array(new_clusters)
    
    i=0
    for cl in new_clusters:
        cl_idx=cl[1]
        if cl_idx == mother_matched_node:
            for c in cells:
                if c[0] == cl_idx or c[0] in families_matched_node: c[4] = i
        else: 
            for c in cells:
                if c[0] == cl_idx: c[4] = i
        i+=1
        
    new_sorted_cells=cells[cells[:, 1].argsort()]
    return new_sorted_cells
