import math, sys, ROOT
import numpy as np

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
import matplotlib.axes as axs
import networkx as nx

from HDBSCAN_CLSelector_Helper import *

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

if __name__ == "__main__":
    events = np.linspace(1,100,100)
    
    inputFileDir='../test/'
    
    inputFileName="h_step2_eta2p5_eta2p5_pt10_pt20_p211_PU200_100.root"
    
    inputFile=inputFileDir+inputFileName
    
    fil=ROOT.TFile(inputFile)
    
    isOwnSelection=False
    isRecHits=True
    isAllLC=True
    whichSide='positive' ### Don't change this
    
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

    responses_eom=[]
    responses_before=[]
    responses_after=[]
    lambdas_diff=[]
    responses=[]
    for event in events:
        print("\n")
        print("******** Start to process event:",int(event),"********")
        #if not event in [10, 27, 43, 44]: continue
        #if not event in [68, 98]: continue
    
        if isRecHits:
            rechits_in_event=RHs[RHs[:,0]==event]
            mask = (rechits_in_event[:, 3]>0) & (2.2<rechits_in_event[:,5]) & (rechits_in_event[:,5]<2.8) & (-0.3<rechits_in_event[:,6]) & (rechits_in_event[:,6]<0.3)
            RH=rechits_in_event[mask]
            outfilename="out_hdbscan_200PU_pt10_pt20_allRecHits_"+whichSide+"_evt"+str(int(event))+"_v3.npz"
            print(outfilename)
            data = np.load(outfilename, allow_pickle=True)
            sample=RH
            g=data['Graph']
            d=data['Dendrogram']
            D=np.zeros((len(d), 5), dtype=np.float32)
            D[:,:-1]=d
            D[:,-1]=-1
            cells=D[D[:,3]==1]
            clusters=D[D[:,3]!=1]
            #print(clusters)
            mothers_idx = clusters[:,0]
            daughters_idx = clusters[:,1]
            #print(mothers_idx)
            #print(daughters_idx)
            cl_stability_descendants = calc_stability_dlambda(cells, clusters, mothers_idx, daughters_idx)
            leaf_clusters = select_leaf_clusters(cl_stability_descendants, daughters_idx)
            eom_clusters = select_eom_clusters(cl_stability_descendants, daughters_idx)
            #eom_clusters = select_eom_clusters(cl_stability_descendants, cells[:,0])
            leaf_labels = assign_cluster_labels(leaf_clusters, cells)
            eom_labels = assign_cluster_labels(eom_clusters, cells)
            labels_leaf = leaf_labels[:, 4]
            labels_eom = eom_labels[:, 4]
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
            sample=LC
            g=data['Graph']
            d=data['Dendrogram']
            D=np.zeros((len(d), 5), dtype=np.float32)
            D[:,:-1]=d
            D[:,-1]=-1
            cells=D[D[:,3]==1]
            clusters=D[D[:,3]!=1]
            mothers_idx = clusters[:,0]
            daughters_idx = clusters[:,1]
            cl_stability_descendants = calc_stability(cells, clusters, mothers_idx, daughters_idx)
            leaf_clusters = select_leaf_clusters(cl_stability_descendants, daughters_idx)
            eom_clusters = select_eom_clusters(cl_stability_descendants, daughters_idx)
            leaf_labels = assign_cluster_labels(leaf_clusters, cells)
            eom_labels = assign_cluster_labels(eom_clusters, cells)
            labels_leaf = leaf_labels[:, 4]
            labels_eom = eom_labels[:, 4]
    
        if not len(sample)==len(labels_leaf) or not len(sample)==len(labels_eom):
            print("------ERROR------!!!")

        sc_in_event=SCs[SCs[:,0]==event]
        sc_energy_cut=20
        cl_energy_cut=2
        distance_cut=0.1

        for icp in [0,1]:
            cpmask = CPs[:,0]==event
            CP = CPs[cpmask]
            cp=CP[CP[:,1]==icp][0]
            cp_energy=cp[2]
            cp_rec_energy=cp[3]
            cp_eta=cp[4]
            cp_phi=cp[5]
            if cp_eta<0: continue
            sc_in_cp=sc_in_event[sc_in_event[:,7]==icp]
            scList=set(sc_in_cp[:,8])
            if len(scList) > 1: continue
            for isc in scList:
                sc_mask = (sc_in_cp[:,8]==isc)
                sc=sc_in_cp[sc_mask]
                sc_energy=sum(sc[:,4])
                if sc_energy <= sc_energy_cut: continue
                matched = do_matching(sample, labels_leaf, sc, sc_energy_cut, cl_energy_cut, distance_cut)
                matched_eom = do_matching(sample, labels_eom, sc, sc_energy_cut, cl_energy_cut, distance_cut)
                print("matched_eom:", matched_eom, len(matched_eom))
                print("matched:", matched, len(matched))

                if len(matched_eom) == 1:
                    response = 0
                    responses_eom += [response]
                else:
                    response=matched_eom[2]/matched_eom[3]
                    if response > 2:
                        response=2
                    responses_eom+=[response]
                print("response_eom", response)
                
                if len(matched)==1:
                    response=0
                    responses_before+=[response]
                    responses_after+=[response]
                    print("response_before", response)
                    continue
                else:
                    matched_idx = matched[1]
                    response=matched[2]/matched[3]
                    if response > 2:
                        response=2
                    responses_before+=[response]
                    print("response_before", response)
                    
                    leaf_clusters = np.array(leaf_clusters)
                    matched_mask = leaf_labels[:, 4] == matched_idx
                    matched_cluster = leaf_labels[matched_mask]
    
                    matched_node = np.unique(matched_cluster[:,0])[0]
                    lambda_value = clusters[clusters[:, 1] == matched_node][0][2]
                    print("matched_node lambda_value:",matched_node, lambda_value)
                    
                    #lambda_value_cut = 120
                    #if response<0.4: lambda_value_cut=120
                    #else: lambda_value_cut = lambda_value+1
                    
                    #if lambda_value < lambda_value_cut:
                    #    responses_after+=[response]
                    #    continue

                    if response > 0.5:
                        responses_after+=[response]
                        continue

                    matched_node_mask = clusters[:, 1] == matched_node
                    mother_matched_node = clusters[matched_node_mask][0][0]
                    mother_matched_node_mask = clusters[:, 1] == mother_matched_node
                    lambda_value_new = clusters[mother_matched_node_mask][0][2]
                    matched_node = mother_matched_node

                    print("lambda_value_new:",lambda_value_new)

                    lambdas_diff+=[lambda_value-lambda_value_new]
                    responses+=[response]

                    #if lambda_value-lambda_value_new < 20:
                    #    responses_after+=[response]
                    #    continue
                        
                    if isOwnSelection:
#                        while lambda_value >= lambda_value_cut:
#                            matched_node_mask = clusters[:, 1] == matched_node
#                            mother_matched_node = clusters[matched_node_mask][0][0]
#                            mother_matched_node_mask = clusters[:, 1] == mother_matched_node
#                            lambda_value = clusters[mother_matched_node_mask][0][2]
#                            matched_node = mother_matched_node
    
                        families_matched_node_mask = clusters[:, 0] == mother_matched_node
                        families_matched_node = clusters[families_matched_node_mask][:,1]
                        families_matched_node=find_descendants(families_matched_node, clusters, mothers_idx)
                        print("families_matched_node:", families_matched_node)
    
                        new_clusters = select_new_clusters(cells, clusters, mothers_idx, families_matched_node, mother_matched_node)
                        new_labels = new_clusters[:,4]
                        new_matched = do_matching(sample, new_labels, sc, sc_energy_cut, cl_energy_cut, distance_cut)
    
                        if len(new_matched)==1:
                            response=0
                            responses_after+=[response]
                            print("response_after", response)
                        else:
                            response=new_matched[2]/new_matched[3]      
                            if response > 2:
                                response=2
                            responses_after+=[response]
                            print("response_after",response)
        
        
    bins=np.linspace(0, 20, 21)*0.1

    presponse=plt.figure(1)
    plt.hist(responses_before, bins, histtype='step')
    plt.hist(responses_eom, bins, histtype='step')
    if isOwnSelection:
        plt.hist(responses_after, bins, histtype='step', fill=False)
    plt.xlim(0,2)
    if isRecHits:
        saveName="response_recHits_200PU_v3.pdf"
    else:
        if isAllLC:
            saveName="response_lcs_200PU_v3.pdf"
        else:
            saveName="response_nonMipLCs_200PU_v3.pdf"
    plt.savefig(saveName)
    presponse.clear()

    #print(responses_before, lambdas_diff)
    plambda_diff = plt.figure(2)
    plt.scatter(responses, lambdas_diff, color='b', s=30)
    plt.savefig("lambdas_diff.pdf")
    plambda_diff.clear()
