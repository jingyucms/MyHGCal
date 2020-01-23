import math, sys, ROOT
import numpy as np

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
import matplotlib.axes as axs
import networkx as nx
from mpl_toolkits.mplot3d import Axes3D
import itertools

from HDBSCAN_CLSelector_Helper import *
from HDBSCAN_PerfStudy_Helper import *

ROOT.gROOT.SetBatch()

def correction(energy, offset, slope):
    corrected_energy = offset / slope + energy / slope
    return corrected_energy

energy=300
radius=50

if radius == 100:
    offset = 2.4
    slope = 0.81
elif radius == 50:
    offset = 4.4
    slope = 0.72

events = np.linspace(1,100,100)
    
inputFileDir='../test/'
    
inputFileName="h_step2_r"+str(radius)+"_z320_e"+str(energy)+"_p211_PU200_100.root"
    
inputFile=inputFileDir+inputFileName
    
fil=ROOT.TFile(inputFile)
    
isOwnSelection=False
isRecHits=True
isAllLC=True
whichSide='positive' ### Don't change this

if __name__ == "__main__":
    
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

    h_response = ROOT.TH1F("Response", "Response", 20, 0, 2)
    h_response_rec = ROOT.TH1F("Response_rec", "Response_rec", 20, 0, 2)
    h_response_before = ROOT.TH1F("Response_before", "Response_before", 20, 0, 2)

    h_energy_reco = ROOT.TH1F("Energy", "Energy", 2000, 0, 1000)
    h_energy_corr = ROOT.TH1F("Energy_Corr", "Energy_Corr", 2000, 0, 1000)
    h_energy_rec = ROOT.TH1F("Energy_rec", "Energy_rec", 2000, 0, 1000)
    
    for event in events:
        print("\n")
        print("******** Start to process event:",int(event),"********")
        #if not event in [10, 27, 43, 44]: continue
        #if not event in [68, 98]: continue
        #if not event in [61]: continue

        cps_in_event=CPs[CPs[:,0]==event]
        signal_cp1 = cps_in_event[0]
        signal_cp2 = cps_in_event[1]
        
        if signal_cp1[4]>0:
            signal_cp = signal_cp1
        else:
            signal_cp = signal_cp2

        cp_eta = signal_cp[4]
        cp_phi = signal_cp[5]
        
        if isRecHits:
            rechits_in_event=RHs[RHs[:,0]==event]
            #mask = (rechits_in_event[:, 3]>0) & (2.2<rechits_in_event[:,5]) & (rechits_in_event[:,5]<2.8) & (-0.3<rechits_in_event[:,6]) & (rechits_in_event[:,6]<0.3)
            #RH=rechits_in_event[mask]
            RH=rechits_in_event
            outfilename="./data/out_hdbscan_allRecHits_r"+str(radius)+"_z320_e"+str(energy)+"_p211_PU200_100_evt"+str(int(event))+"_v4.npz"
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
            mothers_idx = clusters[:,0]
            daughters_idx = clusters[:,1]
            #cl_stability_descendants = calc_stability_dlambda(cells, clusters, mothers_idx, daughters_idx)
            cl_stability_descendants = calc_stability(cells, clusters, mothers_idx, daughters_idx)
            leaf_clusters = select_leaf_clusters(cl_stability_descendants, daughters_idx)
            eom_clusters = select_eom_clusters(cl_stability_descendants, daughters_idx)
            leaf_labels = assign_cluster_labels(leaf_clusters, cells)
            eom_labels = assign_cluster_labels(eom_clusters, cells)
            labels_leaf = leaf_labels[:, 4]
            labels_eom = eom_labels[:, 4]
        else:
            lcs_in_event=LCs[LCs[:,0]==event]
            if isAllLC:
                #lcmask = (lcs_in_event[:, 3]>0) & (2.2<lcs_in_event[:,5]) & (lcs_in_event[:,5]<2.8) & (-0.5<lcs_in_event[:,6]) & (lcs_in_event[:,6]<0.5)
                #outfilename="out_hdbscan_200PU_pt10_pt20_evt"+str(int(event))+"_"+whichSide+"_v3.npz"
                lcmask = (lcs_in_event[:, 3]>0) & (cp_eta-0.3<lcs_in_event[:,5]) & (lcs_in_event[:,5]<cp_eta+0.3) & (cp_phi-0.3<lcs_in_event[:,6]) & (lcs_in_event[:,6]<cp_phi+0.3)
                outfilename="./data/out_hdbscan_allLCs_r"+str(radius)+"_z320_e"+str(energy)+"_p211_PU200_evt"+str(int(event))+"_v4.npz"
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
            print(len(sample), len(labels_leaf), len(labels_eom))
            print("------ERROR------!!!")


        labels = labels_eom
        unique_labels = set(labels)
            
        if True:
            clusters3d = plt.figure(0, figsize=(40,30))
            ax_clusters = Axes3D(clusters3d)
            colors = itertools.cycle(['r', 'b', 'g','m','c','y','k'])
        
            for k in unique_labels:
                color=next(colors)
                size=100
                if k == -1:
                    # Black used for noise.
                    color = [0, 0, 0, 1]
                    size = 10
                if not len(labels)==0:
                    class_member_mask = (labels == k)
                    xyz = sample[class_member_mask]
                    ax_clusters.scatter(xyz[:, 2], xyz[:, 3], xyz[:, 1], 'o', color=color, s=size)
        
            ax_clusters.set_ylim(300, 500)
            ax_clusters.set_zlim(-150, 150)
            ax_clusters.set_xlim(-150, 150)
            ax_clusters.set_ylabel("z", size=30)
            ax_clusters.set_xlabel("y", size=30)
            ax_clusters.set_zlabel("x", size=30)
            ax_clusters.view_init(30,30)
            plt.savefig(outfilename.replace(".npz", "_3d_evt"+str(int(event))+".pdf").replace('data', 'plots'))
            clusters3d.clear()
            ax_clusters.clear()
    
        if True:
            clusters2d = plt.figure(0, figsize=(10,10))
            colors = itertools.cycle(['r', 'b', 'g','m','c','y','k'])
    
            for k in unique_labels:
                color=next(colors)
                size=10
                if k == -1:
                    # Black used for noise.
                    color = [0, 0, 0, 1]
                    size = 1
                if not len(labels)==0:
                    class_member_mask = (labels == k)
                    xyz = sample[class_member_mask]
                    plt.scatter(xyz[:, 5], xyz[:, 6], color=color, s=size)
     
            plt.xlim(cp_eta-0.3, cp_eta+0.3)
            plt.ylim(cp_phi-0.3, cp_phi+0.3)
            plt.savefig(outfilename.replace(".npz", "_2d_evt"+str(int(event))+".pdf").replace('data', 'plots'))
            clusters2d.clear()

        sc_in_event=SCs[SCs[:,0]==event]
        #sc_energy_cut=20
        #cl_energy_cut=2
        distance_cut=0.1

        for icp in [0]:
            cpmask = CPs[:,0]==event
            CP = CPs[cpmask]
            cp=CP[CP[:,1]==icp][0]
            cp_energy=cp[2]
            cp_energy_rec=cp[3]
            h_energy_rec.Fill(cp_energy_rec)
            cp_eta=cp[4]
            cp_phi=cp[5]
            if cp_eta<0: continue
            sc_in_cp=sc_in_event[sc_in_event[:,7]==icp]
            scList=set(sc_in_cp[:,8])
            if len(scList) > 1: continue
##             for isc in scList:
##                 sc_mask = (sc_in_cp[:,8]==isc)
##                 sc=sc_in_cp[sc_mask]
##                 sc_energy=sum(sc[:,4])
##                 if sc_energy <= sc_energy_cut: continue
            #matched = do_matching(sample, labels_leaf, sc, sc_energy_cut, cl_energy_cut, distance_cut)
            #matched_eom = do_matching(sample, labels_eom, sc, sc_energy_cut, cl_energy_cut, distance_cut)
            matched = do_matching_cp(sample, labels_leaf, cp_eta, cp_phi, cp_energy_rec, distance_cut)
            matched_eom = do_matching_cp(sample, labels_eom, cp_eta, cp_phi, cp_energy_rec, distance_cut)
            print("matched_eom:", matched_eom, len(matched_eom))
            print("matched:", matched, len(matched))

            if len(matched_eom) == 1:
                response = 0
                responses_eom += [0]
                h_response.Fill(0)
                h_response_rec.Fill(0)
                h_energy_reco.Fill(0)
                h_energy_corr.Fill(correction(0, offset, slope))
            else:
                response_rec=matched_eom[2]/cp_energy_rec
                response=matched_eom[2]/cp_energy

                h_energy_reco.Fill(matched_eom[2])
                h_energy_corr.Fill(correction(matched_eom[2], offset, slope))
                
                if response > 2:
                    response=1.9999
                if response_rec > 2:
                    response_rec=1.9999
                responses_eom+=[response_rec]
                h_response_rec.Fill(response_rec)
                h_response.Fill(response)

                if response_rec > 1.6 or response_rec < 0.2:
                    print("DEBUG!!!")
                
            print("response_eom", response_rec)
            
            if len(matched)==1:
                response=0
                responses_before+=[response]
                responses_after+=[response]
                print("response_before", response)
                h_response_before.Fill(0)
                continue
            else:
                matched_idx = matched[1]
                response=matched[2]/cp_energy_rec
                if response > 2:
                    response=1.9999
                responses_before+=[response]
                
                print("response_before", response)
                h_response_before.Fill(response)
                leaf_clusters = np.array(leaf_clusters)
                matched_mask = leaf_labels[:, 4] == matched_idx
                matched_cluster = leaf_labels[matched_mask]

                matched_node = np.unique(matched_cluster[:,0])[0]
                lambda_value = clusters[clusters[:, 1] == matched_node][0][2]
                print("matched_node lambda_value:",matched_node, lambda_value)
                
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
                    
                if isOwnSelection:
                    while lambda_value >= lambda_value_cut:
                        matched_node_mask = clusters[:, 1] == matched_node
                        mother_matched_node = clusters[matched_node_mask][0][0]
                        mother_matched_node_mask = clusters[:, 1] == mother_matched_node
                        lambda_value = clusters[mother_matched_node_mask][0][2]
                        matched_node = mother_matched_node

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
        saveName="response_r"+str(radius)+"_z320_e"+str(energy)+"_p211_PU200_rechits_v4_v5.pdf"
    else:
        if isAllLC:
            saveName="response_r"+str(radius)+"_z320_e"+str(energy)+"_p211_PU200_lcs_v4.pdf"
        else:
            saveName="response_nonMipLCs_200PU_v3.pdf"
    #plt.savefig(saveName)
    presponse.clear()

    #print(responses_before, lambdas_diff)
    #plambda_diff = plt.figure(2)
    #plt.scatter(responses, lambdas_diff, color='b', s=30)
    #plt.savefig("lambdas_diff.pdf")
    #plambda_diff.clear()

    outfile=ROOT.TFile(saveName.replace('.pdf','.root'), 'recreate')

    outfile.cd()
    h_energy_rec.Write()
    h_energy_reco.Write()
    h_energy_corr.Write()
    h_response.Write()
    h_response_rec.Write()
    h_response_before.Write()

    outfile.Close()

    fitter = ROOT.TF1('fitter','gaus', 0, 1.2)
    
    ROOT.gROOT.ForceStyle()
    ROOT.gStyle.SetOptStat(1)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptFit(1111)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadRightMargin(0.04)
    ROOT.gStyle.SetPadLeftMargin(0.09)

    c=ROOT.TCanvas("myCanvas", "myCanvas")
    h_response.Fit('fitter', 'R')
    
    ROOT.gPad.Update()
    st=h_response.FindObject("stats")
    st.SetX1NDC(0.65)
    st.SetX2NDC(0.96)
    st.SetY1NDC(0.95)
    st.SetY2NDC(0.5)
    #st.SetTitle("Response")
    st.SetOptStat(1)
    st.SetOptFit(1111)
    h_response.Draw()
    h_response.SetMaximum(h_response_rec.GetMaximum()*1.25)
    h_response.SetLineColor(ROOT.kGreen+3)
    fitter.Draw("same")
    fitter.SetLineColor(ROOT.kGreen+3)
    #fitter.SetMaximum(60)
    h_response_rec.Draw("same")
    h_response_rec.SetLineColor(ROOT.kOrange+10)
    h_response_rec.SetLineWidth(3)
    h_response_rec.SetLineStyle(7)
    h_response_before.SetLineColor(1)
    h_response_before.Draw("same")
    h_response_before.SetLineWidth(3)
    h_response_before.SetLineStyle(7)

    lg = ROOT.TLegend(0.5, 0.6, 0.65, 0.95, str(energy)+" GeV")
    lg.SetTextSize(0.04)
    lg.SetFillStyle(0)
    lg.AddEntry(h_response, "Gen Energy")
    lg.AddEntry(h_response, "Reconstructable Energy")
    
    c.SaveAs(saveName)
