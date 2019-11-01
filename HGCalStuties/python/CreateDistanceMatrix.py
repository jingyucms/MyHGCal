import numpy as np
import ROOT, math, sys, os

def distanceMatrices(x, y, isLC):
    dl = abs(x[8]-y[8])
    deta=x[5]-y[5]
    dphi=x[6]-y[6]
    if dphi > math.pi:
        dphi -= 2*math.pi
    if dphi < -math.pi:
        dphi += 2*math.pi
    dist_angle=math.sqrt(deta**2+dphi**2)

    if isLC:
        w=math.sqrt(x[7])
    else:
        w=1+x[4]

    if dl>=10:
        distance = 9999
    elif x[8] > 28 and x[8] <= 40:
        distance =  dist_angle/(2*w)
    elif x[8] > 40:
        distance =  dist_angle/(4*w)
    else:
        distance = dist_angle/(w)

    return distance

isLCs=False
whichSide='positive'

inputFileDir='../test/'
inputFileName="h_step2_eta2p5_eta2p5_pt10_pt20_p211_PU200_100.root"
inputFile=inputFileDir+inputFileName

fil=ROOT.TFile(inputFile)

if isLCs:
    lcTree=fil.Get("hgcalHitNtuple/LCTree")
    lcInfo=lcTree.GetBranch("lcInfo")

    whichSide='positive'   ### For now, only use positive
    #whichSide='negative'

    LCs = np.asarray([[lcInfo.GetLeaf("event").GetValue(),
                       lcInfo.GetLeaf("x").GetValue(),
                       lcInfo.GetLeaf("y").GetValue(),
                       lcInfo.GetLeaf("z").GetValue(),
                       lcInfo.GetLeaf("energy").GetValue(),
                       lcInfo.GetLeaf("eta").GetValue(),
                       lcInfo.GetLeaf("phi").GetValue(),
                       lcInfo.GetLeaf("size").GetValue(),
                       lcInfo.GetLeaf("layer").GetValue()] for lc in lcTree])
else:
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

#print("Size of the layer cluster collection:",len(L))

#for event in range(19, 26):
#for event in range(41, 51):
#for event in range(63, 76):
#for event in range(87, 101):
for event in [75, 100]:

    print("Starting to Process Event:", event)

    if isLCs:
        lcs_in_event=LCs[LCs[:,0]==event]
        mask = (lcs_in_event[:, 3]>0) & (2.2<lcs_in_event[:,5]) & (lcs_in_event[:,5]<2.8) & (-0.5<lcs_in_event[:,6]) & (lcs_in_event[:,6]<0.5)
        L=lcs_in_event[mask]
        X=L
    else:
        rechits_in_event=RHs[RHs[:,0]==event]
        mask = (rechits_in_event[:, 3]>0) & (2.2<rechits_in_event[:,5]) & (rechits_in_event[:,5]<2.8) & (-0.3<rechits_in_event[:,6]) & (rechits_in_event[:,6]<0.3)
        R=rechits_in_event[mask]
        X=R

    n=len(X)
    print("Number of Inputs:", n)
    matrix=np.zeros((n,n), dtype=np.float32)

    #sys.exit()

    #continue
    
    for i in range(n):
        for j in range(n):
            x=X[i]
            y=X[j]
            matrix[i][j]=distanceMatrices(x, y, isLCs)

    saveDir='root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/HGC/HDBSCAN/'
    
    if isLCs:
        saveName='distance_matrix_hdbscan_allLCs_'+whichSide+'_evt'+str(int(event))+'_v3.npy'
    else:
        saveName='distance_matrix_hdbscan_allRecHits_'+whichSide+'_evt'+str(int(event))+'_v3.npy'

    print("Output Transfer To:", saveDir+saveName)
    np.save(saveName, matrix)

    os.system("xrdcp "+saveName+" "+saveDir+saveName)

    os.system("rm "+saveName)
