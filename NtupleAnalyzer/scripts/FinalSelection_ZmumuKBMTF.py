from ROOT import RDataFrame, TFile, TChain, TTree, TFile, TH1D, TLorentzVector,TCanvas
import numpy as np
import sys
from math import cos,sin,sqrt,pi
import ROOT
import time as timer
time_start=timer.time()
ROOT.gInterpreter.AddIncludePath('/eos/user/f/flarover/HSCP_2025/HSCPanalysis/CMSSW_15_0_10/src/L1ScoutingAnalysisRDataFrame/NtupleAnalyzer/lib')

sim = True

if (sim):
    isdata = False
    ROOT.gInterpreter.Declare('#include "MODbasic_sel.h"')
    ROOT.gInterpreter.ProcessLine('.L /eos/user/f/flarover/HSCP_2025/HSCPanalysis/CMSSW_15_0_10/src/L1ScoutingAnalysisRDataFrame/NtupleAnalyzer/lib/MODbasic_sel.cpp+')

else:
    isdata = True
    ROOT.gInterpreter.Declare('#include "MODbasic_sel.h"')
    ROOT.gInterpreter.ProcessLine('.L /eos/user/f/flarover/HSCP_2025/HSCPanalysis/CMSSW_15_0_10/src/L1ScoutingAnalysisRDataFrame/NtupleAnalyzer/lib/MODbasic_sel.cpp+')
#ROOT.gSystem.Load('/eos/user/f/flarover/HSCP_2025/HSCPanalysis/CMSSW_15_0_10/src/L1ScoutingAnalysisRDataFrame/NtupleAnalyzer/lib/basic_sel_cpp.so')


### sample: output root file name
input_file = sys.argv[1]
output_file = sys.argv[2]

df = RDataFrame(0)

weight=1.0


#df = RDataFrame("Events",input_file)

#print(df.Count().GetValue())

# if ("Scouting" not in input_file):
#     isdata=False
#     df = RDataFrame("Events",input_file)
#     weight=(6346.0/df.Count().GetValue())

#     if "DYMM" in input_file: weight=((6346.0*0.3664)/df.Count().GetValue())
# else:
df = RDataFrame("Events",input_file)

nentries = df.Count().GetValue()

print ("isdata ", isdata)

if (isdata):
    weight = 1.0

print ("Before selection total entries", nentries)

df = df.Filter("nL1KBMTFSkimmed>1")

df_var=df.Define("idx1","GetIndex_nostub_hwK(1,nL1KBMTFSkimmed, L1KBMTFSkimmed_pt, L1KBMTFSkimmed_eta, L1KBMTFSkimmed_phi)").Define("idx2","GetIndex_nostub_hwK(2,nL1KBMTFSkimmed, L1KBMTFSkimmed_pt, L1KBMTFSkimmed_eta, L1KBMTFSkimmed_phi)")

df_var=df_var.Define("my_mu1","GetLepVector_hwK(idx1,L1KBMTFSkimmed_eta,L1KBMTFSkimmed_phi,L1KBMTFSkimmed_pt)").Define("my_mu2","GetLepVector_hwK(idx2,L1KBMTFSkimmed_eta,L1KBMTFSkimmed_phi,L1KBMTFSkimmed_pt)").Define("isOS","L1KBMTFSkimmed_hwCharge[idx1]*L1KBMTFSkimmed_hwCharge[idx2]<0")

df_var = df_var.Filter("idx1!=idx2 && my_mu1.Pt()>20 && my_mu2.Pt()>20 && fabs(my_mu1.Eta())<0.83 && fabs(my_mu2.Eta())<0.83")

df_var = df_var.Define("mmumu","(my_mu1+my_mu2).M()").Define("DRmumu","my_mu1.DeltaR(my_mu2)")

df = df_var.Filter("mmumu>50 && DRmumu>0.3").Define("xsweight","{}".format(weight))

df = df.Define("bxspread1", "GetBxSpread(nL1KBMTFSkimmed, idx1, L1KBMTFSkimmed_nStub, L1KBMTFSkimmed_s1Bx, L1KBMTFSkimmed_s2Bx, L1KBMTFSkimmed_s3Bx, L1KBMTFSkimmed_s4Bx)")\
        .Define("bxspread2", "GetBxSpread(nL1KBMTFSkimmed, idx2, L1KBMTFSkimmed_nStub, L1KBMTFSkimmed_s1Bx, L1KBMTFSkimmed_s2Bx, L1KBMTFSkimmed_s3Bx, L1KBMTFSkimmed_s4Bx)")\
        .Define("isL1MuMatched1", "IsL1MuMatched(nL1KBMTFSkimmed, idx1,L1KBMTFSkimmed_pt,L1KBMTFSkimmed_eta,L1KBMTFSkimmed_phi,nSkimmedL1Mu,SkimmedL1Mu_pt,SkimmedL1Mu_eta,SkimmedL1Mu_phi)")\
        .Define("isL1MuMatched2", "IsL1MuMatched(nL1KBMTFSkimmed, idx2,L1KBMTFSkimmed_pt,L1KBMTFSkimmed_eta,L1KBMTFSkimmed_phi,nSkimmedL1Mu,SkimmedL1Mu_pt,SkimmedL1Mu_eta,SkimmedL1Mu_phi)")\
        .Define("nstub1", "GetNstub(nL1KBMTFSkimmed, idx1, L1KBMTFSkimmed_nStub)")\
        .Define("nstub2", "GetNstub(nL1KBMTFSkimmed, idx2, L1KBMTFSkimmed_nStub)")\
        .Define("pt1","my_mu1.Pt()").Define("eta1","my_mu1.Eta()").Define("phi1","my_mu1.Phi()") \
        .Define("beta1", "L1KBMTFSkimmed_beta[idx1]").Define("beta2", "L1KBMTFSkimmed_beta[idx2]") \
        .Define("charge1","L1KBMTFSkimmed_hwCharge[idx1]").Define("qual1","L1KBMTFSkimmed_hwQual[idx1]") \
        .Define("stub1Bx1","L1KBMTFSkimmed_s1Bx[idx1]").Define("stub2Bx1","L1KBMTFSkimmed_s2Bx[idx1]").Define("stub3Bx1","L1KBMTFSkimmed_s3Bx[idx1]").Define("stub4Bx1","L1KBMTFSkimmed_s4Bx[idx1]") \
        .Define("stub1Station1","L1KBMTFSkimmed_s1Station[idx1]").Define("stub2Station1","L1KBMTFSkimmed_s2Station[idx1]").Define("stub3Station1","L1KBMTFSkimmed_s3Station[idx1]").Define("stub4Station1","L1KBMTFSkimmed_s4Station[idx1]") \
        .Define("stub1Sector1","L1KBMTFSkimmed_s1Sector[idx1]").Define("stub2Sector1","L1KBMTFSkimmed_s2Sector[idx1]").Define("stub3Sector1","L1KBMTFSkimmed_s3Sector[idx1]").Define("stub4Sector1","L1KBMTFSkimmed_s4Sector[idx1]") \
        .Define("stub1Wheel1","L1KBMTFSkimmed_s1Wheel[idx1]").Define("stub2Wheel1","L1KBMTFSkimmed_s2Wheel[idx1]").Define("stub3Wheel1","L1KBMTFSkimmed_s3Wheel[idx1]").Define("stub4Wheel1","L1KBMTFSkimmed_s4Wheel[idx1]") \
        .Define("stub1Bx2","L1KBMTFSkimmed_s1Bx[idx2]").Define("stub2Bx2","L1KBMTFSkimmed_s2Bx[idx2]").Define("stub3Bx2","L1KBMTFSkimmed_s3Bx[idx2]").Define("stub4Bx2","L1KBMTFSkimmed_s4Bx[idx2]") \
        .Define("stub1Station2","L1KBMTFSkimmed_s1Station[idx2]").Define("stub2Station2","L1KBMTFSkimmed_s2Station[idx2]").Define("stub3Station2","L1KBMTFSkimmed_s3Station[idx2]").Define("stub4Station2","L1KBMTFSkimmed_s4Station[idx2]") \
        .Define("dxy1","L1KBMTFSkimmed_hwDXY[idx1]") \
        .Define("pt2","my_mu2.Pt()").Define("eta2","my_mu2.Eta()").Define("phi2","my_mu2.Phi()") \
        .Define("charge2","L1KBMTFSkimmed_hwCharge[idx2]").Define("qual2","L1KBMTFSkimmed_hwQual[idx2]") \
        .Define("dxy2","L1KBMTFSkimmed_hwDXY[idx2]") \
        .Define("hwK1","L1KBMTFSkimmed_hwK[idx1]") \
        .Define("hwK2","L1KBMTFSkimmed_hwK[idx2]") \
        .Define("phiB_S1","L1KBMTFSkimmed_s1PhiB[idx1]") \
        .Define("phiB_S2","L1KBMTFSkimmed_s2PhiB[idx1]") \
        .Define("phiB_S3","L1KBMTFSkimmed_s3PhiB[idx1]") \
        .Define("phiB_S4","L1KBMTFSkimmed_s4PhiB[idx1]") \
        

if (isdata):
    df = df.Define("met","L1KBMTFSkimmed_met_bx0[idx1]") 
else:
    df = df.Define("met","L1MET_pt[0]")

if isdata:
   df = df.Define("is_colliding", "IsColliding(run,bunchCrossing)")
   df = df.Define("is_earlier_colliding", "IsEarlierColliding(run,bunchCrossing,1,is_colliding)")
   df = df.Define("genbeta1","{}".format(weight)) \
          .Define("genpt1","{}".format(weight)) \
          .Define("genbeta2","{}".format(weight)) \
          .Define("genpt2","{}".format(weight))

else:
   df = df.Define("is_colliding", "true").Define("is_earlier_colliding", "true")
   df = df.Define("genbeta1","Get_genbeta(my_mu1.Eta(), my_mu1.Phi(), nGen, Gen_eta, Gen_phi, Gen_pdgid, Gen_beta)") \
          .Define("genpt1","Get_genbeta(my_mu1.Eta(), my_mu1.Phi(), nGen, Gen_eta, Gen_phi, Gen_pdgid, Gen_pt)") \
          .Define("genbeta2","Get_genbeta(my_mu2.Eta(), my_mu2.Phi(), nGen, Gen_eta, Gen_phi, Gen_pdgid, Gen_beta)") \
          .Define("genpt2","Get_genbeta(my_mu2.Eta(), my_mu2.Phi(), nGen, Gen_eta, Gen_phi, Gen_pdgid, Gen_pt)") \
        .Define("genK1","Get_genbeta(my_mu1.Eta(), my_mu1.Phi(), nGen, Gen_eta, Gen_phi, Gen_pdgid, Gen_K)") \
        .Define("genK2","Get_genbeta(my_mu2.Eta(), my_mu2.Phi(), nGen, Gen_eta, Gen_phi, Gen_pdgid, Gen_K)") \
        .Define("geneta1","Get_geneta(my_mu1.Eta(), my_mu1.Phi(), nGen, Gen_eta, Gen_phi, Gen_pdgid)") \
        .Define("geneta2","Get_geneta(my_mu2.Eta(), my_mu2.Phi(), nGen, Gen_eta, Gen_phi, Gen_pdgid)") \
       .Define("genCharge1","Get_genbeta(eta1, phi1, nGen, Gen_eta, Gen_phi, Gen_pdgid, Gen_charge)") \
       .Define("genCharge2","Get_genbeta(eta2, phi2, nGen, Gen_eta, Gen_phi, Gen_pdgid, Gen_charge)") \



columns = ROOT.std.vector("string")()
for c in ("run", "luminosityBlock", "bunchCrossing", "orbitNumber", \
        "mmumu", "isOS", "DRmumu", "xsweight", "met", "is_colliding", "is_earlier_colliding", \
        "bxspread1", "bxspread2", "isL1MuMatched1", "isL1MuMatched2", "nstub1", "nstub2",\
        "genbeta1","genpt1","genbeta2","genpt2", \
        "stub1Bx1", "stub2Bx1", "stub3Bx1", "stub4Bx1", "stub1Bx2", "stub2Bx2", "stub3Bx2", "stub4Bx2", \
        "stub1Station1", "stub2Station1", "stub3Station1", "stub4Station1", "stub1Station2", "stub2Station2", "stub3Station2", "stub4Station2", \
        "stub1Sector1", "stub2Sector1", "stub3Sector1", "stub4Sector1", "stub1Wheel1", "stub2Wheel1", "stub3Wheel1", "stub4Wheel1", \
        "pt1","eta1","phi1","charge1","qual1","dxy1","pt2","eta2","phi2","charge2","qual2","dxy2", "beta1", "beta2", "genK1", "genK2", "hwK1", "hwK2",\
         "phiB_S1", "phiB_S2", "phiB_S3", "phiB_S4", "geneta1", "geneta2", "genCharge1", "genCharge2"):
    columns.push_back(c)

#df.Snapshot("Events","/eos/cms/store/cmst3/group/taug2/AnalysisXuelong/ntuples_mutau_2018_basicsel/{}.root".format(sample),columns)
df.Snapshot("Events",output_file,columns)

nentries = df.Count().GetValue()
print ("After selection entries", nentries)

time_end=timer.time()
print('totally cost',time_end-time_start)
