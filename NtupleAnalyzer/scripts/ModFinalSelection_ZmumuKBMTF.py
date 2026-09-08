from ROOT import RDataFrame, TFile, TChain, TTree, TFile, TH1D, TLorentzVector,TCanvas
import numpy as np
import sys
from math import cos,sin,sqrt,pi
import ROOT
import time as timer
time_start=timer.time()
ROOT.gInterpreter.AddIncludePath('/eos/user/f/flarover/HSCP_2025/HSCPanalysis/CMSSW_15_0_10/src/L1ScoutingAnalysisRDataFrame/NtupleAnalyzer/lib')

sim = False

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

if sim: mode = 'Skimmed'
else: mode = 'Slow'

df = df.Filter(f"nL1KBMTF{mode}>1")

#df_var=df.Define("idx1",f"GetIndex_nostub_hwK(1,nL1KBMTF{mode}, L1KBMTF{mode}_hwK, L1KBMTF{mode}_eta, L1KBMTF{mode}_phi)").Define("idx2",f"GetIndex_nostub_hwK(2,nL1KBMTF{mode}, L1KBMTF{mode}_hwK, L1KBMTF{mode}_eta, L1KBMTF{mode}_phi)")
df_var=df.Define("idx1",f"GetIndexByCharge(+1,nL1KBMTF{mode}, L1KBMTF{mode}_hwK, L1KBMTF{mode}_eta, L1KBMTF{mode}_phi, L1KBMTF{mode}_nStub)").Define("idx2",f"GetIndexByCharge(-1,nL1KBMTF{mode}, L1KBMTF{mode}_hwK, L1KBMTF{mode}_eta, L1KBMTF{mode}_phi, L1KBMTF{mode}_nStub)")

#df_var=df_var.Define("my_mu1",f"GetLepVector_hwK(idx1,L1KBMTF{mode}_eta,L1KBMTF{mode}_phi,L1KBMTF{mode}_hwK)").Define("my_mu2",f"GetLepVector_hwK(idx2,L1KBMTF{mode}_eta,L1KBMTF{mode}_phi,L1KBMTF{mode}_hwK)").Define("isOS",f"L1KBMTF{mode}_hwCharge[idx1]*L1KBMTF{mode}_hwCharge[idx2]<0")
df_var=df_var.Define("my_mu1",f"GetLepVector_hwK(idx1,L1KBMTF{mode}_eta,L1KBMTF{mode}_phi,L1KBMTF{mode}_pt)").Define("my_mu2",f"GetLepVector_hwK(idx2,L1KBMTF{mode}_eta,L1KBMTF{mode}_phi,L1KBMTF{mode}_pt)").Define("isOS",f"L1KBMTF{mode}_hwCharge[idx1]*L1KBMTF{mode}_hwCharge[idx2]<0")

df_var = df_var.Filter("idx1!=idx2 && my_mu1.Pt()>20 && my_mu2.Pt()>20 && fabs(my_mu1.Eta())<0.83 && fabs(my_mu2.Eta())<0.83")

df_var = df_var.Define("mmumu","(my_mu1+my_mu2).M()").Define("DRmumu","my_mu1.DeltaR(my_mu2)")
df_var = df_var.Filter("fabs(my_mu1.Eta()-my_mu2.Eta())>0.05")


df = df_var.Filter("mmumu>50 && DRmumu>0.3").Define("xsweight","{}".format(weight))


df = df.Define("bxspread1", f"GetBxSpread(nL1KBMTF{mode}, idx1, L1KBMTF{mode}_nStub, L1KBMTF{mode}_s1Bx, L1KBMTF{mode}_s2Bx, L1KBMTF{mode}_s3Bx, L1KBMTF{mode}_s4Bx)")\
        .Define("bxspread2", f"GetBxSpread(nL1KBMTF{mode}, idx2, L1KBMTF{mode}_nStub, L1KBMTF{mode}_s1Bx, L1KBMTF{mode}_s2Bx, L1KBMTF{mode}_s3Bx, L1KBMTF{mode}_s4Bx)")\
        .Define("nstub1", f"GetNstub(nL1KBMTF{mode}, idx1, L1KBMTF{mode}_nStub)")\
        .Define("nstub2", f"GetNstub(nL1KBMTF{mode}, idx2, L1KBMTF{mode}_nStub)")\
        .Define("new_pT1", f"ptLUT(L1KBMTF{mode}_hwK[idx1], L1KBMTF{mode}_nStub[idx1], L1KBMTF{mode}_hwCharge[idx1])")\
        .Define("new_pT2", f"ptLUT(L1KBMTF{mode}_hwK[idx2], L1KBMTF{mode}_nStub[idx2], L1KBMTF{mode}_hwCharge[idx2])")\
        .Define("pt1","my_mu1.Pt()").Define("eta1","my_mu1.Eta()").Define("phi1","my_mu1.Phi()") \
        .Define("beta1", f"L1KBMTF{mode}_beta[idx1]").Define("beta2", f"L1KBMTF{mode}_beta[idx2]") \
        .Define("charge1",f"L1KBMTF{mode}_hwCharge[idx1]").Define("qual1",f"L1KBMTF{mode}_hwQual[idx1]") \
        .Define("stub1Bx1",f"L1KBMTF{mode}_s1Bx[idx1]").Define("stub2Bx1",f"L1KBMTF{mode}_s2Bx[idx1]").Define("stub3Bx1",f"L1KBMTF{mode}_s3Bx[idx1]").Define("stub4Bx1",f"L1KBMTF{mode}_s4Bx[idx1]") \
        .Define("stub1Station1",f"L1KBMTF{mode}_s1Station[idx1]").Define("stub2Station1",f"L1KBMTF{mode}_s2Station[idx1]").Define("stub3Station1",f"L1KBMTF{mode}_s3Station[idx1]").Define("stub4Station1",f"L1KBMTF{mode}_s4Station[idx1]") \
        .Define("stub1Sector1",f"L1KBMTF{mode}_s1Sector[idx1]").Define("stub2Sector1",f"L1KBMTF{mode}_s2Sector[idx1]").Define("stub3Sector1",f"L1KBMTF{mode}_s3Sector[idx1]").Define("stub4Sector1",f"L1KBMTF{mode}_s4Sector[idx1]") \
        .Define("stub1Wheel1",f"L1KBMTF{mode}_s1Wheel[idx1]").Define("stub2Wheel1",f"L1KBMTF{mode}_s2Wheel[idx1]").Define("stub3Wheel1",f"L1KBMTF{mode}_s3Wheel[idx1]").Define("stub4Wheel1",f"L1KBMTF{mode}_s4Wheel[idx1]") \
        .Define("stub1Bx2",f"L1KBMTF{mode}_s1Bx[idx2]").Define("stub2Bx2",f"L1KBMTF{mode}_s2Bx[idx2]").Define("stub3Bx2",f"L1KBMTF{mode}_s3Bx[idx2]").Define("stub4Bx2",f"L1KBMTF{mode}_s4Bx[idx2]") \
        .Define("stub1Station2",f"L1KBMTF{mode}_s1Station[idx2]").Define("stub2Station2",f"L1KBMTF{mode}_s2Station[idx2]").Define("stub3Station2",f"L1KBMTF{mode}_s3Station[idx2]").Define("stub4Station2",f"L1KBMTF{mode}_s4Station[idx2]") \
        .Define("dxy1",f"L1KBMTF{mode}_hwDXY[idx1]") \
        .Define("pt2","my_mu2.Pt()").Define("eta2","my_mu2.Eta()").Define("phi2","my_mu2.Phi()") \
        .Define("charge2",f"L1KBMTF{mode}_hwCharge[idx2]").Define("qual2",f"L1KBMTF{mode}_hwQual[idx2]") \
        .Define("dxy2",f"L1KBMTF{mode}_hwDXY[idx2]") \
        .Define("hwK1",f"L1KBMTF{mode}_hwK[idx1]") \
        .Define("hwK2",f"L1KBMTF{mode}_hwK[idx2]") \
        .Define("ngen","{}".format(nentries))
        

if (isdata):
    df = df.Define("met",f"L1KBMTF{mode}_met_bx0[idx1]") 
else:
    df = df.Define("met","L1MET_pt[0]")

if isdata:
   df = df.Define("is_colliding", "IsColliding(run,bunchCrossing)")
   df = df.Define("is_earlier_colliding", "IsEarlierColliding(run,bunchCrossing,1,is_colliding)")
   df = df.Define("genbeta1","{}".format(weight)) \
          .Define("genpt1","{}".format(weight)) \
          .Define("genbeta2","{}".format(weight)) \
          .Define("genpt2","{}".format(weight)) \
          .Define("genCharge1","{}".format(weight)) \
          .Define("genCharge2","{}".format(weight)) \
          .Define("geneta2","{}".format(weight)) \
          .Define("geneta1","{}".format(weight)) \
          .Define("genK1","{}".format(weight)) \
          .Define("genK2","{}".format(weight)) \

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
        "bxspread1", "bxspread2", "nstub1", "nstub2",\
        "genbeta1","genpt1","genbeta2","genpt2", \
        "stub1Bx1", "stub2Bx1", "stub3Bx1", "stub4Bx1", "stub1Bx2", "stub2Bx2", "stub3Bx2", "stub4Bx2", \
        "stub1Station1", "stub2Station1", "stub3Station1", "stub4Station1", "stub1Station2", "stub2Station2", "stub3Station2", "stub4Station2", \
        "stub1Sector1", "stub2Sector1", "stub3Sector1", "stub4Sector1", "stub1Wheel1", "stub2Wheel1", "stub3Wheel1", "stub4Wheel1", \
        "pt1","eta1","phi1","charge1","qual1","dxy1","pt2","eta2","phi2","charge2","qual2","dxy2", "beta1", "beta2", "genK1", "genK2", "hwK1", "hwK2",\
         "geneta1", "geneta2", "genCharge1", "genCharge2", "ngen", "new_pT1", "new_pT2"):
    columns.push_back(c)

#df.Snapshot("Events","/eos/cms/store/cmst3/group/taug2/AnalysisXuelong/ntuples_mutau_2018_basicsel/{}.root".format(sample),columns)
df.Snapshot("Events",output_file,columns)

nentries = df.Count().GetValue()
print ("After selection entries", nentries)

time_end=timer.time()
print('totally cost',time_end-time_start)
