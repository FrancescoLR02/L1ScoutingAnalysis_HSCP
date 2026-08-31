from ROOT import RDataFrame, TFile, TChain, TTree, TFile, TH1D, TLorentzVector,TCanvas
import numpy as np
import sys
from math import cos,sin,sqrt,pi
import ROOT
import time as timer
time_start=timer.time()

#TODO compile the script using        root -l -b -q -e '.L MODbasic_sel.cpp++'      in the lib folder

####---------------------- UP TO RUN 398392 ---------------------
ROOT.gInterpreter.AddIncludePath('/eos/user/f/flarover/HSCP_2025/HSCPanalysis/CMSSW_15_0_10/src/L1ScoutingAnalysisRDataFrame/NtupleAnalyzer/lib')
ROOT.gInterpreter.Declare('#include "MODbasic_sel.h"')
ROOT.gInterpreter.ProcessLine('.L /eos/user/f/flarover/HSCP_2025/HSCPanalysis/CMSSW_15_0_10/src/L1ScoutingAnalysisRDataFrame/NtupleAnalyzer/lib/MODbasic_sel.cpp+')
#ROOT.gSystem.Load('/eos/user/f/flarover/HSCP_2025/HSCPanalysis/CMSSW_15_0_10/src/L1ScoutingAnalysisRDataFrame/NtupleAnalyzer/lib/basic_sel_cpp.so')

### sample: output root file name
input_file = sys.argv[1]
output_file = sys.argv[2]

isdata = True
ngen=1.


#? Decide whether to use "Skimmed" or "Slow" (first or second pass in the chain)
Mode = 'Slow'

print ("isdata ", isdata)

weight = 1.0

if (isdata):
    weight = 1.0

df = RDataFrame(0)
df = RDataFrame("Events",input_file)
nentries = df.Count().GetValue()

print ("Before selection total entries", nentries)


df = df.Filter(f"nL1KBMTF{Mode}>0")


df = df.Define("idx1", f"GetIndex(1, nL1KBMTF{Mode}, L1KBMTF{Mode}_pt, L1KBMTF{Mode}_eta, L1KBMTF{Mode}_phi, L1KBMTF{Mode}_nStub, L1KBMTF{Mode}_s1Bx, L1KBMTF{Mode}_s2Bx, L1KBMTF{Mode}_s3Bx, L1KBMTF{Mode}_s4Bx)").Define("idx2", f"GetIndex(2, nL1KBMTF{Mode}, L1KBMTF{Mode}_pt, L1KBMTF{Mode}_eta, L1KBMTF{Mode}_phi, L1KBMTF{Mode}_nStub, L1KBMTF{Mode}_s1Bx, L1KBMTF{Mode}_s2Bx, L1KBMTF{Mode}_s3Bx, L1KBMTF{Mode}_s4Bx)")

df = df.Define("nstub1",f"L1KBMTF{Mode}_nStub[idx1]").Define("nstub2",f"L1KBMTF{Mode}_nStub[idx2]") 



df = df.Define("bxspread1", f"GetBxSpread(nL1KBMTF{Mode}, idx1, L1KBMTF{Mode}_nStub, L1KBMTF{Mode}_s1Bx, L1KBMTF{Mode}_s2Bx, L1KBMTF{Mode}_s3Bx, L1KBMTF{Mode}_s4Bx)") \
       .Define("bxspread2", f"GetBxSpread(nL1KBMTF{Mode}, idx2, L1KBMTF{Mode}_nStub, L1KBMTF{Mode}_s1Bx, L1KBMTF{Mode}_s2Bx, L1KBMTF{Mode}_s3Bx, L1KBMTF{Mode}_s4Bx)") \
       .Define("stationspread1", f"GetStationSpread(nL1KBMTF{Mode}, idx1, L1KBMTF{Mode}_nStub, L1KBMTF{Mode}_s1Station, L1KBMTF{Mode}_s2Station, L1KBMTF{Mode}_s3Station, L1KBMTF{Mode}_s4Station)") \
       .Define("stationspread2", f"GetStationSpread(nL1KBMTF{Mode}, idx2, L1KBMTF{Mode}_nStub, L1KBMTF{Mode}_s1Station, L1KBMTF{Mode}_s2Station, L1KBMTF{Mode}_s3Station, L1KBMTF{Mode}_s4Station)") \
       .Define("firstbx1", f"GetFirstBx(nL1KBMTF{Mode}, idx1, L1KBMTF{Mode}_nStub, L1KBMTF{Mode}_s1Bx, L1KBMTF{Mode}_s2Bx, L1KBMTF{Mode}_s3Bx, L1KBMTF{Mode}_s4Bx)") \
       .Define("firstbx2", f"GetFirstBx(nL1KBMTF{Mode}, idx2, L1KBMTF{Mode}_nStub, L1KBMTF{Mode}_s1Bx, L1KBMTF{Mode}_s2Bx, L1KBMTF{Mode}_s3Bx, L1KBMTF{Mode}_s4Bx)") \
       .Define("pt1",f"L1KBMTF{Mode}_pt[idx1]").Define("pt2",f"L1KBMTF{Mode}_pt[idx2]") \
       .Define("recobeta1", f"L1KBMTF{Mode}_beta[idx1]").Define("recobeta2", f"L1KBMTF{Mode}_beta[idx2]") \
       .Define("HwK1", f"L1KBMTF{Mode}_hwK[idx1]").Define("HwK2", f"L1KBMTF{Mode}_hwK[idx2]") \
       .Define("eta1",f"L1KBMTF{Mode}_eta[idx1]").Define("eta2",f"L1KBMTF{Mode}_eta[idx2]") \
       .Define("phi1",f"L1KBMTF{Mode}_phi[idx1]").Define("phi2",f"L1KBMTF{Mode}_phi[idx2]") \
       .Define("dxy1",f"L1KBMTF{Mode}_hwDXY[idx1]").Define("dxy2",f"L1KBMTF{Mode}_hwDXY[idx2]") \
       .Define("qual1",f"L1KBMTF{Mode}_hwQual[idx1]").Define("qual2",f"L1KBMTF{Mode}_hwQual[idx2]") \
       .Define("charge1",f"L1KBMTF{Mode}_hwCharge[idx1]").Define("charge2",f"L1KBMTF{Mode}_hwCharge[idx2]") \
       .Define("met_bx0", f"GetMET(0,firstbx1,L1KBMTF{Mode}_bx[idx1],L1KBMTF{Mode}_met_bx0[idx1],L1KBMTF{Mode}_met_bxm1[idx1],L1KBMTF{Mode}_met_bxm2[idx1],L1KBMTF{Mode}_met_bxm3[idx1],L1KBMTF{Mode}_met_bxm4[idx1],L1KBMTF{Mode}_met_bxm5[idx1],L1KBMTF{Mode}_met_bxm6[idx1],L1KBMTF{Mode}_met_bxm7[idx1],L1KBMTF{Mode}_met_bxm8[idx1],L1KBMTF{Mode}_met_bxm9[idx1])") \
       .Define("met_bxm1", f"GetMET(1,firstbx1,L1KBMTF{Mode}_bx[idx1],L1KBMTF{Mode}_met_bx0[idx1],L1KBMTF{Mode}_met_bxm1[idx1],L1KBMTF{Mode}_met_bxm2[idx1],L1KBMTF{Mode}_met_bxm3[idx1],L1KBMTF{Mode}_met_bxm4[idx1],L1KBMTF{Mode}_met_bxm5[idx1],L1KBMTF{Mode}_met_bxm6[idx1],L1KBMTF{Mode}_met_bxm7[idx1],L1KBMTF{Mode}_met_bxm8[idx1],L1KBMTF{Mode}_met_bxm9[idx1])") \
       .Define("met_bxm2", f"GetMET(2,firstbx1,L1KBMTF{Mode}_bx[idx1],L1KBMTF{Mode}_met_bx0[idx1],L1KBMTF{Mode}_met_bxm1[idx1],L1KBMTF{Mode}_met_bxm2[idx1],L1KBMTF{Mode}_met_bxm3[idx1],L1KBMTF{Mode}_met_bxm4[idx1],L1KBMTF{Mode}_met_bxm5[idx1],L1KBMTF{Mode}_met_bxm6[idx1],L1KBMTF{Mode}_met_bxm7[idx1],L1KBMTF{Mode}_met_bxm8[idx1],L1KBMTF{Mode}_met_bxm9[idx1])") \
       .Define("met_bxm3", f"GetMET(3,firstbx1,L1KBMTF{Mode}_bx[idx1],L1KBMTF{Mode}_met_bx0[idx1],L1KBMTF{Mode}_met_bxm1[idx1],L1KBMTF{Mode}_met_bxm2[idx1],L1KBMTF{Mode}_met_bxm3[idx1],L1KBMTF{Mode}_met_bxm4[idx1],L1KBMTF{Mode}_met_bxm5[idx1],L1KBMTF{Mode}_met_bxm6[idx1],L1KBMTF{Mode}_met_bxm7[idx1],L1KBMTF{Mode}_met_bxm8[idx1],L1KBMTF{Mode}_met_bxm9[idx1])") \
       .Define("met_bxm4", f"GetMET(4,firstbx1,L1KBMTF{Mode}_bx[idx1],L1KBMTF{Mode}_met_bx0[idx1],L1KBMTF{Mode}_met_bxm1[idx1],L1KBMTF{Mode}_met_bxm2[idx1],L1KBMTF{Mode}_met_bxm3[idx1],L1KBMTF{Mode}_met_bxm4[idx1],L1KBMTF{Mode}_met_bxm5[idx1],L1KBMTF{Mode}_met_bxm6[idx1],L1KBMTF{Mode}_met_bxm7[idx1],L1KBMTF{Mode}_met_bxm8[idx1],L1KBMTF{Mode}_met_bxm9[idx1])") \
       .Define("met_bxm5", f"GetMET(5,firstbx1,L1KBMTF{Mode}_bx[idx1],L1KBMTF{Mode}_met_bx0[idx1],L1KBMTF{Mode}_met_bxm1[idx1],L1KBMTF{Mode}_met_bxm2[idx1],L1KBMTF{Mode}_met_bxm3[idx1],L1KBMTF{Mode}_met_bxm4[idx1],L1KBMTF{Mode}_met_bxm5[idx1],L1KBMTF{Mode}_met_bxm6[idx1],L1KBMTF{Mode}_met_bxm7[idx1],L1KBMTF{Mode}_met_bxm8[idx1],L1KBMTF{Mode}_met_bxm9[idx1])")


#df = df.Filter("(bxspread1>0 || pt1>500) || (bxspread2>0 || pt2>500)")
#df = df.Filter("(bxspread1>0 || bxspread2>0)")

df = df.Define("is_colliding", "IsColliding(run,firstbx1)")

df = df.Define("is_earlier_colliding", "IsEarlierColliding(run,firstbx1,4,is_colliding)")

columns = ROOT.std.vector("string")()
for c in ("run", "luminosityBlock", "bunchCrossing", "orbitNumber", "is_colliding", "is_earlier_colliding", \
        #"nL1KBMTF{Mode}", "L1KBMTF{Mode}_hwCharge", "L1KBMTF{Mode}_hwQual", \
        "idx1", "idx2", \
        "bxspread1", "bxspread2", "stationspread1", "stationspread2", "nstub1", "nstub2", \
        "firstbx1","firstbx2", "pt1", "pt2", "eta1", "eta2", "phi1", "phi2", "dxy1", "dxy2", "qual1", "qual2", "charge1", "charge2", "recobeta1", "recobeta2", "HwK1", "HwK2", \
        "met_bx0","met_bxm1","met_bxm2","met_bxm3","met_bxm4","met_bxm5",):

    columns.push_back(c)

df.Snapshot("Events",output_file,columns)

nentries = df.Count().GetValue()
print ("After selection entries", nentries)

time_end=timer.time()
print('totally cost',time_end-time_start)
