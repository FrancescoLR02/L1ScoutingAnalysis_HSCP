#include <TH2.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include "TMultiGraph.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <stdio.h>
#include <TF1.h>
#include <TDirectoryFile.h>
#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TKey.h"
#include "THashList.h"
#include "THStack.h"
#include "TPaveLabel.h"
#include "TFile.h"
#include "TTree.h"
#include <TRandom3.h>
#include <algorithm>
#include "TMath.h"
#include "modzmumu_Tree.h"


// g++ -O3 CorrectionFactorsZmumu.cc -o CorrectionFactorsZmumu.exe $(root-config --cflags --glibs)

using namespace std;

// binning of the charge-odd curvature map
const int    NPHI      = 12;    // 30 deg bins
const int    NETA      = 5;     // over |eta| < 0.83
const int    NPHI_FINE = 12;    // 5 deg bins, for the charge count asymmetry
const double ETAMAX    = 0.83;
const double D = 1;

//maps phi into the correct geometrical sector
int phiBin(double phi, int nbins){
   const double w = 2*TMath::Pi()/nbins;
   double x = std::fmod(phi + 0.5*w, 2*TMath::Pi());
   if (x < 0) x += 2*TMath::Pi();
   int b = int(x / w);
   return std::min(std::max(b, 0), nbins-1);
}

int etaBin(double eta){
   int b = int((eta + ETAMAX) / (2*ETAMAX) * NETA);
   return std::min(std::max(b, 0), NETA-1);
}

//computes pT with the misalignment correction factor D
double Get_pTfromK(double K){
   int charge = (K >= 0) ? +1 : -1;
   const double lsb = 1.25 / float(1 << 13);
 
   double FK = fabs(K * lsb);

   //FK = .8569 * FK / (1.0 + 0.1144 * FK);

   double pt = 1/FK;
   if (pt > 5000) pt = 5000;

   return pt;
}

//Original implementation of the pT in the Kalman Filter
double ptLUT(double K) { 
  int charge = (K >= 0) ? +1 : -1;
  float lsb = 1.25 / float(1 << 13);
  double FK = fabs(K);

  if (FK > 2047) FK = 2047.; 
  if (FK < 9) FK = 9.; 

  FK = FK * lsb;
  //step 1 -material and B-field
  FK = .8569 * FK / (1.0 + 0.1144 * FK);
  //step 2 - misalignment
  FK = FK - charge * 1.23e-03;
  //Get to BMTF scale
  FK = FK / 1.17;

  double pt = 0;

  if (FK != 0) pt = 1 / FK;
  if (pt < 4) pt = 4;

  return pt;
}

//Delta values in each sector
const double DELTA_LUT[12] = {
    -4.517551982392351,  15.357480990344975,  16.890352083459184,
     4.568931798089349,  -8.445488610821378,  -3.2565913333840792,
     6.753763570952405,  17.763404130940582,  14.926302197139007,
     0.5210699539132264, -14.19312748508228,  -18.101331219109678
};

//Mow we apply eta dependent corrections! 
//Delta values in each eta bin, measured on the phi-corrected K 
const double DELTA_ETA_LUT[NETA] = {
   -5.031174684878439, -1.7392444643104927, 0.9836903473143631,
    2.591548000863815,  2.8348744112647695
};

double deltaKEta(double eta){ return DELTA_ETA_LUT[etaBin(eta)]; }


int deltaBin(double phi){ return phiBin(phi, 12); }
double deltaK(double phi){ return DELTA_LUT[deltaBin(phi)]; }


//Delta values vs number of stubs, measured on the phi+eta corrected K
const double DELTA_NSTUB_LUT[3] = { 0.970, -0.471, 0.084 };  

double deltaKNStub(int nstub){
   int k = nstub - 2;
   if (k < 0 || k > 2) return 0.;
   return DELTA_NSTUB_LUT[k];
}



int main(int argc, char** argv) {

   std::string input = *(argv + 1);
   std::string output = *(argv + 2);

   TFile *f_Double = new TFile(input.c_str());
   cout<<"XXXXXXXXXXXXX "<<input.c_str()<<" XXXXXXXXXXXX"<<endl;
   TTree *arbre1 = (TTree*) f_Double->Get("Events");


   /*arbre1->SetBranchAddress("run", &run);
   arbre1->SetBranchAddress("luminosityBlock", &luminosityBlock);
   arbre1->SetBranchAddress("bunchCrossing", &bunchCrossing);
   arbre1->SetBranchAddress("orbitNumber", &orbitNumber);*/

   arbre1->SetBranchAddress("mmumu", &mmumu);
   arbre1->SetBranchAddress("DRmumu", &DRmumu);
   arbre1->SetBranchAddress("xsweight", &xsweight);
   arbre1->SetBranchAddress("met", &met);
   arbre1->SetBranchAddress("bxspread1", &bxspread1);
   arbre1->SetBranchAddress("bxspread2", &bxspread2);
   arbre1->SetBranchAddress("isL1MuMatched1", &isL1MuMatched1);
   arbre1->SetBranchAddress("isL1MuMatched2", &isL1MuMatched2);
   arbre1->SetBranchAddress("nstub1", &nstub1);
   arbre1->SetBranchAddress("nstub2", &nstub2);
   arbre1->SetBranchAddress("pt1", &pt1);
   arbre1->SetBranchAddress("eta1", &eta1);
   arbre1->SetBranchAddress("phi1", &phi1);
   arbre1->SetBranchAddress("pt2", &pt2);
   arbre1->SetBranchAddress("eta2", &eta2);
   arbre1->SetBranchAddress("phi2", &phi2);
   arbre1->SetBranchAddress("stub1Bx1", &stub1Bx1);
   arbre1->SetBranchAddress("stub2Bx1", &stub2Bx1);
   arbre1->SetBranchAddress("stub3Bx1", &stub3Bx1);
   arbre1->SetBranchAddress("stub4Bx1", &stub4Bx1);
   arbre1->SetBranchAddress("stub1Bx2", &stub1Bx2);
   arbre1->SetBranchAddress("stub2Bx2", &stub2Bx2);
   arbre1->SetBranchAddress("stub3Bx2", &stub3Bx2);
   arbre1->SetBranchAddress("stub4Bx2", &stub4Bx2);
   arbre1->SetBranchAddress("hwK1", &hwK1);
   arbre1->SetBranchAddress("hwK2", &hwK2);
   arbre1->SetBranchAddress("beta1", &beta1);
   arbre1->SetBranchAddress("stub1Station1", &stub1Station1);
   arbre1->SetBranchAddress("stub2Station1", &stub2Station1);
   arbre1->SetBranchAddress("stub3Station1", &stub3Station1);
   arbre1->SetBranchAddress("stub4Station1", &stub4Station1);
   arbre1->SetBranchAddress("stub1Station2", &stub1Station2);
   arbre1->SetBranchAddress("stub2Station2", &stub2Station2);
   arbre1->SetBranchAddress("stub3Station2", &stub3Station2);
   arbre1->SetBranchAddress("stub4Station2", &stub4Station2);

   // charge/qual/dxy change scalar type from one ntuple version to the next
   // (dxy: Int_t in the 2024 skims, Double_t in simulation, Float_t in data),
   // so they are connected through ScalarBranch, which adapts to the type on file.
   charge1.connect(arbre1, "charge1");
   qual1.connect(arbre1, "qual1");
   dxy1.connect(arbre1, "dxy1");
   charge2.connect(arbre1, "charge2");
   qual2.connect(arbre1, "qual2");
   dxy2.connect(arbre1, "dxy2");
   

   const int    NK    = 2400;
   const double KMIN  = -400.0, KMAX = 400.0;  
   const double lsb = 1.25 / float(1 << 13);

   TH1F* h_mmumu_OS=new TH1F("h_mmumu_OS", "h_mmumu_OS", 50,50,160); h_mmumu_OS->Sumw2();
   TH1F* h_mmumu_SS=new TH1F("h_mmumu_SS", "h_mmumu_SS", 50,50,160); h_mmumu_SS->Sumw2();
   TH1F* h_mmumu_OS_corr = new TH1F("h_mmumu_OS_corr","h_mmumu_OS_corr", 50,50,160); h_mmumu_OS_corr->Sumw2();

   TH1F* h_pt_OS=new TH1F("h_pt_OS", "h_pt_OS", 13,0,260); h_pt_OS->Sumw2();
   TH1F* h_pt_SS=new TH1F("h_pt_SS", "h_pt_SS", 13,0,260); h_pt_SS->Sumw2();


   TH1F* h_K  = new TH1F("h_K", "hw curvature K", 1000, KMIN, KMAX); h_K->Sumw2();
   TH1F* h_lowK = new TH1F("h_lowK","h_lowK", NK/2, -150, 150); h_lowK->Sumw2();
   TH1F* h_D_pT = new TH1F("h_D_pT","h_D_pT", 300, 12.5, 1050); h_D_pT->Sumw2();
   TH1F* h_D_pT_plus = new TH1F("h_D_pT_plus","h_D_pT_plus", 300, 12.5, 1050); h_D_pT_plus->Sumw2();
   TH1F* h_D_pT_minus = new TH1F("h_D_pT_minus","h_D_pT_minus", 300, 12.5, 1050); h_D_pT_minus->Sumw2();
   
   TH1F* h_pTLUT = new TH1F("h_pTLUT","h_pTLUT", 300, 12.5, 1000); h_pTLUT->Sumw2();


   //! Phi-Eta dependent
   TH1F* h_K_plus_phieta[NPHI_FINE][NETA];
   TH1F* h_K_minus_phieta[NPHI_FINE][NETA];
   for (int i = 0; i < NPHI_FINE; ++i){
      for (int j = 0; j < NETA; ++j){
         h_K_plus_phieta[i][j] = new TH1F(Form("h_K_plus_phi%d_eta%d" ,i,j), "|K|, positive muons", 200, 0, 500);
         h_K_minus_phieta[i][j] = new TH1F(Form("h_K_minus_phi%d_eta%d",i,j), "|K|, negative muons", 200, 0, 500);
         h_K_plus_phieta[i][j]->Sumw2(); h_K_minus_phieta[i][j]->Sumw2();
      }
   }

   //! Phi-Eta dependent, after the Delta(phi) LUT correction
   TH1::AddDirectory(kFALSE); 
   TH1F* h_K_plus_phieta_corr[NPHI_FINE][NETA];
   TH1F* h_K_minus_phieta_corr[NPHI_FINE][NETA];
   for (int i = 0; i < NPHI_FINE; ++i){
      for (int j = 0; j < NETA; ++j){
         h_K_plus_phieta_corr[i][j]  = new TH1F(Form("h_K_plus_phi%d_eta%d" ,i,j), "|K| corrected, positive muons", 200, 0, 500);
         h_K_minus_phieta_corr[i][j] = new TH1F(Form("h_K_minus_phi%d_eta%d",i,j), "|K| corrected, negative muons", 200, 0, 500);
         // same names as the uncorrected map, so detach them from gDirectory to avoid clashes
         h_K_plus_phieta_corr[i][j]->SetDirectory(nullptr);
         h_K_minus_phieta_corr[i][j]->SetDirectory(nullptr);
         h_K_plus_phieta_corr[i][j]->Sumw2(); h_K_minus_phieta_corr[i][j]->Sumw2();
      }
   }

   // the LUT itself, for bookkeeping
   TH1F* h_deltaLUT = new TH1F("h_deltaLUT", "Delta(phi) LUT;phi bin;#Delta [LSB]", 12, -0.5, 11.5);
   h_deltaLUT->SetDirectory(nullptr);
   for (int k = 0; k < 12; ++k) h_deltaLUT->SetBinContent(k+1, DELTA_LUT[k]);


   //! Phi-Eta dependent, after the phi LUT and then the eta LUT
   TH1F* h_K_plus_phieta_corr2[NPHI_FINE][NETA];
   TH1F* h_K_minus_phieta_corr2[NPHI_FINE][NETA];
   for (int i = 0; i < NPHI_FINE; ++i){
      for (int j = 0; j < NETA; ++j){
         h_K_plus_phieta_corr2[i][j]  = new TH1F(Form("h_K_plus_phi%d_eta%d" ,i,j), "|K| phi+eta corrected, positive muons", 200, 0, 500);
         h_K_minus_phieta_corr2[i][j] = new TH1F(Form("h_K_minus_phi%d_eta%d",i,j), "|K| phi+eta corrected, negative muons", 200, 0, 500);
         h_K_plus_phieta_corr2[i][j]->Sumw2(); h_K_minus_phieta_corr2[i][j]->Sumw2();
      }
   } 

   TH1F* h_deltaEtaLUT = new TH1F("h_deltaEtaLUT", "Delta(eta) LUT after phi correction;eta bin;#Delta [LSB]", NETA, -0.5, NETA-0.5);
   for (int j = 0; j < NETA; ++j) h_deltaEtaLUT->SetBinContent(j+1, DELTA_ETA_LUT[j]);


   //! Phi-Eta dependent, after the phi, eta and nStub LUTs
   TH1F* h_K_plus_phieta_corr3[NPHI_FINE][NETA];
   TH1F* h_K_minus_phieta_corr3[NPHI_FINE][NETA];
   for (int i = 0; i < NPHI_FINE; ++i){
      for (int j = 0; j < NETA; ++j){
         h_K_plus_phieta_corr3[i][j]  = new TH1F(Form("h_K_plus_phi%d_eta%d" ,i,j), "|K| phi+eta+nStub corrected, positive muons", 200, 0, 500);
         h_K_minus_phieta_corr3[i][j] = new TH1F(Form("h_K_minus_phi%d_eta%d",i,j), "|K| phi+eta+nStub corrected, negative muons", 200, 0, 500);
         h_K_plus_phieta_corr3[i][j]->Sumw2(); h_K_minus_phieta_corr3[i][j]->Sumw2();
      }
   }

   //! nStub dependent, after all three corrections (for validation)
   TH1F* h_K_plus_nStub_corr[3];
   TH1F* h_K_minus_nStub_corr[3];
   for (int i = 0; i < 3; ++i){
      h_K_plus_nStub_corr[i]  = new TH1F(Form("h_K_plus_nStub%d" ,i+2), "|K| fully corrected, positive muons", 200, 0, 500);
      h_K_minus_nStub_corr[i] = new TH1F(Form("h_K_minus_nStub%d",i+2), "|K| fully corrected, negative muons", 200, 0, 500);
      h_K_plus_nStub_corr[i]->Sumw2(); h_K_minus_nStub_corr[i]->Sumw2();
   }

   TH1F* h_deltaNStubLUT = new TH1F("h_deltaNStubLUT", "Delta(nStub) LUT after phi+eta correction;nStub;#Delta [LSB]", 3, 1.5, 4.5);
   for (int k = 0; k < 3; ++k) h_deltaNStubLUT->SetBinContent(k+1, DELTA_NSTUB_LUT[k]);


   //! Phi dependent
   TH1F* h_K_plus_phi[NPHI];
   TH1F* h_K_minus_phi[NPHI];
   for (int i = 0; i < NPHI; ++i){
         h_K_plus_phi[i] = new TH1F(Form("h_K_plus_phi%d" ,i), "|K|, positive muons", 200, 0, 500);
         h_K_minus_phi[i] = new TH1F(Form("h_K_minus_phi%d",i), "|K|, negative muons", 200, 0, 500);
         h_K_plus_phi[i]->Sumw2(); h_K_minus_phi[i]->Sumw2();
   }

   //! nStub dependent
   TH1F* h_K_plus_nStub[3];
   TH1F* h_K_minus_nStub[3];
   for (int i = 0; i < 3; ++i){
         h_K_plus_nStub[i] = new TH1F(Form("h_K_plus_nStub%d" ,i+2), "|K|, positive muons", 200, 0, 500);
         h_K_minus_nStub[i] = new TH1F(Form("h_K_minus_nStub%d",i+2), "|K|, negative muons", 200, 0, 500);
         h_K_plus_nStub[i]->Sumw2(); h_K_minus_nStub[i]->Sumw2();
   }


   // the same spectra integrated over the map, for the count-above-threshold cross-check
   TH1F* h_K_plus_all  = new TH1F("h_K_plus_all" ,"|K|, positive muons", 400, 0, 150); h_K_plus_all->Sumw2();
   TH1F* h_K_minus_all = new TH1F("h_K_minus_all","|K|, negative muons", 400, 0, 150); h_K_minus_all->Sumw2();

   // charge count asymmetry vs phi, fine bins: this is the defect map, not a bias map
   TH1F* h_phi_plus  = new TH1F("h_phi_plus" ,"phi, positive muons", NPHI, -TMath::Pi(), TMath::Pi()); h_phi_plus->Sumw2();
   TH1F* h_phi_minus = new TH1F("h_phi_minus","phi, negative muons", NPHI, -TMath::Pi(), TMath::Pi()); h_phi_minus->Sumw2();
   TH1F* h_dxy  = new TH1F("h_dxy", "h_dxy", 100, 0, 1); h_dxy->Sumw2();
   TH1F* h_nstub  = new TH1F("h_nstub", "h_nstub", 3, 2, 5); h_nstub->Sumw2();


   TH1F* h_charge = new TH1F("h_charge","h_charge", 3, -1.5, 1.5); h_charge->Sumw2();
   TH1F* h_pt = new TH1F("h_pt","h_pt",300, 12.5, 1000); h_pt->Sumw2();

   TH1F* h_beta = new TH1F("h_beta", "h_beta", 50, 0, 1); h_beta->Sumw2();
   TH1F* misID_pt = new TH1F("misID_pt", "misID_pt", 100, 12.5, 1000); misID_pt->Sumw2();
   TH1F* misID_dxy = new TH1F("misID_dxy", "misID_dxy", 100, 0, 1); misID_dxy->Sumw2();
   TH1F* misID_nstub = new TH1F("misID_nstub", "misID_nstub", 3, 2, 5); misID_nstub->Sumw2();
   TH1F* misID_K = new TH1F("misID_K", "misID_K", 60, -0.04, 0.04); misID_K->Sumw2();
   TH1F* misID_invpT = new TH1F("misID_invpT", "misID_invpT", 60, -0.03, 0.03); misID_invpT->Sumw2();
   TH1F* misID_mmumuOS = new TH1F("misID_mmumuOS", "misID_mmumuOS", 50, 50, 160); misID_mmumuOS->Sumw2();
   TH1F* misID_mmumuSS = new TH1F("misID_mmumuSS", "misID_mmumuSS", 50, 50, 160); misID_mmumuSS->Sumw2();
   TH1F* misID_nStubProb = new TH1F("misID_nStubProb", "misID_nStubProb", 3, 2, 5); misID_nStubProb->Sumw2();




   static TRandom3 randGen(1234);
   float ptmin=0.0;
   float ptmax=10000.0;


   int nStub2_misID = 0; int nStub3_misID = 0; int nStub4_misID = 0;

   Int_t nentries_wtn = (Int_t) arbre1->GetEntries();
   for (Int_t i = 0; i < nentries_wtn; i++) {
   	arbre1->GetEntry(i);
      if (i % 100000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
      fflush(stdout);

      

      TLorentzVector my_mu1; my_mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.105);
      TLorentzVector my_mu2; my_mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.105);

      TLorentzVector my_mu1_corr; my_mu1_corr.SetPtEtaPhiM(Get_pTfromK(hwK1), eta1,phi1,0.105);
      TLorentzVector my_mu2_corr; my_mu2_corr.SetPtEtaPhiM(Get_pTfromK(hwK2), eta2,phi2,0.105);


      if (dxy1>=1 or dxy2>=1) continue;
      if (pt1<15 or pt2<15) continue;
      if (qual1<12 or qual2<12) continue;

      if (nstub1==4 and qual1<14) continue;
      if (nstub1==3 and qual1<13) continue;
      if (nstub1==2 and qual1<12) continue;
      if (nstub2==4 and qual2<14) continue;
      if (nstub2==3 and qual2<13) continue;
      if (nstub2==2 and qual2<12) continue;


      h_K->Fill(hwK1); h_K->Fill(hwK2);
      h_lowK->Fill(hwK1); h_lowK->Fill(hwK2);
      h_D_pT->Fill(Get_pTfromK(hwK1)); h_D_pT->Fill(Get_pTfromK(hwK2));
      h_pTLUT->Fill(ptLUT(hwK1)); h_pTLUT->Fill(ptLUT(hwK2));

      if(charge1 > 0) h_D_pT_plus->Fill(Get_pTfromK(hwK1));
      if(charge1 < 0) h_D_pT_minus->Fill(Get_pTfromK(hwK1));
      if(charge2 > 0) h_D_pT_plus->Fill(Get_pTfromK(hwK2));
      if(charge2 < 0) h_D_pT_minus->Fill(Get_pTfromK(hwK2));

      int i1f = phiBin(phi1, NPHI_FINE), i1 = phiBin(phi1, NPHI), j1 = etaBin(eta1);
      int i2f = phiBin(phi2, NPHI_FINE), i2 = phiBin(phi2, NPHI), j2 = etaBin(eta2);

      
      //Compute the corrected hwK with the LUT values defined before      
      double hwK1_lut = hwK1 - deltaK(phi1);
      double hwK2_lut = hwK2 - deltaK(phi2);
      
      if (charge1 > 0) h_K_plus_phieta_corr [i1f][j1]->Fill(fabs(hwK1_lut));
      else h_K_minus_phieta_corr[i1f][j1]->Fill(fabs(hwK1_lut));
      if (charge2 > 0) h_K_plus_phieta_corr [i2f][j2]->Fill(fabs(hwK2_lut));
      else h_K_minus_phieta_corr[i2f][j2]->Fill(fabs(hwK2_lut));


      //Compute the corrected hwK with the LUT values for ETA defined before
      double hwK1_lut2 = hwK1_lut - deltaKEta(eta1);
      double hwK2_lut2 = hwK2_lut - deltaKEta(eta2);

      if (charge1 > 0) h_K_plus_phieta_corr2 [i1f][j1]->Fill(fabs(hwK1_lut2));
      else             h_K_minus_phieta_corr2[i1f][j1]->Fill(fabs(hwK1_lut2));
      if (charge2 > 0) h_K_plus_phieta_corr2 [i2f][j2]->Fill(fabs(hwK2_lut2));
      else             h_K_minus_phieta_corr2[i2f][j2]->Fill(fabs(hwK2_lut2));
      


      if (charge1 > 0){ 
         h_K_plus_phieta [i1f][j1]->Fill(fabs(hwK1)); 
         h_K_plus_all ->Fill(fabs(hwK1)); h_phi_plus ->Fill(phi1); 
         h_K_plus_phi[i1]->Fill(fabs(hwK1));

         if(nstub1 == 2) h_K_plus_nStub[0]->Fill(fabs(hwK1_lut2));
         if(nstub1 == 3) h_K_plus_nStub[1]->Fill(fabs(hwK1_lut2));
         if(nstub1 == 4) h_K_plus_nStub[2]->Fill(fabs(hwK1_lut2));
      }
      else{ 
         h_K_minus_phieta[i1f][j1]->Fill(fabs(hwK1)); 
         h_K_minus_all->Fill(fabs(hwK1)); h_phi_minus->Fill(phi1); 
         h_K_minus_phi[i1]->Fill(fabs(hwK1));

         if(nstub1 == 2) h_K_minus_nStub[0]->Fill(fabs(hwK1_lut2));
         if(nstub1 == 3) h_K_minus_nStub[1]->Fill(fabs(hwK1_lut2));
         if(nstub1 == 4) h_K_minus_nStub[2]->Fill(fabs(hwK1_lut2));
      }

      if (charge2 > 0){ 
         h_K_plus_phieta [i2f][j2]->Fill(fabs(hwK2)); 
         h_K_plus_all ->Fill(fabs(hwK2)); h_phi_plus ->Fill(phi2); 
         h_K_plus_phi[i2]->Fill(fabs(hwK2));

         if(nstub2 == 2) h_K_plus_nStub[0]->Fill(fabs(hwK2_lut2));
         if(nstub2 == 3) h_K_plus_nStub[1]->Fill(fabs(hwK2_lut2));
         if(nstub2 == 4) h_K_plus_nStub[2]->Fill(fabs(hwK2_lut2));
      }
      else{
         h_K_minus_phieta[i2f][j2]->Fill(fabs(hwK2)); 
         h_K_minus_all->Fill(fabs(hwK2)); h_phi_minus->Fill(phi2); 
         h_K_minus_phi[i2]->Fill(fabs(hwK2));

         if(nstub2 == 2) h_K_minus_nStub[0]->Fill(fabs(hwK2_lut2));
         if(nstub2 == 3) h_K_minus_nStub[1]->Fill(fabs(hwK2_lut2));
         if(nstub2 == 4) h_K_minus_nStub[2]->Fill(fabs(hwK2_lut2));
      }


      //Apply the nStub correction on top of the phi+eta corrected K
      double hwK1_lut3 = hwK1_lut2 - deltaKNStub(int(nstub1));
      double hwK2_lut3 = hwK2_lut2 - deltaKNStub(int(nstub2));

      int n1 = int(nstub1) - 2, n2 = int(nstub2) - 2;

      if (charge1 > 0){
         h_K_plus_phieta_corr3[i1f][j1]->Fill(fabs(hwK1_lut3));
         if (n1 >= 0 && n1 < 3) h_K_plus_nStub_corr[n1]->Fill(fabs(hwK1_lut3));
      } else {
         h_K_minus_phieta_corr3[i1f][j1]->Fill(fabs(hwK1_lut3));
         if (n1 >= 0 && n1 < 3) h_K_minus_nStub_corr[n1]->Fill(fabs(hwK1_lut3));
      }
      if (charge2 > 0){
         h_K_plus_phieta_corr3[i2f][j2]->Fill(fabs(hwK2_lut3));
         if (n2 >= 0 && n2 < 3) h_K_plus_nStub_corr[n2]->Fill(fabs(hwK2_lut3));
      } else {
         h_K_minus_phieta_corr3[i2f][j2]->Fill(fabs(hwK2_lut3));
         if (n2 >= 0 && n2 < 3) h_K_minus_nStub_corr[n2]->Fill(fabs(hwK2_lut3));
      }
      
      h_charge->Fill(charge1); h_charge->Fill(charge2);
      h_pt->Fill(pt1); h_pt->Fill(pt2);
      h_dxy->Fill(dxy1); h_dxy->Fill(dxy2);
      h_nstub->Fill(nstub1); h_nstub->Fill(nstub2);
      h_mmumu_OS->Fill(mmumu);

      h_mmumu_OS_corr->Fill((my_mu1_corr + my_mu2_corr).M());


      if(beta1 < 0.95){
         h_beta->Fill(beta1);
         misID_pt->Fill(pt1); misID_pt->Fill(pt2);
         misID_dxy->Fill(dxy1); misID_dxy->Fill(dxy2);
         misID_nstub->Fill(nstub1); misID_nstub->Fill(nstub2);
         misID_K->Fill(hwK1*lsb); misID_K->Fill(hwK2*lsb);
         misID_invpT->Fill(charge1/pt1); misID_invpT->Fill(charge2/pt2);
         
         if(charge1*charge2<0) misID_mmumuOS->Fill((my_mu1_corr+my_mu2_corr).M());

      }
   } // end of loop over events

   

   TFile *fout = TFile::Open(output.c_str(), "RECREATE");
   fout->cd();

   TDirectory* dir1=fout->mkdir("Zmumu");
   dir1->cd();
   //h_mmumu_OS->SetName(name.c_str());
   h_mmumu_OS->Write();
   h_mmumu_OS_corr->Write();
   //h_mmumu_SS->Write();
   h_K->Write();
   h_lowK->Write();
   h_D_pT->Write();
   h_D_pT_plus->Write();
   h_D_pT_minus->Write();
   h_pTLUT->Write();

   h_charge->Write();
   h_pt->Write();
   //h_pt2->Write();
   h_dxy->Write();
   h_nstub->Write();

   h_K_plus_all->Write();
   h_K_minus_all->Write();
   h_phi_plus->Write();
   h_phi_minus->Write();

   TDirectory* dir3=fout->mkdir("KmapPhiEta");
   dir3->cd();
   for (int i = 0; i < NPHI_FINE; ++i){
      for (int j = 0; j < NETA; ++j){ h_K_plus_phieta[i][j]->Write(); h_K_minus_phieta[i][j]->Write(); }
   }

   TDirectory* dir4=fout->mkdir("KmapPhi");
   dir4->cd();
   for (int i = 0; i < NPHI; ++i){
      h_K_plus_phi[i]->Write(); h_K_minus_phi[i]->Write();
   }

   TDirectory* dir5=fout->mkdir("KmapnStub");
   dir5->cd();
   for (int i = 0; i < 3; ++i){
      h_K_plus_nStub[i]->Write(); h_K_minus_nStub[i]->Write();
   }

   TDirectory* dir6 = fout->mkdir("KmapPhiEtaCorr");
   dir6->cd();
   for (int i = 0; i < NPHI_FINE; ++i){
      for (int j = 0; j < NETA; ++j){ h_K_plus_phieta_corr[i][j]->Write(); h_K_minus_phieta_corr[i][j]->Write(); }
   }
   h_deltaLUT->Write();


   TDirectory* dir7 = fout->mkdir("KmapPhiEtaCorrEta");
   dir7->cd();
   for (int i = 0; i < NPHI_FINE; ++i){
      for (int j = 0; j < NETA; ++j){ h_K_plus_phieta_corr2[i][j]->Write(); h_K_minus_phieta_corr2[i][j]->Write(); }
   }
   h_deltaLUT->Write();
   h_deltaEtaLUT->Write();


   TDirectory* dir8 = fout->mkdir("KmapPhiEtaCorrEtaNStub");
   dir8->cd();
   for (int i = 0; i < NPHI_FINE; ++i){
      for (int j = 0; j < NETA; ++j){ h_K_plus_phieta_corr3[i][j]->Write(); h_K_minus_phieta_corr3[i][j]->Write(); }
   }
   for (int i = 0; i < 3; ++i){ h_K_plus_nStub_corr[i]->Write(); h_K_minus_nStub_corr[i]->Write(); }
   h_deltaLUT->Write();
   h_deltaEtaLUT->Write();
   h_deltaNStubLUT->Write();

   TDirectory* dir2=fout->mkdir("misID_BX");
   dir2->cd();
   h_beta->Write();
   misID_pt->Write();
   misID_dxy->Write();
   misID_nstub->Write();
   misID_K->Write();
   misID_invpT->Write();
   misID_mmumuOS->Write();



   fout->Close();
}

