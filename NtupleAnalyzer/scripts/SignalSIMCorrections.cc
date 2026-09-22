#include <iostream>
#include <string>
#include <stdio.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Correction_tr_Tree.h"
#include "KCorrection.h"
#include "SIM_KCorrection.h"


using namespace std;

//!  g++ -O3 SignalSIMCorrections.cc -o SignalSIMCorrections.exe $(root-config --cflags --glibs)
//!  ./SignalSIMCorrections.exe input.root output.root <sample_name>

bool passSelection(int qual, int nstub, double dxy, double K){
   if (dxy >= 1)                return false;
   if (qual < 12)               return false;
   if (nstub == 4 && qual < 14) return false;
   if (nstub == 3 && qual < 13) return false;
   //if (std::fabs(K) > 1000)     return false;
   if (nstub <= 2)              return false;
   return true;
}

const int    NGENPT    = 150;      // 20 GeV bins
const double GENPT_MIN = 0.;
const double GENPT_MAX = 3000.;


const int NSTAGE = 3;
const char* STAGE_NAME[NSTAGE] = {"an_2024", "an_2025", "LUT_applied"};

int main(int argc, char** argv) {

   std::string input  = argv[1];
   std::string output = argv[2];
   std::string name   = argv[3];

   TFile *fin = TFile::Open(input.c_str(), "READ");
   TTree *arbre = (TTree*) fin->Get("Events");

   // K branch name differs between the DY and HSCP skims
   //const char* kname = arbre->GetBranch("hwK1") ? "hwK1" : "HwK1";
   //if (!arbre->GetBranch(kname)) { cerr << "No hwK1/HwK1 branch in " << input << endl; return 1; }

   // ScalarBranch reads whatever type is on file (Short_t/Int_t, Float_t/Double_t)
   ScalarBranch b_K, b_nstub, b_genpt, b_charge, b_qual, b_eta, b_phi, b_pt, b_dxy, b_ntrk;
   ScalarBranch b_K2, b_nstub2, b_genpt2, b_charge2, b_qual2, b_eta2, b_phi2, b_pt2, b_dxy2, b_ntrk;
   b_K.connect(arbre, kname);
   b_nstub.connect(arbre, "nstub1");
   b_genpt.connect(arbre, "genpt1");
   b_charge.connect(arbre, "charge1");
   b_qual.connect(arbre, "qual1");
   b_eta.connect(arbre, "eta1");
   b_phi.connect(arbre, "phi1");
   b_pt.connect(arbre, "pt1");
   b_dxy.connect(arbre, "dxy1");
   arbre->SetBranchAddress("genCharge1", &genCharge1);

   b_K2.connect(arbre, kname);
   b_nstub.connect(arbre, "nstub1");
   b_genpt.connect(arbre, "genpt1");
   b_charge.connect(arbre, "charge1");
   b_qual.connect(arbre, "qual1");
   b_eta.connect(arbre, "eta1");
   b_phi.connect(arbre, "phi1");
   b_pt.connect(arbre, "pt1");
   b_dxy.connect(arbre, "dxy1");
   arbre->SetBranchAddress("genCharge1", &genCharge1);
   bool hasNTrk = arbre->GetBranch("nL1KBMTFSkimmed");
   if (hasNTrk) b_ntrk.connect(arbre, "nL1KBMTFSkimmed");

   TH1::AddDirectory(kFALSE);
   //TH1::SetDefaultSumw2(kTRUE);

   TH1D* h_den = new TH1D("h_den_genpt", ";gen p_{T} [GeV];tracks", NGENPT, GENPT_MIN, GENPT_MAX);
   TH1D* h_genK = new TH1D("h_genK", "h_genK", 1000, -80, 80);

   TH1D* h_hwK[NSTAGE], *h_pt[NSTAGE], *h_resp[NSTAGE], *h_Kresp[NSTAGE];


   const char* RLAB = "p_{T}^{L1}/p_{T}^{gen} - 1";
   const char* RKLAB = "K^{reco}/K^{gen} - 1";
   for (int s = 0; s < NSTAGE; ++s){
      h_hwK[s]         = new TH1D("h_hwK", ";signed hwK [LSB];tracks", 400, -80, 80);
      h_pt[s]          = new TH1D("h_pt", ";L1 p_{T} [GeV];tracks", 300, 0, 1500);
      h_resp[s]        = new TH1D("h_resp", Form(";%s;tracks", RLAB), 200, -5, 10);
      h_Kresp[s]        = new TH1D("h_Kresp", Form(";%s;Ktracks", RKLAB), 200, -5, 10);
   }

   long nSel = 0;

   Long64_t nentries = arbre->GetEntries();
   for (Long64_t ev = 0; ev < nentries; ++ev) {
      arbre->GetEntry(ev);
      if (ev % 10000 == 0) { fprintf(stdout, "\r  Processed events: %8lld of %8lld ", ev, nentries); fflush(stdout); }

      double K      = double(b_K);
      double phi    = double(b_phi);
      double eta    = double(b_eta);
      double dxy    = double(b_dxy);
      double genpt  = double(b_genpt);
      double ptreco = double(b_pt);
      int    nstub  = int(b_nstub);
      int    qual   = int(b_qual);
      int    charge = int(b_charge);

      //bool valid = hasNTrk ? (int(b_ntrk) >= 1) : (nstub >= 2);
      if (genpt <= 0 || ptreco < 0) continue;
      if (!passSelection(qual, nstub, dxy, K)) continue;
      ++nSel;

      double Kcorr = sim_kcorr::correctK(K, phi, eta, nstub);

      double Ks[NSTAGE] = {K, K, Kcorr};
      double pt[NSTAGE];
      pt[0] = sim_kcorr::Get_newpt(K);
      pt[1] = sim_kcorr::ptLUT_corr(K, true);
      pt[2] = sim_kcorr::ptLUT_corr(Kcorr, true);

      h_den->Fill(genpt);
      double Kgen = genCharge1/(genpt)/sim_kcorr::LSB;
      h_genK->Fill(Kgen);

      for (int s = 0; s < NSTAGE; ++s) {
         double r = pt[s] / genpt - 1.;
         double Kresp = Ks[s] / Kgen -1.;
         h_Kresp[s]->Fill(Kresp);
         h_hwK[s]->Fill(Ks[s]);
         h_pt[s]->Fill(pt[s]);
         h_resp[s]->Fill(r);
      }
   }

   printf("\n  selected tracks: %ld\n", nSel);

   TFile *fout = TFile::Open(output.c_str(), "RECREATE");
   TDirectory* top = fout->mkdir(name.c_str());

   top->mkdir("common")->cd();
   h_den->Write();
   h_genK->Write();

   for (int s = 0; s < NSTAGE; ++s) {
      top->mkdir(STAGE_NAME[s])->cd();
      h_hwK[s]->Write();
      h_Kresp[s]->Write();
      h_pt[s]->Write();
      h_resp[s]->Write();
   }

   fout->Close();
   return 0;
}