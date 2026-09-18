#include <iostream>
#include <string>
#include <stdio.h>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "Correction_tr_Tree.h"
#include "KCorrection.h"

using namespace std;

//!  g++ -O3 DataCorrectionsCheck.cc -o DataCorrectionsCheck.exe $(root-config --cflags --glibs)
//!  ./DataCorrectionsCheck.exe input.root output.root <sample_name> [all|even|odd]

const double PT_MIN = 15.;   // same pT cut used to derive the LUTs
const double KCUT = 700;

bool passSelection(int qual, int nstub, double dxy, double ptreco){
   if (nstub <= 2)              return false;
   //if (ptreco < PT_MIN)         return false;
   if (dxy >= 1)                return false;
   if (qual < 12)               return false;
   if (nstub == 4 && qual < 14) return false;
   if (nstub == 3 && qual < 13) return false;
   return true;
}

const int    NTHR = 5;
const double PT_THR[NTHR] = {5., 20., 22., 50., 100.};

const int NSTAGE = 3;
const char* STAGE_NAME[NSTAGE] = {"orig", "an_2024", "corr"};

const int NKV = 2;
const char* KV_NAME[NKV] = {"KmapRaw", "KmapCorr"};
const char* QN[2] = {"plus", "minus"};

// full hwK range, so no track ends up in the overflow
const int    NKBIN = 700;
const double KMAX  = 700.;

// ---- model-dependent categories
const int NCAT = 14;
const char* CAT_NAME[NCAT] = {
   "stub4_bx124", "stub4_bx123", "stub4_bx1234",
   "stub4_bx1222_2tracks", "stub4_bx1122_2tracks", "stub4_bx1112_2tracks",
   "stub3_bx124_slow", "stub3_bx124_fast", "stub3_bx123_slow", "stub3_bx123_fast",
   "stub3_bx122_slow_2tracks", "stub3_bx122_fast_2tracks",
   "stub3_bx112_slow_2tracks", "stub3_bx112_fast_2tracks"
};

// true = the category also requires a second good track
const bool NEEDS_2TRK[NCAT] = {
   false, false, false,
   true,  true,  true,
   false, false, false, false,
   true,  true,  true,  true
};

// returns the category index, or -1 if the track is in none
int category(int nstub, int bx, bool slow){
   if (nstub == 4){
      if (bx==3100 || bx==3110 || bx==3310 || bx==3200 || bx==3220 || bx==3320) return 0;  // bx124
      if (bx==2100 || bx==2110 || bx==2210) return 1;                                     // bx123
      if (bx==3210) return 2;                                                             // bx1234
      if (bx==1110) return 3;                                                             // bx1222
      if (bx==1100) return 4;                                                             // bx1122
      if (bx==1000) return 5;                                                             // bx1112
   }
   if (nstub == 3){
      int k = -1;
      if      (bx==3200 || bx==3100) k = 6;    // bx124
      else if (bx==2100)             k = 8;    // bx123
      else if (bx==1100)             k = 10;   // bx122
      else if (bx==1000)             k = 12;   // bx112
      if (k >= 0) return slow ? k : k + 1;
   }
   return -1;
}

double catPtMin(int c){ return (c == 2) ? 15. : 50.; }   // bx1234 has a looser cut

const int NMOD = 2;
const char* MOD_DIR[NMOD]   = {"model_an2024", "model_corr"};
const int   MOD_STAGE[NMOD] = {1, 2};                     // index in ptS[]

int main(int argc, char** argv) {

   std::string input  = argv[1];
   std::string output = argv[2];
   std::string name   = argv[3];
   std::string parity = (argc > 4) ? argv[4] : "all";

   TFile *fin = TFile::Open(input.c_str(), "READ");
   TTree *arbre = (TTree*) fin->Get("Events");

   const char* kname = arbre->GetBranch("hwK1") ? "hwK1" : "HwK1";
   if (!arbre->GetBranch(kname)) { cerr << "No hwK1/HwK1 branch" << endl; return 1; }

   // ---- track 1
   ScalarBranch b_K, b_nstub, b_charge, b_qual, b_eta, b_phi, b_pt, b_dxy, b_ntrk, b_bxspread, b_stationspread;
   b_K.connect(arbre, kname);
   b_nstub.connect(arbre, "nstub1");
   b_charge.connect(arbre, "charge1");
   b_qual.connect(arbre, "qual1");
   b_eta.connect(arbre, "eta1");
   b_phi.connect(arbre, "phi1");
   b_pt.connect(arbre, "pt1");
   b_dxy.connect(arbre, "dxy1");
   b_bxspread.connect(arbre, "bxspread1");
   b_stationspread.connect(arbre, "stationspread1");
   bool hasNTrk = arbre->GetBranch("nL1KBMTFSkimmed");
   if (hasNTrk) b_ntrk.connect(arbre, "nL1KBMTFSkimmed");

   // ---- track 2, used only for the 2-track categories
   ScalarBranch b_idx2, b_pt2, b_dxy2, b_qual2, b_nstub2;
   bool hasTrk2 = arbre->GetBranch("idx2");
   if (hasTrk2) {
      b_idx2.connect(arbre, "idx2");
      b_pt2.connect(arbre, "pt2");
      b_dxy2.connect(arbre, "dxy2");
      b_qual2.connect(arbre, "qual2");
      b_nstub2.connect(arbre, "nstub2");
   } else {
      cerr << "No idx2 branch: the _2tracks categories will stay empty." << endl;
   }

   printf("  K branch: %s   events used: %s\n", kname, parity.c_str());

   TH1::AddDirectory(kFALSE);
   TH1::SetDefaultSumw2(kTRUE);

   // ---- charge-odd |K| maps, raw and corrected K
   TH1D* h_Kmap[NKV][2][kcorr::NSEC][kcorr::NETA];
   TH1D* h_Kns[NKV][2][3];
   TH1D* h_Kall[NKV][2];
   for (int v = 0; v < NKV; ++v)
      for (int q = 0; q < 2; ++q){
         for (int i = 0; i < kcorr::NSEC; ++i)
            for (int j = 0; j < kcorr::NETA; ++j)
               h_Kmap[v][q][i][j] = new TH1D(Form("h_K_%s_phi%d_eta%d", QN[q], i, j), "|hwK| [LSB]", NKBIN, 0, KMAX);
         for (int n = 0; n < 3; ++n)
            h_Kns[v][q][n] = new TH1D(Form("h_K_%s_nStub%d", QN[q], n+2), "|hwK| [LSB]", NKBIN, 0, KMAX);
         h_Kall[v][q] = new TH1D(Form("h_K_%s_all", QN[q]), "|hwK| [LSB]", NKBIN, 0, KMAX);
      }

   // ---- per pT stage
   TH1D* h_K[NSTAGE], *h_pt[NSTAGE], *h_ptq[NSTAGE][2], *h_npass[NSTAGE];
   for (int s = 0; s < NSTAGE; ++s){
      h_K[s]  = new TH1D("h_K",  "signed hwK [LSB]", 2*NKBIN, -KMAX, KMAX);
      h_pt[s] = new TH1D("h_pt", "L1 p_{T} [GeV]", 150, 0, 1050);
      for (int q = 0; q < 2; ++q)
         h_ptq[s][q] = new TH1D(Form("h_pt_%s", QN[q]), "L1 p_{T} [GeV]", 150, 0, 1050);
      h_npass[s] = new TH1D("h_npass", "threshold", NTHR, -0.5, NTHR - 0.5);
      for (int k = 0; k < NTHR; ++k) h_npass[s]->GetXaxis()->SetBinLabel(k+1, Form("p_{T} > %.0f", PT_THR[k]));
   }

   TH1D* h_flip = new TH1D("h_flip_pt", "charge sign changed by the correction;L1 p_{T} (orig) [GeV]", 150, 0, 1050);
   TH1D* h_charge = new TH1D("h_charge", "h_charge", 3, -1, 2);

   // ---- pT per category, for an_2024 and corrected pT
   TH1D* h_cat[NMOD][NCAT];
   for (int m = 0; m < NMOD; ++m)
      for (int c = 0; c < NCAT; ++c) {
         if (c == 2) h_cat[m][c] = new TH1D(CAT_NAME[c], ";L1 p_{T} [GeV]", 20, 15, 100);
         else        h_cat[m][c] = new TH1D(CAT_NAME[c], ";L1 p_{T} [GeV]", 100, 50, 1050);
      }

   long nTrk = 0, nSignOK = 0, nFlip = 0;

   Long64_t nentries = arbre->GetEntries();
   for (Long64_t ev = 0; ev < nentries; ++ev) {
      if (parity == "even" && ev % 2 != 0) continue;
      if (parity == "odd"  && ev % 2 == 0) continue;

      arbre->GetEntry(ev);
      if (ev % 100000 == 0) { fprintf(stdout, "\r  Processed events: %8lld of %8lld ", ev, nentries); fflush(stdout); }

      if (hasNTrk && int(b_ntrk) < 1) continue;

      double K             = double(b_K);
      double eta           = double(b_eta);
      double phi           = double(b_phi);
      double ptreco        = double(b_pt);
      double dxy           = double(b_dxy);
      int    nstub         = int(b_nstub);
      int    qual          = int(b_qual);
      int    charge        = int(b_charge);
      int    bxspread      = int(b_bxspread);
      int    stationspread = int(b_stationspread);

      if (!passSelection(qual, nstub, dxy, ptreco)) continue;
      ++nTrk;
      if ((K >= 0) == (charge > 0)) ++nSignOK;

      h_charge->Fill(charge);

      int    Kint  = static_cast<int>(K);
      double Kcorr = kcorr::correctK(K, phi, eta, nstub);

      double Ks[NSTAGE]  = {double(Kint), double(Kint), Kcorr};
      double ptS[NSTAGE] = {kcorr::ptLUT_orig(Kint),
                            kcorr::Get_newpt(Kint),
                            kcorr::ptLUT_corr(Kcorr, true)};

      int sec = kcorr::sectorBin(phi);
      int je  = kcorr::etaBin(eta);
      int ns  = nstub - 2;

      // charge labels: raw K keeps the ntuple charge, corrected K uses its own sign
      int q     = (charge > 0) ? 0 : 1;
      int qcorr = (kcorr::chargeFromK(Kcorr) > 0) ? 0 : 1;

      // charge-odd maps: raw K and fully corrected K
      //double Kv[NKV] = {K, kcorr::applyMaterialMap(Kcorr)};
      double Kv[NKV] = {K, Kcorr};
      int    Qv[NKV] = {q, qcorr};
      for (int v = 0; v < NKV; ++v) {
         if (std::fabs(Kv[v]) > KCUT) continue;
         h_Kmap[v][Qv[v]][sec][je]->Fill(fabs(Kv[v]));
         if (ns >= 0 && ns < 3) h_Kns[v][Qv[v]][ns]->Fill(fabs(Kv[v]));
         h_Kall[v][Qv[v]]->Fill(fabs(Kv[v]));
      }

      if (qcorr != q) { ++nFlip; h_flip->Fill(ptS[0]); }

      // per stage: orig and an_2024 use the raw charge, corr the corrected one
      int Qs[NSTAGE] = {q, q, qcorr};
      for (int s = 0; s < NSTAGE; ++s) {
         h_K[s]->Fill(Ks[s]);
         h_pt[s]->Fill(ptS[s]);
         h_ptq[s][Qs[s]]->Fill(ptS[s]);
         for (int k = 0; k < NTHR; ++k)
            if (ptS[s] > PT_THR[k]) h_npass[s]->Fill(k);
      }

      // second good track
      bool has_2goodtracks = false;
      if (hasTrk2) {
         int    idx2   = int(b_idx2);
         double pt2    = double(b_pt2);
         double dxy2   = double(b_dxy2);
         int    qual2  = int(b_qual2);
         int    nstub2 = int(b_nstub2);
         if (idx2 < 90 && pt2 > 100 && dxy2 < 1 && qual2 > 12) {
            if (nstub2 == 2 && qual2 > 12) has_2goodtracks = true;
            if (nstub2 == 3 && qual2 > 13) has_2goodtracks = true;
            if (nstub2 == 4 && qual2 > 14) has_2goodtracks = true;
         }
      }

      // model-dependent categories
      bool is_slow = (stationspread == 4320 || stationspread == 3210);
      int  cat     = (bxspread > 0) ? category(nstub, bxspread, is_slow) : -1;
      if (cat >= 0 && (!NEEDS_2TRK[cat] || has_2goodtracks))
         for (int m = 0; m < NMOD; ++m) {
            double pt = ptS[MOD_STAGE[m]];
            if (pt > catPtMin(cat)) h_cat[m][cat]->Fill(pt);
         }
   } // end of event loop

   auto frac = [&](long n){ return nTrk ? double(n)/nTrk : 0.; };
   auto ferr = [&](long n){ double p = frac(n); return nTrk ? std::sqrt(p*(1-p)/nTrk) : 0.; };
   printf("\n  selected tracks: %ld\n", nTrk);
   printf("  sign(K) == charge for %.4f +/- %.4f of tracks\n", frac(nSignOK), ferr(nSignOK));
   printf("  charge sign changed by the correction: %.4f +/- %.4f of tracks\n", frac(nFlip), ferr(nFlip));

   TFile *fout = TFile::Open(output.c_str(), "RECREATE");
   TDirectory* top = fout->mkdir(name.c_str());

   top->mkdir("common")->cd();
   h_flip->Write();
   h_charge->Write();

   for (int v = 0; v < NKV; ++v) {
      top->mkdir(KV_NAME[v])->cd();
      for (int q = 0; q < 2; ++q) {
         for (int i = 0; i < kcorr::NSEC; ++i)
            for (int j = 0; j < kcorr::NETA; ++j) h_Kmap[v][q][i][j]->Write();
         for (int n = 0; n < 3; ++n) h_Kns[v][q][n]->Write();
         h_Kall[v][q]->Write();
      }
   }

   for (int s = 0; s < NSTAGE; ++s) {
      top->mkdir(STAGE_NAME[s])->cd();
      h_K[s]->Write();
      h_pt[s]->Write();
      for (int q = 0; q < 2; ++q) h_ptq[s][q]->Write();
      h_npass[s]->Write();
   }

   for (int m = 0; m < NMOD; ++m) {
      top->mkdir(MOD_DIR[m])->cd();
      for (int c = 0; c < NCAT; ++c) h_cat[m][c]->Write();
   }

   fout->Close();
   return 0;
}