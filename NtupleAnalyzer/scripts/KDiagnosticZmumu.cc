// Diagnostic of the charge-odd curvature bias on Z->mumu, AFTER the misalignment correction
// applied inside the KBMTF emulator. hwK on file is already corrected: nothing is applied here,
// the maps only measure what is left (the notebook KDiagnosticZmumu.ipynb reads the output).

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
#include "modzmumu_Tree.h"

//! g++ -O3 KDiagnosticZmumu.cc -o KDiagnosticZmumu.exe $(root-config --cflags --glibs)
//! ./KDiagnosticZmumu.exe input.root output.root <name>     (name = DYSIM fills the gen histograms)

using namespace std;

// same binning as the derivation of the correction and as the emulator
const int    NPHI   = 12;    // 30 deg bins, bin k centred at k*30 deg (= sector k+1)
const int    NETA   = 5;     // over |eta| < 0.83
const double ETAMAX = 0.83;

int phiBin(double phi){
   const double w = 2*TMath::Pi()/NPHI;
   double x = std::fmod(phi + 0.5*w, 2*TMath::Pi());
   if (x < 0) x += 2*TMath::Pi();
   return std::min(std::max(int(x / w), 0), NPHI-1);
}

int etaBin(double eta){
   int b = int((eta + ETAMAX) / (2*ETAMAX) * NETA);
   return std::min(std::max(b, 0), NETA-1);
}

// muon selection
bool passSel(double dxy, int qual, int nstub, double pt){
   if (dxy >= 1)                return false;
   if (pt < 15)                  return false;
   if (nstub == 4 && qual < 14) return false;
   if (nstub == 3 && qual < 13) return false;
   //if (nstub == 2 && qual < 12) return false;
   if (nstub == 2) return false;
   return true;
}


int main(int argc, char** argv) {

   if (argc < 4) { cerr << "usage: " << argv[0] << " input.root output.root <name>" << endl; return 1; }
   std::string input  = argv[1];
   std::string output = argv[2];
   std::string name   = argv[3];
   const bool isSim   = (name == "DYSIM");

   TFile *fin = TFile::Open(input.c_str(), "READ");
   TTree *arbre1 = (TTree*) fin->Get("Events");
   cout << "XXXXXXXXXXXXX " << input << " XXXXXXXXXXXX" << endl;
   cout << "branches on file: " << arbre1->GetListOfBranches()->GetEntries()
        << "   compression: " << fin->GetCompressionAlgorithm() << "  (1=zlib 2=LZMA 4=LZ4 5=ZSTD)" << endl;

   // read ONLY the branches used below: GetEntry otherwise decompresses the whole skim
   std::vector<std::string> used = {"nstub1", "nstub2", "pt1", "eta1", "phi1", "pt2", "eta2", "phi2",
                                    "charge1", "charge2", "qual1", "qual2", "dxy1", "dxy2",
                                    "hwK1", "hwK2", "mmumu"};
   
   if (isSim) { used.push_back("genpt1"); used.push_back("genpt2"); }

   arbre1->SetBranchStatus("*", 0);
   arbre1->SetCacheSize(100 * 1024 * 1024);   // 100 MB read-ahead, matters a lot on EOS
   for (const auto& b : used) {
      arbre1->SetBranchStatus(b.c_str(), 1);
      arbre1->AddBranchToCache(b.c_str(), true);
   }

   nstub1.connect(arbre1, "nstub1");
   nstub2.connect(arbre1, "nstub2");
   pt1.connect(arbre1, "pt1");
   eta1.connect(arbre1, "eta1");
   phi1.connect(arbre1, "phi1");
   pt2.connect(arbre1, "pt2");
   eta2.connect(arbre1, "eta2");
   phi2.connect(arbre1, "phi2");
   mmumu.connect(arbre1, "mmumu");
   if (isSim) {
      genpt1.connect(arbre1, "genpt1");
      genpt2.connect(arbre1, "genpt2");
   }

   // types change between ntuple versions: ScalarBranch adapts to the type on file
   charge1.connect(arbre1, "charge1");
   qual1.connect(arbre1, "qual1");
   dxy1.connect(arbre1, "dxy1");
   hwK1.connect(arbre1, "hwK1");
   charge2.connect(arbre1, "charge2");
   qual2.connect(arbre1, "qual2");
   dxy2.connect(arbre1, "dxy2");
   hwK2.connect(arbre1, "hwK2");

   TH1::SetDefaultSumw2(kTRUE);

   // |K| of mu+ and mu- in each (phi, eta) cell
   TH1F* h_plus[NPHI][NETA];
   TH1F* h_minus[NPHI][NETA];
   for (int i = 0; i < NPHI; ++i){
      for (int j = 0; j < NETA; ++j){
         h_plus[i][j]  = new TH1F(Form("h_K_plus_phi%d_eta%d" ,i,j), "|K|, positive muons;|hwK| [LSB]", 200, 0, 450);
         h_minus[i][j] = new TH1F(Form("h_K_minus_phi%d_eta%d",i,j), "|K|, negative muons;|hwK| [LSB]", 200, 0, 450);
      }
   }

   // |K| of mu+ and mu- vs number of stubs (2, 3, 4)
   TH1F* h_plus_nStub[3];
   TH1F* h_minus_nStub[3];
   for (int n = 0; n < 3; ++n){
      h_plus_nStub[n]  = new TH1F(Form("h_K_plus_nStub%d" ,n+2), "|K|, positive muons;|hwK| [LSB]", 200, 0, 450);
      h_minus_nStub[n] = new TH1F(Form("h_K_minus_nStub%d",n+2), "|K|, negative muons;|hwK| [LSB]", 200, 0, 450);
   }

   TH1F* h_mmumu  = new TH1F("h_mmumu",  ";m_{#mu#mu} [GeV];events",                 100, 50, 150);
   TH1F* h_pt     = new TH1F("h_pt",     ";L1 p_{T} [GeV];muons",                     100, 0, 1000);
   TH1F* h_genpt  = new TH1F("h_genpt",  ";gen p_{T} [GeV];muons",                    100, 0, 1000);
   TH1F* h_ptResp = new TH1F("h_ptResp", ";p_{T}^{L1}/p_{T}^{gen} - 1;muons",         100, -1, 1);
   TH1F* h_hwK     = new TH1F("h_hwK",     "h_hwK", 1000, -600, 600);

   // one muon: |K| into its (phi, eta) cell and its nStub histogram
   auto fillMuon = [&](double K, double charge, double phi, double eta, int nstub){
      int i = phiBin(phi), j = etaBin(eta), n = nstub - 2;
      if (charge > 0) h_plus[i][j]->Fill(fabs(K));
      else            h_minus[i][j]->Fill(fabs(K));
      if (n < 0 || n > 2) return;
      if (charge > 0) h_plus_nStub[n]->Fill(fabs(K));
      else            h_minus_nStub[n]->Fill(fabs(K));
   };

   long nSel = 0;
   Long64_t nentries = arbre1->GetEntries();
   for (Long64_t ev = 0; ev < nentries; ++ev) {
      arbre1->GetEntry(ev);
      if (ev % 2000000 == 0) { fprintf(stdout, "\r  Processed events: %8lld of %8lld ", ev, nentries); fflush(stdout); }

      if (!passSel(double(dxy1), int(qual1), int(nstub1), pt1)) continue;
      if(name == "DY" || name == "SIMDY") if (!passSel(double(dxy2), int(qual2), int(nstub2), pt2)) continue;
      ++nSel;

      // hwK is already corrected by the emulator: use it as it is
      fillMuon(double(hwK1), double(charge1), phi1, eta1, int(nstub1));
      if(name == "DY" || name == "SIMDY") fillMuon(double(hwK2), double(charge2), phi2, eta2, int(nstub2));

      if(name=="DYSIM"){ h_pt->Fill(pt1);}// h_pt->Fill(pt2);}
      else {h_pt->Fill(pt1);}
      h_mmumu->Fill(mmumu);
      h_hwK->Fill(hwK1); //h_hwK->Fill(hwK2); 

      if (isSim) {
         h_genpt->Fill(genpt1); //h_genpt->Fill(genpt2);
         if (genpt1 > 0) h_ptResp->Fill(pt1 / genpt1 - 1);
         //if (genpt2 > 0) h_ptResp->Fill(pt2 / genpt2 - 1);
      }
   }
   printf("\n  selected events: %ld\n", nSel);

   TFile *fout = TFile::Open(output.c_str(), "RECREATE");
   fout->cd();
   h_mmumu->Write();
   h_pt->Write();
   h_genpt->Write();
   h_ptResp->Write();
   h_hwK->Write();

   fout->mkdir("KmapPhiEta")->cd();
   for (int i = 0; i < NPHI; ++i)
      for (int j = 0; j < NETA; ++j){ h_plus[i][j]->Write(); h_minus[i][j]->Write(); }

   fout->mkdir("KmapnStub")->cd();
   for (int n = 0; n < 3; ++n){ h_plus_nStub[n]->Write(); h_minus_nStub[n]->Write(); }

   fout->Close();
   fin->Close();
   return 0;
}