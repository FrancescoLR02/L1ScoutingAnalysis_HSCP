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


// g++ -O3 HistNew_pT.cc -o HistNew_pT.exe $(root-config --cflags --glibs)

using namespace std;

// binning of the charge-odd curvature map
const int    NPHI      = 12;    // 30 deg bins
const int    NETA      = 5;     // over |eta| < 0.83
const int    NPHI_FINE = 72;    // 5 deg bins, for the charge count asymmetry
const double ETAMAX    = 0.83;
const double D = 0;



double ptLUT(double K,int& charge) {
  charge = (K >= 0) ? +1 : -1;
  float lsb = 1.25 / float(1 << 13);
  float FK = fabs(K);

  if (FK > 2047)
    FK = 2047.;
  if (FK < 9)
    FK = 9.;

  FK = FK * lsb;
  //step 1 -material and B-field
  FK = .8569 * FK / (1.0 + 0.1144 * FK);
  //step 2 - misalignment
  FK = FK - charge * 1.23e-03;

   if (FK < 0){
   charge = -charge;
   FK = -FK;
  }

  //Get to BMTF scale
//   FK = FK / 1.17;



  int pt = 0;
  if (FK != 0) pt = int(2.0 / FK);
   if (pt < 8) pt = 8;

  return pt/2;
}




int main(int argc, char** argv) {

   std::string input = *(argv + 1);
   std::string output = *(argv + 2);

   TFile *f_Double = new TFile(input.c_str());
   cout<<"XXXXXXXXXXXXX "<<input.c_str()<<" XXXXXXXXXXXX"<<endl;
   TTree *arbre1 = (TTree*) f_Double->Get("Events");



   arbre1->SetBranchAddress("bxspread1", &bxspread1);
   arbre1->SetBranchAddress("bxspread2", &bxspread2);
   arbre1->SetBranchAddress("nstub1", &nstub1);
   arbre1->SetBranchAddress("nstub2", &nstub2);
   arbre1->SetBranchAddress("pt1", &pt1);
   arbre1->SetBranchAddress("eta1", &eta1);
   arbre1->SetBranchAddress("phi1", &phi1);
   arbre1->SetBranchAddress("pt2", &pt2);
   arbre1->SetBranchAddress("eta2", &eta2);
   arbre1->SetBranchAddress("phi2", &phi2);

   arbre1->SetBranchAddress("HwK1", &hwK1);
   arbre1->SetBranchAddress("HwK2", &hwK2);
   arbre1->SetBranchAddress("recobeta1", &recobeta1);
   arbre1->SetBranchAddress("new_pT1", &new_pT1);



   // charge/qual/dxy change scalar type from one ntuple version to the next
   // (dxy: Int_t in the 2024 skims, Double_t in simulation, Float_t in data),
   // so they are connected through ScalarBranch, which adapts to the type on file.
   charge1.connect(arbre1, "charge1");
   qual1.connect(arbre1, "qual1");
   dxy1.connect(arbre1, "dxy1");
   charge2.connect(arbre1, "charge2");
   qual2.connect(arbre1, "qual2");
   dxy2.connect(arbre1, "dxy2");
   

   const int    NK    = 1200;
   const double KMIN  = -500.0, KMAX = 500.0;  
   const double lsb = 1.25 / float(1 << 13);


   TH1F* h_K  = new TH1F("h_K", "hw curvature K", 1000, KMIN, KMAX); h_K->Sumw2();
   TH1F* h_pt  = new TH1F("h_pt", "h_pt", 200, 0, 1500); h_pt->Sumw2();
   TH1F* h_pt_pos  = new TH1F("h_pt_pos", "h_pt_pos", 200, 12, 1500); h_pt_pos->Sumw2();
   TH1F* h_pt_neg  = new TH1F("h_pt_neg", "h_pt_neg", 200, 12, 1500); h_pt_neg->Sumw2();

   TH1F* h_newpt  = new TH1F("h_newpt", "h_newpt", 200, 0, 1500); h_newpt->Sumw2();
   TH1F* h_newpt_pos  = new TH1F("h_newpt_pos", "h_newpt_pos", 200, 12, 1500); h_newpt_pos->Sumw2();
   TH1F* h_newpt_neg  = new TH1F("h_newpt_neg", "h_newpt_neg", 200, 12, 1500); h_newpt_neg->Sumw2();


   TH1F* h_charge = new TH1F("h_charge","h_charge", 3, -1.5, 1.5); h_charge->Sumw2();

   TH1F* h_beta = new TH1F("h_beta", "h_beta", 50, 0, 1); h_beta->Sumw2();
   TH1F* misID_pt = new TH1F("misID_pt", "misID_pt", 100, 12.5, 1000); misID_pt->Sumw2();




   static TRandom3 randGen(1234);
   float ptmin=0.0;
   float ptmax=10000.0;


   int nStub2_misID = 0; int nStub3_misID = 0; int nStub4_misID = 0;

   Int_t nentries_wtn = (Int_t) arbre1->GetEntries();
   for (Int_t i = 0; i < nentries_wtn; i++) {
   	arbre1->GetEntry(i);
      if (i % 100000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
      fflush(stdout);


      if (dxy1>=1 or dxy2>=1) continue;
      if (qual1<12 or qual2<12) continue;
      if (nstub1==4 and qual1<14) continue;
      if (nstub1==3 and qual1<13) continue;
      if (nstub1==2) continue;
      if (nstub2==4 and qual2<14) continue;
      if (nstub2==3 and qual2<13) continue;
      if (nstub2==2) continue;

      int q = 0;
      double ptLUT1 = ptLUT(hwK1, q);


      h_K->Fill(hwK1);
      h_pt->Fill(pt1);
      h_newpt->Fill(ptLUT1);

      if (charge1 > 0){ 
         h_pt_pos->Fill(pt1);
         h_newpt_pos->Fill(ptLUT1);
      }
      else{ 
         h_pt_neg->Fill(pt1);
         h_newpt_neg->Fill(ptLUT1);
      }
      
      h_charge->Fill(q);// h_charge->Fill(charge2);
      


      if(recobeta1 < 0.95){
         h_beta->Fill(beta1);
         misID_pt->Fill(pt1); //misID_pt->Fill(pt2);


      }
   } // end of loop over events

   

   TFile *fout = TFile::Open(output.c_str(), "RECREATE");
   fout->cd();

   h_K->Write();

   h_charge->Write();
   h_pt->Write();
   h_pt_pos->Write();
   h_pt_neg->Write();

   h_newpt->Write();
   h_newpt_pos->Write();
   h_newpt_neg->Write();

   TDirectory* dir2=fout->mkdir("misID_BX");
   dir2->cd();
   h_beta->Write();
   misID_pt->Write();

   fout->Close();
}

