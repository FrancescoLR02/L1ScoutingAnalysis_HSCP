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
#include "modzmumu_Tree.h"


// g++ -O3 CorrectionFactorsZmumu.cc -o CorrectionFactorsZmumu.exe $(root-config --cflags --glibs)

using namespace std;

int main(int argc, char** argv) {

   std::string input = *(argv + 1);
   std::string output = *(argv + 2);
   std::string name = *(argv + 3);

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
//  arbre1->SetBranchAddress("isL1MuMatched1", &isL1MuMatched1);
//  arbre1->SetBranchAddress("isL1MuMatched2", &isL1MuMatched2);
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

   TH1F* h_pt_OS=new TH1F("h_pt_OS", "h_pt_OS", 13,0,260); h_pt_OS->Sumw2();
   TH1F* h_pt_SS=new TH1F("h_pt_SS", "h_pt_SS", 13,0,260); h_pt_SS->Sumw2();


   TH1F* h_K  = new TH1F("h_K", "hw curvature K", NK, KMIN, KMAX); h_K->Sumw2();
   TH1F* h_K1 = new TH1F("h_K1","hw curvature K, mu1", NK, KMIN, KMAX); h_K1->Sumw2();
   TH1F* h_K_plus  = new TH1F("h_K_plus", "hw curvature K, K>0", NK/2, 0, 150); h_K_plus->Sumw2();
   TH1F* h_K_minus = new TH1F("h_K_minus","hw curvature K, K<0", NK/2, 0, 150); h_K_minus->Sumw2();
   TH1F* h_dxy  = new TH1F("h_dxy", "h_dxy", 100, 0, 1); h_dxy->Sumw2();
   TH1F* h_nstub  = new TH1F("h_nstub", "h_nstub", 3, 2, 5); h_nstub->Sumw2();




   TH1F* h_charge = new TH1F("h_charge","h_charge", 3, -1.5, 1.5); h_charge->Sumw2();
   TH1F* h_pt = new TH1F("h_pt","h_pt",400, 12.5, 1100); h_pt->Sumw2();

   TH1F* h_beta = new TH1F("h_beta", "h_beta", 50, 0, 1); h_beta->Sumw2();
   TH1F* misID_pt = new TH1F("misID_pt", "misID_pt", 400, 12.5, 1100); misID_pt->Sumw2();
   TH1F* misID_dxy = new TH1F("misID_dxy", "misID_dxy", 100, 0, 1); misID_dxy->Sumw2();
   TH1F* misID_nstub = new TH1F("misID_nstub", "misID_nstub", 3, 2, 5); misID_nstub->Sumw2();
   TH1F* misID_K = new TH1F("misID_K", "misID_K", 60, -0.03, 0.03); misID_K->Sumw2();
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
      if (i % 10000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
      fflush(stdout);

      

      TLorentzVector my_mu1; my_mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.105);
      TLorentzVector my_mu2; my_mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.105);


      if (dxy1>=1 or dxy2>=1) continue;
      if (pt1<15 or pt2<15) continue;
      //if (pt1<100 or pt2<100) continue;
      if (qual1<12 or qual2<12) continue;

      if (nstub1==4 and qual1<14) continue;
      if (nstub1==3 and qual1<13) continue;
      if (nstub1==2 and qual1<12) continue;
      if (nstub2==4 and qual2<14) continue;
      if (nstub2==3 and qual2<13) continue;
      if (nstub2==2 and qual2<12) continue;

      h_K->Fill(hwK1);// h_K->Fill(hwK2);
      //h_K1->Fill(hwK1); //h_K2->Fill(hwK2);
      if (hwK1 > 0) h_K_plus->Fill(hwK1);
      //if (hwK2 > 0) h_K_plus->Fill(hwK2);
      if (hwK1 < 0) h_K_minus->Fill(-hwK1);
      //if (hwK2 < 0) h_K_minus->Fill(hwK2);
      
      h_charge->Fill(charge1); //h_charge->Fill(charge2);
      h_pt->Fill(pt1); //h_pt->Fill(pt2);
      h_dxy->Fill(dxy1);
      h_nstub->Fill(nstub1);


      if(beta1 < 0.95){
         h_beta->Fill(beta1);
         misID_pt->Fill(pt1);
         misID_dxy->Fill(dxy1);
         misID_nstub->Fill(nstub1);
         misID_K->Fill(hwK1*lsb);
         misID_invpT->Fill(charge1/pt1);
         
         if(charge1*charge2<0) misID_mmumuOS->Fill((my_mu1+my_mu2).M());
         if(charge1*charge2>0) misID_mmumuSS->Fill((my_mu1+my_mu2).M());

         if(nstub1 == 2) nStub2_misID += 1;
         if(nstub1 == 3) nStub3_misID += 1;
         if(nstub1 == 4) nStub4_misID += 1;
      }
   } // end of loop over events

   

   TFile *fout = TFile::Open(output.c_str(), "RECREATE");
   fout->cd();

   TDirectory* dir1=fout->mkdir("Zmumu");
   dir1->cd();
   //h_mmumu_OS->SetName(name.c_str());
   h_mmumu_OS->Write();
   h_mmumu_SS->Write();
   h_K->Write();
   //h_K2->Write();
   h_K_plus->Write();
   h_K_minus->Write();

   h_charge->Write();
   h_pt->Write();
   h_dxy->Write();
   h_nstub->Write();

   TDirectory* dir2=fout->mkdir("misID_BX");
   dir2->cd();
   //h_beta->Write();
   misID_pt->Write();
   misID_dxy->Write();
   misID_nstub->Write();
   misID_K->Write();
   misID_invpT->Write();
   misID_mmumuOS->Write();
   misID_mmumuSS->Write();

   fout->Close();
}

