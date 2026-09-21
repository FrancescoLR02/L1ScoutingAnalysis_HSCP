
#include <stdio.h>
#include <TF1.h>
#include <TDirectoryFile.h>
#include <TRandom3.h>
#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TH1D.h"
#include "TKey.h"
#include "THashList.h"
#include "THStack.h"
#include "TPaveLabel.h"
#include "TFile.h"
#include "TTree.h"
//#include "zmumu_Tree.h"
#include "Correction_tr_Tree.h"
#include "KCorrection.h"
#include <TRandom3.h>
#include <algorithm>

//! g++ -O3 Validate_dXY.cc -o Validate_dXY.exe $(root-config --cflags --glibs)

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

   arbre1->SetBranchAddress("bxspread1", &bxspread1);
   arbre1->SetBranchAddress("bxspread2", &bxspread2);
   arbre1->SetBranchAddress("pt1", &pt1);
   arbre1->SetBranchAddress("charge1", &charge1);
   arbre1->SetBranchAddress("qual1", &qual1);
   arbre1->SetBranchAddress("pt2", &pt2);

   arbre1->SetBranchAddress("charge2", &charge2);
   arbre1->SetBranchAddress("qual2", &qual2);


   nstub1.connect(arbre1, "nstub1");  nstub2.connect(arbre1, "nstub2");
   eta1.connect(arbre1, "eta1");      eta2.connect(arbre1, "eta2");
   phi1.connect(arbre1, "phi1");      phi2.connect(arbre1, "phi2");
   hwK1.connect(arbre1, "hwK1");      hwK2.connect(arbre1, "hwK2");
   dxy1.connect(arbre1, "dxy1");      dxy2.connect(arbre1, "dxy2");

   TH1D* h_mmumu_OS=new TH1D("h_mmumu_OS", "h_mmumu_OS", 50,50,150); h_mmumu_OS->Sumw2();
   TH1D* h_mmumu_OS_corr=new TH1D("h_mmumu_OS_corr", "h_mmumu_OS_corr", 50,50,150); h_mmumu_OS_corr->Sumw2();
   TH1D* h_mmumu_SS=new TH1D("h_mmumu_SS", "h_mmumu_SS", 50,50,150); h_mmumu_SS->Sumw2();


   TH1D* h_dxy = new TH1D("h_dxy","h_dxy",80,0,4); h_dxy->Sumw2();
   TH1D* h_dxy_3stub = new TH1D("h_dxy_3stub","h_dxy_3stub",80,0,4); h_dxy_3stub->Sumw2();
   TH1D* h_dxy_4stub = new TH1D("h_dxy_4stub","h_dxy_4stub",80,0,4); h_dxy_4stub->Sumw2();
   TH1D* h_dxy_SS = new TH1D("h_dxy_SS","h_dxy_SS",80,0,4);



   static TRandom3 randGen(1234);
   float ptmin=0.0;
   float ptmax=10000.0;

   Int_t nentries_wtn = (Int_t) arbre1->GetEntries();
   for (Int_t i = 0; i < nentries_wtn; i++) {
   	arbre1->GetEntry(i);
      if (i % 10000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
      fflush(stdout);

      float aweight=1;//xsweight*lumiweight;

      if (dxy1>=4 or dxy2>=4) continue;
      if (pt1<15 or pt2<15) continue;
      //if (pt1<100 or pt2<100) continue;
      if (qual1<12 or qual2<12) continue;

      if (nstub1==4 and qual1<14) continue;
      if (nstub1==3 and qual1<13) continue;
      if (nstub1==2) continue;
      if (nstub2==4 and qual2<14) continue;
      if (nstub2==3 and qual2<13) continue;
      if (nstub2==2) continue;



      TLorentzVector my_mu1; my_mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.105);
      TLorentzVector my_mu2; my_mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.105);

      double Kcorr1 = kcorr::correctK(hwK1, phi1, eta1, nstub1);
      double Kcorr2 = kcorr::correctK(hwK2, phi2, eta2, nstub2);

      double pt1_corr = kcorr::ptLUT_corr(Kcorr1, true);
      double pt2_corr = kcorr::ptLUT_corr(Kcorr2, true);


      TLorentzVector my_mu1_corr; my_mu1_corr.SetPtEtaPhiM(pt1_corr, eta1, phi1, 0.105);
      TLorentzVector my_mu2_corr; my_mu2_corr.SetPtEtaPhiM(pt2_corr, eta2, phi2, 0.105);


      if (name=="data_obs") aweight=1.0;

      float new_mmumu=(my_mu1+my_mu2).M();
      float mmumu_corr = (my_mu1_corr + my_mu2_corr).M();
      if (charge1*charge2<0){
         h_mmumu_OS->Fill(new_mmumu, aweight);
         h_mmumu_OS_corr->Fill(mmumu_corr, aweight);
         if (new_mmumu>80 and new_mmumu < 100 ){

            h_dxy->Fill(dxy1); h_dxy->Fill(dxy2);

            if(nstub1 == 3) h_dxy_3stub->Fill(dxy1); 
            if(nstub2 == 3) h_dxy_3stub->Fill(dxy2); 
            if(nstub1 == 4) h_dxy_4stub->Fill(dxy1); 
            if(nstub2 == 4) h_dxy_4stub->Fill(dxy2); 

         }
      }
      if(charge1*charge2>0){
         h_dxy_SS->Fill(dxy1); h_dxy_SS->Fill(dxy2);
      }


   } // end of loop over events

   TFile *fout = TFile::Open(output.c_str(), "RECREATE");
   fout->cd();


   TDirectory* dir1=fout->mkdir("OS");
   dir1->cd();
   h_mmumu_OS->SetName(name.c_str());
   if (name=="DY") h_mmumu_OS->SetName("DY");
   h_mmumu_OS->Write();
   h_mmumu_OS_corr->Write();

   h_dxy->Write();
   h_dxy_3stub->Write();
   h_dxy_4stub->Write();
   h_dxy_SS->Write();


   fout->Close();
}

