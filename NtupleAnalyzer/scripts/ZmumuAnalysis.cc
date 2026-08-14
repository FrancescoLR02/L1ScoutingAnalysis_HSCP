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

using namespace std;

int main(int argc, char** argv) {

   std::string input = *(argv + 1);
   std::string output = *(argv + 2);
   std::string name = *(argv + 3);
   std::string muobject = *(argv + 4);

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
   arbre1->SetBranchAddress("genpt1", &genpt1);
   arbre1->SetBranchAddress("genpt2", &genpt2);
   arbre1->SetBranchAddress("stub1Bx1", &stub1Bx1);
   arbre1->SetBranchAddress("stub2Bx1", &stub2Bx1);
   arbre1->SetBranchAddress("stub3Bx1", &stub3Bx1);
   arbre1->SetBranchAddress("stub4Bx1", &stub4Bx1);
   arbre1->SetBranchAddress("stub1Bx2", &stub1Bx2);
   arbre1->SetBranchAddress("stub2Bx2", &stub2Bx2);
   arbre1->SetBranchAddress("stub3Bx2", &stub3Bx2);
   arbre1->SetBranchAddress("stub4Bx2", &stub4Bx2);
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

   TH1F* h_mmumu_OS=new TH1F("h_mmumu_OS", "h_mmumu_OS", 50,50,150); h_mmumu_OS->Sumw2();
   TH1F* h_mmumu_SS=new TH1F("h_mmumu_SS", "h_mmumu_SS", 50,50,150); h_mmumu_SS->Sumw2();

   TH1F* h_nstub_OS=new TH1F("h_nstub_OS", "h_nstub_OS", 3,2,5); h_nstub_OS->Sumw2();
   TH1F* h_nstub_SS=new TH1F("h_nstub_SS", "h_nstub_SS", 3,2,5); h_nstub_SS->Sumw2();

   TH1F* h_pt_OS=new TH1F("h_pt_OS", "h_pt_OS", 13,0,260); h_pt_OS->Sumw2();
   TH1F* h_pt_SS=new TH1F("h_pt_SS", "h_pt_SS", 13,0,260); h_pt_SS->Sumw2();

   TH1F* h_met_OS=new TH1F("h_met_OS", "h_met_OS", 20,0,100); h_met_OS->Sumw2();
   TH1F* h_met_SS=new TH1F("h_met_SS", "h_met_SS", 20,0,100); h_met_SS->Sumw2();

   TH1F* h_bxspan = new TH1F("h_bxspan","h_bxspan",6,0,6); h_bxspan->Sumw2();
   TH1F* h_bx1 = new TH1F("h_bx1","h_bx1",11,-5,6); h_bx1->Sumw2();
   TH1F* h_bx2 = new TH1F("h_bx2","h_bx2",11,-5,6); h_bx2->Sumw2();
   TH1F* h_bx3 = new TH1F("h_bx3","h_bx3",11,-5,6); h_bx3->Sumw2();
   TH1F* h_bx4 = new TH1F("h_bx4","h_bx4",11,-5,6); h_bx4->Sumw2();

   TH1F* h_reso = new TH1F("h_reso","h_reso",20,-1,1); h_reso->Sumw2();
   TH1F* h_reso_pt15to25 = new TH1F("h_reso_pt15to25","h_reso_pt15to25",20,-1,1); h_reso_pt15to25->Sumw2();
   TH1F* h_reso_pt25to35 = new TH1F("h_reso_pt25to35","h_reso_pt25to35",20,-1,1); h_reso_pt25to35->Sumw2();
   TH1F* h_reso_pt35to45 = new TH1F("h_reso_pt35to45","h_reso_pt35to45",20,-1,1); h_reso_pt35to45->Sumw2();
   TH1F* h_reso_pt45to55 = new TH1F("h_reso_pt45to55","h_reso_pt45to55",20,-1,1); h_reso_pt45to55->Sumw2();
   TH1F* h_reso_ptgt55 = new TH1F("h_reso_ptgt55","h_reso_ptgt55",20,-1,1); h_reso_ptgt55->Sumw2();

   TH1F* h_station = new TH1F("h_station","h_station",5,0,5); h_station->Sumw2();


   //float lumiweight= (40400.0/15.046) + (5820.0/15.060) + (12320.0/15.048);
   float lumiweight = 68.165642 + 35.022010 + 458.565962;
   //float lumiweight = 3918;

   static TRandom3 randGen(1234);
   float ptmin=0.0;
   float ptmax=10000.0;

   Int_t nentries_wtn = (Int_t) arbre1->GetEntries();
   for (Int_t i = 0; i < nentries_wtn; i++) {
   	arbre1->GetEntry(i);
      if (i % 10000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
      fflush(stdout);

      float aweight=6346 * 0.3664 * lumiweight/1478803;//xsweight*lumiweight;

      
      // // Muon energy calibration
      // if (name=="data_obs"){
      //    pt1=pt1;
      //    pt2=pt2;
      // }
      // else{
      //    pt1=0.99*pt1;
      //    pt2=0.99*pt2;
      //    }

      // Muon energy resolution
      double smearFactor1 = randGen.Gaus(1.0, 0.08);
      double smearFactor2 = randGen.Gaus(1.0, 0.08);
      if (name!="data_obs"){
         pt1 = pt1 * smearFactor1;
         pt2 = pt2 * smearFactor2;
      }

      // Muon reconstruction efficiency
      //float musf=1.00;
      //aweight = aweight*musf*musf;
      // if (nstub1==2) aweight = aweight*1.05*1.025;
      // else if (nstub1==3) aweight = aweight*1.01*1.025;
      // else if (nstub1==4) aweight = aweight*0.86*1.025;
      // if (nstub2==2) aweight = aweight*1.05*1.025;
      // else if (nstub2==3) aweight = aweight*1.01*1.025;
      // else if (nstub2==4) aweight = aweight*0.86*1.025;

      if (nstub1==2) aweight = aweight*1.14*1.025;
      else if (nstub1==3) aweight = aweight*1.04*1.025;
      else if (nstub1==4) aweight = aweight*0.78*1.025;
      if (nstub2==2) aweight = aweight*1.14*1.025;
      else if (nstub2==3) aweight = aweight*1.04*1.025;
      else if (nstub2==4) aweight = aweight*0.78*1.025;
      

      TLorentzVector my_mu1; my_mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.105);
      TLorentzVector my_mu2; my_mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.105);

      // Prompt muons. In the 2024 skims dxy was the 2-bit hardware index (0 = prompt),
      // in the modified KBMTF ntuples it is a continuous displacement peaking below 1;
      // dxy<1 is the prompt selection in both conventions.
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

      //if (nstub1!=2 or nstub2!=2) continue;

      if (name=="data_obs") aweight=1.0;

      if (name!="data_obs"){ 
         //cout<<(pt1-genpt1)/genpt1<<" "<<pt1<<" "<<genpt1<<endl;
         h_reso->Fill((pt1-genpt1)/genpt1);
         h_reso->Fill((pt2-genpt2)/genpt2);
         if (pt1>=15 and pt1<25) h_reso_pt15to25->Fill((pt1-genpt1)/genpt1);
         else if (pt1>=25 and pt1<35) h_reso_pt25to35->Fill((pt1-genpt1)/genpt1);
         else if (pt1>=35 and pt1<45) h_reso_pt35to45->Fill((pt1-genpt1)/genpt1);
         else if (pt1>=45 and pt1<55) h_reso_pt45to55->Fill((pt1-genpt1)/genpt1);
         else if (pt1>=55 and pt1<65) h_reso_ptgt55->Fill((pt1-genpt1)/genpt1);

         if (pt2>=15 and pt2<25) h_reso_pt15to25->Fill((pt2-genpt2)/genpt2);
         else if (pt2>=25 and pt2<35) h_reso_pt25to35->Fill((pt2-genpt2)/genpt2);
         else if (pt2>=35 and pt2<45) h_reso_pt35to45->Fill((pt2-genpt2)/genpt2);
         else if (pt2>=45 and pt2<55) h_reso_pt45to55->Fill((pt2-genpt2)/genpt2);
         else if (pt2>=55 and pt2<65) h_reso_ptgt55->Fill((pt2-genpt2)/genpt2);
      }

      float new_mmumu=(my_mu1+my_mu2).M();
      if (charge1*charge2<0){
         h_mmumu_OS->Fill(new_mmumu, aweight);
         if (new_mmumu>70){
            h_bxspan->Fill(std::max({stub1Bx1, stub2Bx1, stub3Bx1, stub4Bx1})-std::min({stub1Bx1, stub2Bx1, stub3Bx1, stub4Bx1}),aweight);
            h_bxspan->Fill(std::max({stub1Bx2, stub2Bx2, stub3Bx2, stub4Bx2})-std::min({stub1Bx2, stub2Bx2, stub3Bx2, stub4Bx2}),aweight);
            if (stub1Station1==1) h_bx1->Fill(stub1Bx1,aweight);
            else if (stub2Station1==1) h_bx1->Fill(stub2Bx1,aweight);
            else if (stub3Station1==1) h_bx1->Fill(stub3Bx1,aweight);
            else if (stub4Station1==1) h_bx1->Fill(stub4Bx1,aweight);
            if (stub1Station1==2) h_bx2->Fill(stub1Bx1,aweight);
            else if (stub2Station1==2) h_bx2->Fill(stub2Bx1,aweight);
            else if (stub3Station1==2) h_bx2->Fill(stub3Bx1,aweight);
            else if (stub4Station1==2) h_bx2->Fill(stub4Bx1,aweight);
            if (stub1Station1==3) h_bx3->Fill(stub1Bx1,aweight);
            else if (stub2Station1==3) h_bx3->Fill(stub2Bx1,aweight);
            else if (stub3Station1==3) h_bx3->Fill(stub3Bx1,aweight);
            else if (stub4Station1==3) h_bx3->Fill(stub4Bx1,aweight);
            if (stub1Station1==4) h_bx4->Fill(stub1Bx1,aweight);
            else if (stub2Station1==4) h_bx4->Fill(stub2Bx1,aweight);
            else if (stub3Station1==4) h_bx4->Fill(stub3Bx1,aweight);
            else if (stub4Station1==4) h_bx4->Fill(stub4Bx1,aweight);

            if (stub1Station2==1) h_bx1->Fill(stub1Bx1,aweight);
            else if (stub2Station2==1) h_bx1->Fill(stub2Bx1,aweight);
            else if (stub3Station2==1) h_bx1->Fill(stub3Bx1,aweight);
            else if (stub4Station2==1) h_bx1->Fill(stub4Bx1,aweight);
            if (stub1Station2==2) h_bx2->Fill(stub1Bx1,aweight);
            else if (stub2Station2==2) h_bx2->Fill(stub2Bx1,aweight);
            else if (stub3Station2==2) h_bx2->Fill(stub3Bx1,aweight);
            else if (stub4Station2==2) h_bx2->Fill(stub4Bx1,aweight);
            if (stub1Station2==3) h_bx3->Fill(stub1Bx1,aweight);
            else if (stub2Station2==3) h_bx3->Fill(stub2Bx1,aweight);
            else if (stub3Station2==3) h_bx3->Fill(stub3Bx1,aweight);
            else if (stub4Station2==3) h_bx3->Fill(stub4Bx1,aweight);
            if (stub1Station2==4) h_bx4->Fill(stub1Bx1,aweight);
            else if (stub2Station2==4) h_bx4->Fill(stub2Bx1,aweight);
            else if (stub3Station2==4) h_bx4->Fill(stub3Bx1,aweight);
            else if (stub4Station2==4) h_bx4->Fill(stub4Bx1,aweight);

            if (nstub1==4) h_station->Fill(0.,aweight);
            if (nstub1==3 and stub1Station1==4 and stub2Station1==3 and stub3Station1==2) h_station->Fill(1,aweight);
            if (nstub1==3 and stub1Station1==4 and stub2Station1==3 and stub3Station1==1) h_station->Fill(2,aweight);
            if (nstub1==3 and stub1Station1==4 and stub2Station1==2 and stub3Station1==1) h_station->Fill(3,aweight);
            if (nstub1==3 and stub1Station1==3 and stub2Station1==2 and stub3Station1==1) h_station->Fill(4,aweight);

            if (nstub2==4) h_station->Fill(0.,aweight);
            if (nstub2==3 and stub1Station2==4 and stub2Station2==3 and stub3Station2==2) h_station->Fill(1,aweight);
            if (nstub2==3 and stub1Station2==4 and stub2Station2==3 and stub3Station2==1) h_station->Fill(2,aweight);
            if (nstub2==3 and stub1Station2==4 and stub2Station2==2 and stub3Station2==1) h_station->Fill(3,aweight);
            if (nstub2==3 and stub1Station2==3 and stub2Station2==2 and stub3Station2==1) h_station->Fill(4,aweight);

            if (pt1>ptmin and pt1<=ptmax) h_nstub_OS->Fill(nstub1, aweight);
            if (pt2>ptmin and pt2<=ptmax) h_nstub_OS->Fill(nstub2, aweight);

            h_pt_OS->Fill(pt1, aweight);
            h_pt_OS->Fill(pt2, aweight);

            h_met_OS->Fill(met, aweight);
         }
      }

      else{
         h_mmumu_SS->Fill(new_mmumu, aweight);
         if (new_mmumu>70){
            if (pt1>ptmin and pt2<=ptmax) h_nstub_SS->Fill(nstub1, aweight);
            if (pt2>ptmin and pt2<=ptmax) h_nstub_SS->Fill(nstub2, aweight);
            h_pt_SS->Fill(pt1, aweight);
            h_pt_SS->Fill(pt2, aweight);
            h_met_SS->Fill(met, aweight);
         }
      }
   } // end of loop over events

   TFile *fout = TFile::Open(output.c_str(), "RECREATE");
   fout->cd();

   h_bxspan->Write();
   h_bx1->Write();
   h_bx2->Write();
   h_bx3->Write();
   h_bx4->Write();
   h_station->Write();
   h_reso->Write();
   h_reso_pt15to25->Write();
   h_reso_pt25to35->Write();
   h_reso_pt35to45->Write();
   h_reso_pt45to55->Write();
   h_reso_ptgt55->Write();

   TDirectory* dir1=fout->mkdir("OS");
   dir1->cd();
   h_mmumu_OS->SetName(name.c_str());
   if (name=="DY") h_mmumu_OS->SetName("DY");
   h_mmumu_OS->Write();

   TDirectory* dir2=fout->mkdir("SS");
   dir2->cd();
   h_mmumu_SS->SetName(name.c_str());
   if (name=="DY") h_mmumu_SS->SetName("DY");
   h_mmumu_SS->Write();

   if (muobject=="KBMTF"){
      TDirectory* dir3=fout->mkdir("met_OS");
      dir3->cd();
      h_met_OS->SetName(name.c_str());
      h_met_OS->Write();

      TDirectory* dir4=fout->mkdir("met_SS");
      dir4->cd();
      h_met_SS->SetName(name.c_str());
      h_met_SS->Write();

      TDirectory* dir5=fout->mkdir("nstub_OS");
      dir5->cd();
      h_nstub_OS->SetName(name.c_str());
      h_nstub_OS->Write();

      TDirectory* dir6=fout->mkdir("nstub_SS");
      dir6->cd();
      h_nstub_SS->SetName(name.c_str());
      h_nstub_SS->Write();

      TDirectory* dir7=fout->mkdir("pt_OS");
      dir7->cd();
      h_pt_OS->SetName(name.c_str());
      h_pt_OS->Write();

      TDirectory* dir8=fout->mkdir("pt_SS");
      dir8->cd();
      h_pt_SS->SetName(name.c_str());
      h_pt_SS->Write();
   }

   cout<< h_mmumu_OS->Integral()<<endl;

   fout->Close();
}

