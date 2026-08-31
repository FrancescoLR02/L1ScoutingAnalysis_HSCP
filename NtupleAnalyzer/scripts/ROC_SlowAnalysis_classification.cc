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
#include <TRandom3.h>
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
#include "tr_Tree.h"
#include "myHelper.h"
#include "fiducial_weight.h"
#include <TProfile.h>

using namespace std;

//!  g++ -O3 ROC_SlowAnalysis_classification.cc -o ROC_SlowAnalysis_classification.exe $(root-config --cflags --glibs)

int main(int argc, char** argv) {

    std::string input = *(argv + 1);
    std::string output = *(argv + 2);
    std::string name = *(argv + 3);
    std::string region = *(argv + 4);

    TFile *f_Double = TFile::Open(input.c_str(), "READ");
    cout<<"XXXXXXXXXXXXX "<<input.c_str()<<" XXXXXXXXXXXX"<<endl;
    TTree *arbre = (TTree*) f_Double->Get("Events");

    arbre->SetBranchAddress("bunchCrossing", &bunchCrossing);
    if(name != "data_obs") arbre->SetBranchAddress("ngen", &ngen);
    if(name == "data_obs") arbre->SetBranchAddress("is_colliding", &is_colliding);
    if(name == "data_obs") arbre->SetBranchAddress("is_earlier_colliding", &is_earlier_colliding);
    arbre->SetBranchAddress("nL1KBMTFSkimmed", &nL1KBMTFSkimmed);
    arbre->SetBranchAddress("nstub1", &nstub1);
    arbre->SetBranchAddress("genbeta1", &genbeta1);
    arbre->SetBranchAddress("geneta1", &geneta1);
    arbre->SetBranchAddress("genpt1", &genpt1);
    arbre->SetBranchAddress("nstub2", &nstub2);
    arbre->SetBranchAddress("genbeta2", &genbeta2);
    arbre->SetBranchAddress("geneta2", &geneta2);
    arbre->SetBranchAddress("genpt2", &genpt2);
    arbre->SetBranchAddress("idx2", &idx2);
    arbre->SetBranchAddress("bxspread1", &bxspread1);
    arbre->SetBranchAddress("bxspread2", &bxspread2);
    arbre->SetBranchAddress("stationspread1", &stationspread1);
    arbre->SetBranchAddress("stationspread2", &stationspread2);
    arbre->SetBranchAddress("firstbx1", &firstbx1);
    arbre->SetBranchAddress("firstbx2", &firstbx2);
    arbre->SetBranchAddress("charge1", &charge1);
    arbre->SetBranchAddress("charge2", &charge2);
    arbre->SetBranchAddress("dxy1", &dxy1);
    arbre->SetBranchAddress("dxy2", &dxy2);
    arbre->SetBranchAddress("qual1", &qual1);
    arbre->SetBranchAddress("qual2", &qual2);
    arbre->SetBranchAddress("pt1", &pt1);
    arbre->SetBranchAddress("pt2", &pt2);
    arbre->SetBranchAddress("eta1", &eta1);
    arbre->SetBranchAddress("eta2", &eta2);
    arbre->SetBranchAddress("phi1", &phi1);
    arbre->SetBranchAddress("processor1", &processor1);
    arbre->SetBranchAddress("phi2", &phi2);
    arbre->SetBranchAddress("HwK1", &HwK1);
    if(name != "data_obs") arbre->SetBranchAddress("L1MET_pt", &L1MET_pt);


   if(name == "data_obs"){
      arbre->SetBranchAddress("met_bx0", &met_bx0);
      arbre->SetBranchAddress("met_bxm1", &met_bxm1);
      arbre->SetBranchAddress("met_bxm2", &met_bxm2);
      arbre->SetBranchAddress("met_bxm3", &met_bxm3);
      arbre->SetBranchAddress("met_bxm4", &met_bxm4);
      arbre->SetBranchAddress("met_bxm5", &met_bxm5);
    }

    TH1F* h_metbx0_mu0 = new TH1F("h_metbx0_mu0", "h_metbx0_mu0", 40,0,200); h_metbx0_mu0->Sumw2();
    TH1F* h_metbxm1_mu0 = new TH1F("h_metbxm1_mu0", "h_metbxm1_mu0", 40,0,200); h_metbxm1_mu0->Sumw2();
    TH1F* h_metbxm2_mu0 = new TH1F("h_metbxm2_mu0", "h_metbxm2_mu0", 40,0,200); h_metbxm2_mu0->Sumw2();
    TH1F* h_metbxm3_mu0 = new TH1F("h_metbxm3_mu0", "h_metbxm3_mu0", 40,0,200); h_metbxm3_mu0->Sumw2();
    TH1F* h_metbxm4_mu0 = new TH1F("h_metbxm4_mu0", "h_metbxm4_mu0", 40,0,200); h_metbxm4_mu0->Sumw2();
    TH1F* h_metbxm5_mu0 = new TH1F("h_metbxm5_mu0", "h_metbxm5_mu0", 40,0,200); h_metbxm5_mu0->Sumw2();

    TH1F* h_metbx0_mu0fail = new TH1F("h_metbx0_mu0fail", "h_metbx0_mu0fail", 40,0,200); h_metbx0_mu0fail->Sumw2();
    TH1F* h_metbxm1_mu0fail = new TH1F("h_metbxm1_mu0fail", "h_metbxm1_mu0fail", 40,0,200); h_metbxm1_mu0fail->Sumw2();
    TH1F* h_metbxm2_mu0fail = new TH1F("h_metbxm2_mu0fail", "h_metbxm2_mu0fail", 40,0,200); h_metbxm2_mu0fail->Sumw2();
    TH1F* h_metbxm3_mu0fail = new TH1F("h_metbxm3_mu0fail", "h_metbxm3_mu0fail", 40,0,200); h_metbxm3_mu0fail->Sumw2();
    TH1F* h_metbxm4_mu0fail = new TH1F("h_metbxm4_mu0fail", "h_metbxm4_mu0fail", 40,0,200); h_metbxm4_mu0fail->Sumw2();
    TH1F* h_metbxm5_mu0fail = new TH1F("h_metbxm5_mu0fail", "h_metbxm5_mu0fail", 40,0,200); h_metbxm5_mu0fail->Sumw2();


   float bins_lowpt[] = {15,20,25, 30,40, 50};
   int  binnum_lowpt = sizeof(bins_lowpt)/sizeof(Float_t) - 1;

   float bins_mediumpt[] = {50, 100,  150, 200, 300, 400, 500};
   int  binnum_mediumpt = sizeof(bins_mediumpt)/sizeof(Float_t) - 1;

   float bins_highpt[] = {150, 200, 250, 300, 350, 450, 550, 650, 750};
   int  binnum_highpt = sizeof(bins_highpt)/sizeof(Float_t) - 1;

   float bins_bx123[] = {100, 200, 300, 400, 500};                       // cats 2, 9
   int  binnum_bx123 = sizeof(bins_bx123)/sizeof(Float_t) - 1;


   TH1F* h_stub4_bx1234 = new TH1F("h_stub4_bx1234", "h_stub4_bx1234", binnum_lowpt, bins_lowpt); h_stub4_bx1234->Sumw2();
   TH1F* h_stub4_bx1234_wrong = new TH1F("h_stub4_bx1234_wrong", "h_stub4_bx1234_wrong", binnum_lowpt, bins_lowpt); h_stub4_bx1234_wrong->Sumw2();
   TH1F* h_stub4_bx1234_wrongU = new TH1F("h_stub4_bx1234_wrongU", "h_stub4_bx1234_wrongU", binnum_lowpt, bins_lowpt); h_stub4_bx1234_wrongU->Sumw2();
   TH1F* h_stub4_bx1234_wrongD = new TH1F("h_stub4_bx1234_wrongD", "h_stub4_bx1234_wrongD", binnum_lowpt, bins_lowpt); h_stub4_bx1234_wrongD->Sumw2();
   TH1F* h_stub4_bx1234_fail = new TH1F("h_stub4_bx1234_fail", "h_stub4_bx1234_fail", binnum_lowpt, bins_lowpt); h_stub4_bx1234_fail->Sumw2();
   TH1F* h_stub4_bx1234_wrong_fail = new TH1F("h_stub4_bx1234_wrong_fail", "h_stub4_bx1234_wrong_fail", binnum_lowpt, bins_lowpt); h_stub4_bx1234_wrong_fail->Sumw2();
   TH1F* h_stub4_bx1234_wrongU_fail = new TH1F("h_stub4_bx1234_wrongU_fail", "h_stub4_bx1234_wrongU_fail", binnum_lowpt, bins_lowpt); h_stub4_bx1234_wrongU_fail->Sumw2();
   TH1F* h_stub4_bx1234_wrongD_fail = new TH1F("h_stub4_bx1234_wrongD_fail", "h_stub4_bx1234_wrongD_fail", binnum_lowpt, bins_lowpt); h_stub4_bx1234_wrongD_fail->Sumw2();

    TH1F* shape_stub4_bx1234_3120 = new TH1F("shape_stub4_bx1234_3120", "shape_stub4_bx1234_3120", binnum_lowpt, bins_lowpt); shape_stub4_bx1234_3120->Sumw2();
    TH1F* shape_stub4_bx1234_3201 = new TH1F("shape_stub4_bx1234_3201", "shape_stub4_bx1234_3201", binnum_lowpt, bins_lowpt); shape_stub4_bx1234_3201->Sumw2();
    TH1F* shape_stub4_bx1234_3102 = new TH1F("shape_stub4_bx1234_3102", "shape_stub4_bx1234_3102", binnum_lowpt, bins_lowpt); shape_stub4_bx1234_3102->Sumw2();
    TH1F* shape_stub4_bx1234_3012 = new TH1F("shape_stub4_bx1234_3012", "shape_stub4_bx1234_3012", binnum_lowpt, bins_lowpt); shape_stub4_bx1234_3012->Sumw2();
    TH1F* shape_stub4_bx1234_3021 = new TH1F("shape_stub4_bx1234_3021", "shape_stub4_bx1234_3021", binnum_lowpt, bins_lowpt); shape_stub4_bx1234_3021->Sumw2();
    TH1F* shape_stub4_bx1234_1230 = new TH1F("shape_stub4_bx1234_1230", "shape_stub4_bx1234_1230", binnum_lowpt, bins_lowpt); shape_stub4_bx1234_1230->Sumw2();
    TH1F* shape_stub4_bx1234_1320 = new TH1F("shape_stub4_bx1234_1320", "shape_stub4_bx1234_1320", binnum_lowpt, bins_lowpt); shape_stub4_bx1234_1320->Sumw2();
    TH1F* shape_stub4_bx1234_2130 = new TH1F("shape_stub4_bx1234_2130", "shape_stub4_bx1234_2130", binnum_lowpt, bins_lowpt); shape_stub4_bx1234_2130->Sumw2();
    TH1F* shape_stub4_bx1234_2310 = new TH1F("shape_stub4_bx1234_2310", "shape_stub4_bx1234_2310", binnum_lowpt, bins_lowpt); shape_stub4_bx1234_2310->Sumw2();


   TH1F* h_stub3_bx124_slow = new TH1F("h_stub3_bx124_slow", "h_stub3_bx124_slow", binnum_bx123, bins_bx123); h_stub3_bx124_slow->Sumw2();
   TH1F* h_stub3_bx124_slow_wrong = new TH1F("h_stub3_bx124_slow_wrong", "h_stub3_bx124_slow_wrong", binnum_bx123, bins_bx123); h_stub3_bx124_slow_wrong->Sumw2();
   TH1F* h_stub3_bx124_slow_wrongU = new TH1F("h_stub3_bx124_slow_wrongU", "h_stub3_bx124_slow_wrongU", binnum_bx123, bins_bx123); h_stub3_bx124_slow_wrongU->Sumw2();
   TH1F* h_stub3_bx124_slow_wrongD = new TH1F("h_stub3_bx124_slow_wrongD", "h_stub3_bx124_slow_wrongD", binnum_bx123, bins_bx123); h_stub3_bx124_slow_wrongD->Sumw2();
   TH1F* h_stub3_bx124_slow_fail = new TH1F("h_stub3_bx124_slow_fail", "h_stub3_bx124_slow_fail", binnum_bx123, bins_bx123); h_stub3_bx124_slow_fail->Sumw2();
   TH1F* h_stub3_bx124_slow_wrong_fail = new TH1F("h_stub3_bx124_slow_wrong_fail", "h_stub3_bx124_slow_wrong_fail", binnum_bx123, bins_bx123); h_stub3_bx124_slow_wrong_fail->Sumw2();
   TH1F* h_stub3_bx124_slow_wrongU_fail = new TH1F("h_stub3_bx124_slow_wrongU_fail", "h_stub3_bx124_slow_wrongU_fail", binnum_bx123, bins_bx123); h_stub3_bx124_slow_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx124_slow_wrongD_fail = new TH1F("h_stub3_bx124_slow_wrongD_fail", "h_stub3_bx124_slow_wrongD_fail", binnum_bx123, bins_bx123); h_stub3_bx124_slow_wrongD_fail->Sumw2();

   TH1F* h_stub3_bx124_fast = new TH1F("h_stub3_bx124_fast", "h_stub3_bx124_fast", binnum_mediumpt, bins_mediumpt); h_stub3_bx124_fast->Sumw2();
   TH1F* h_stub3_bx124_fast_wrong = new TH1F("h_stub3_bx124_fast_wrong", "h_stub3_bx124_fast_wrong", binnum_mediumpt, bins_mediumpt); h_stub3_bx124_fast_wrong->Sumw2();
   TH1F* h_stub3_bx124_fast_wrongU = new TH1F("h_stub3_bx124_fast_wrongU", "h_stub3_bx124_fast_wrongU", binnum_mediumpt, bins_mediumpt); h_stub3_bx124_fast_wrongU->Sumw2();
   TH1F* h_stub3_bx124_fast_wrongD = new TH1F("h_stub3_bx124_fast_wrongD", "h_stub3_bx124_fast_wrongD", binnum_mediumpt, bins_mediumpt); h_stub3_bx124_fast_wrongD->Sumw2();
   TH1F* h_stub3_bx124_fast_fail = new TH1F("h_stub3_bx124_fast_fail", "h_stub3_bx124_fast_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx124_fast_fail->Sumw2();
   TH1F* h_stub3_bx124_fast_wrong_fail = new TH1F("h_stub3_bx124_fast_wrong_fail", "h_stub3_bx124_fast_wrong_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx124_fast_wrong_fail->Sumw2();
   TH1F* h_stub3_bx124_fast_wrongU_fail = new TH1F("h_stub3_bx124_fast_wrongU_fail", "h_stub3_bx124_fast_wrongU_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx124_fast_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx124_fast_wrongD_fail = new TH1F("h_stub3_bx124_fast_wrongD_fail", "h_stub3_bx124_fast_wrongD_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx124_fast_wrongD_fail->Sumw2();

   TH1F* h_stub3_bx123_slow = new TH1F("h_stub3_bx123_slow", "h_stub3_bx123_slow", binnum_bx123, bins_bx123); h_stub3_bx123_slow->Sumw2();
   TH1F* h_stub3_bx123_slow_wrong = new TH1F("h_stub3_bx123_slow_wrong", "h_stub3_bx123_slow_wrong", binnum_bx123, bins_bx123); h_stub3_bx123_slow_wrong->Sumw2();
   TH1F* h_stub3_bx123_slow_wrongU = new TH1F("h_stub3_bx123_slow_wrongU", "h_stub3_bx123_slow_wrongU", binnum_bx123, bins_bx123); h_stub3_bx123_slow_wrongU->Sumw2();
   TH1F* h_stub3_bx123_slow_wrongD = new TH1F("h_stub3_bx123_slow_wrongD", "h_stub3_bx123_slow_wrongD", binnum_bx123, bins_bx123); h_stub3_bx123_slow_wrongD->Sumw2();
   TH1F* h_stub3_bx123_slow_fail = new TH1F("h_stub3_bx123_slow_fail", "h_stub3_bx123_slow_fail", binnum_bx123, bins_bx123); h_stub3_bx123_slow_fail->Sumw2();
   TH1F* h_stub3_bx123_slow_wrong_fail = new TH1F("h_stub3_bx123_slow_wrong_fail", "h_stub3_bx123_slow_wrong_fail", binnum_bx123, bins_bx123); h_stub3_bx123_slow_wrong_fail->Sumw2();
   TH1F* h_stub3_bx123_slow_wrongU_fail = new TH1F("h_stub3_bx123_slow_wrongU_fail", "h_stub3_bx123_slow_wrongU_fail", binnum_bx123, bins_bx123); h_stub3_bx123_slow_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx123_slow_wrongD_fail = new TH1F("h_stub3_bx123_slow_wrongD_fail", "h_stub3_bx123_slow_wrongD_fail", binnum_bx123, bins_bx123); h_stub3_bx123_slow_wrongD_fail->Sumw2();

   TH1F* h_stub3_bx123_fast = new TH1F("h_stub3_bx123_fast", "h_stub3_bx123_fast", binnum_mediumpt, bins_mediumpt); h_stub3_bx123_fast->Sumw2();
   TH1F* h_stub3_bx123_fast_wrong = new TH1F("h_stub3_bx123_fast_wrong", "h_stub3_bx123_fast_wrong", binnum_mediumpt, bins_mediumpt); h_stub3_bx123_fast_wrong->Sumw2();
   TH1F* h_stub3_bx123_fast_wrongU = new TH1F("h_stub3_bx123_fast_wrongU", "h_stub3_bx123_fast_wrongU", binnum_mediumpt, bins_mediumpt); h_stub3_bx123_fast_wrongU->Sumw2();
   TH1F* h_stub3_bx123_fast_wrongD = new TH1F("h_stub3_bx123_fast_wrongD", "h_stub3_bx123_fast_wrongD", binnum_mediumpt, bins_mediumpt); h_stub3_bx123_fast_wrongD->Sumw2();
   TH1F* h_stub3_bx123_fast_fail = new TH1F("h_stub3_bx123_fast_fail", "h_stub3_bx123_fast_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx123_fast_fail->Sumw2();
   TH1F* h_stub3_bx123_fast_wrong_fail = new TH1F("h_stub3_bx123_fast_wrong_fail", "h_stub3_bx123_fast_wrong_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx123_fast_wrong_fail->Sumw2();
   TH1F* h_stub3_bx123_fast_wrongU_fail = new TH1F("h_stub3_bx123_fast_wrongU_fail", "h_stub3_bx123_fast_wrongU_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx123_fast_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx123_fast_wrongD_fail = new TH1F("h_stub3_bx123_fast_wrongD_fail", "h_stub3_bx123_fast_wrongD_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx123_fast_wrongD_fail->Sumw2();

   TH1F* h_stub4_bx123 = new TH1F("h_stub4_bx123", "h_stub4_bx123", binnum_bx123, bins_bx123); h_stub4_bx123->Sumw2();
   TH1F* h_stub4_bx123_wrong = new TH1F("h_stub4_bx123_wrong", "h_stub4_bx123_wrong", binnum_bx123, bins_bx123); h_stub4_bx123_wrong->Sumw2();
   TH1F* h_stub4_bx123_wrongU = new TH1F("h_stub4_bx123_wrongU", "h_stub4_bx123_wrongU", binnum_bx123, bins_bx123); h_stub4_bx123_wrongU->Sumw2();
   TH1F* h_stub4_bx123_wrongD = new TH1F("h_stub4_bx123_wrongD", "h_stub4_bx123_wrongD", binnum_bx123, bins_bx123); h_stub4_bx123_wrongD->Sumw2();
   TH1F* h_stub4_bx123_neutron = new TH1F("h_stub4_bx123_neutron", "h_stub4_bx123_neutron", binnum_bx123, bins_bx123); h_stub4_bx123_neutron->Sumw2();
   TH1F* h_stub4_bx123_fail = new TH1F("h_stub4_bx123_fail", "h_stub4_bx123_fail", binnum_bx123, bins_bx123); h_stub4_bx123_fail->Sumw2();
   TH1F* h_stub4_bx123_wrong_fail = new TH1F("h_stub4_bx123_wrong_fail", "h_stub4_bx123_wrong_fail", binnum_bx123, bins_bx123); h_stub4_bx123_wrong_fail->Sumw2();
   TH1F* h_stub4_bx123_wrongU_fail = new TH1F("h_stub4_bx123_wrongU_fail", "h_stub4_bx123_wrongU_fail", binnum_bx123, bins_bx123); h_stub4_bx123_wrongU_fail->Sumw2();
   TH1F* h_stub4_bx123_wrongD_fail = new TH1F("h_stub4_bx123_wrongD_fail", "h_stub4_bx123_wrongD_fail", binnum_bx123, bins_bx123); h_stub4_bx123_wrongD_fail->Sumw2();

   TH1F* h_stub4_bx124 = new TH1F("h_stub4_bx124", "h_stub4_bx124", binnum_mediumpt, bins_mediumpt); h_stub4_bx124->Sumw2();
   TH1F* h_stub4_bx124_wrong = new TH1F("h_stub4_bx124_wrong", "h_stub4_bx124_wrong", binnum_mediumpt, bins_mediumpt); h_stub4_bx124_wrong->Sumw2();
   TH1F* h_stub4_bx124_wrongU = new TH1F("h_stub4_bx124_wrongU", "h_stub4_bx124_wrongU", binnum_mediumpt, bins_mediumpt); h_stub4_bx124_wrongU->Sumw2();
   TH1F* h_stub4_bx124_wrongD = new TH1F("h_stub4_bx124_wrongD", "h_stub4_bx124_wrongD", binnum_mediumpt, bins_mediumpt); h_stub4_bx124_wrongD->Sumw2();
   TH1F* h_stub4_bx124_fail = new TH1F("h_stub4_bx124_fail", "h_stub4_bx124_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx124_fail->Sumw2();
   TH1F* h_stub4_bx124_wrong_fail = new TH1F("h_stub4_bx124_wrong_fail", "h_stub4_bx124_wrong_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx124_wrong_fail->Sumw2();
   TH1F* h_stub4_bx124_wrongU_fail = new TH1F("h_stub4_bx124_wrongU_fail", "h_stub4_bx124_wrongU_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx124_wrongU_fail->Sumw2();
   TH1F* h_stub4_bx124_wrongD_fail = new TH1F("h_stub4_bx124_wrongD_fail", "h_stub4_bx124_wrongD_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx124_wrongD_fail->Sumw2();

   TH1F* h_stub4_bx1122 = new TH1F("h_stub4_bx1122", "h_stub4_bx1122", binnum_highpt, bins_highpt); h_stub4_bx1122->Sumw2();
   TH1F* h_stub4_bx1122_wrong = new TH1F("h_stub4_bx1122_wrong", "h_stub4_bx1122_wrong", binnum_highpt, bins_highpt); h_stub4_bx1122_wrong->Sumw2();
   TH1F* h_stub4_bx1122_wrongU = new TH1F("h_stub4_bx1122_wrongU", "h_stub4_bx1122_wrongU", binnum_highpt, bins_highpt); h_stub4_bx1122_wrongU->Sumw2();
   TH1F* h_stub4_bx1122_wrongD = new TH1F("h_stub4_bx1122_wrongD", "h_stub4_bx1122_wrongD", binnum_highpt, bins_highpt); h_stub4_bx1122_wrongD->Sumw2();
   TH1F* h_stub4_bx1122_fail = new TH1F("h_stub4_bx1122_fail", "h_stub4_bx1122_fail", binnum_highpt, bins_highpt); h_stub4_bx1122_fail->Sumw2();
   TH1F* h_stub4_bx1122_wrong_fail = new TH1F("h_stub4_bx1122_wrong_fail", "h_stub4_bx1122_wrong_fail", binnum_highpt, bins_highpt); h_stub4_bx1122_wrong_fail->Sumw2();
   TH1F* h_stub4_bx1122_wrongU_fail = new TH1F("h_stub4_bx1122_wrongU_fail", "h_stub4_bx1122_wrongU_fail", binnum_highpt, bins_highpt); h_stub4_bx1122_wrongU_fail->Sumw2();
   TH1F* h_stub4_bx1122_wrongD_fail = new TH1F("h_stub4_bx1122_wrongD_fail", "h_stub4_bx1122_wrongD_fail", binnum_highpt, bins_highpt); h_stub4_bx1122_wrongD_fail->Sumw2();

   TH1F* h_stub4_bx1122_2tracks = new TH1F("h_stub4_bx1122_2tracks", "h_stub4_bx1122_2tracks", binnum_mediumpt, bins_mediumpt); h_stub4_bx1122_2tracks->Sumw2();
   TH1F* h_stub4_bx1122_2tracks_wrong = new TH1F("h_stub4_bx1122_2tracks_wrong", "h_stub4_bx1122_2tracks_wrong", binnum_mediumpt, bins_mediumpt); h_stub4_bx1122_2tracks_wrong->Sumw2();
   TH1F* h_stub4_bx1122_2tracks_wrongU = new TH1F("h_stub4_bx1122_2tracks_wrongU", "h_stub4_bx1122_2tracks_wrongU", binnum_mediumpt, bins_mediumpt); h_stub4_bx1122_2tracks_wrongU->Sumw2();
   TH1F* h_stub4_bx1122_2tracks_wrongD = new TH1F("h_stub4_bx1122_2tracks_wrongD", "h_stub4_bx1122_2tracks_wrongD", binnum_mediumpt, bins_mediumpt); h_stub4_bx1122_2tracks_wrongD->Sumw2();
   TH1F* h_stub4_bx1122_2tracks_fail = new TH1F("h_stub4_bx1122_2tracks_fail", "h_stub4_bx1122_2tracks_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1122_2tracks_fail->Sumw2();
   TH1F* h_stub4_bx1122_2tracks_wrong_fail = new TH1F("h_stub4_bx1122_2tracks_wrong_fail", "h_stub4_bx1122_2tracks_wrong_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1122_2tracks_wrong_fail->Sumw2();
   TH1F* h_stub4_bx1122_2tracks_wrongU_fail = new TH1F("h_stub4_bx1122_2tracks_wrongU_fail", "h_stub4_bx1122_2tracks_wrongU_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1122_2tracks_wrongU_fail->Sumw2();
   TH1F* h_stub4_bx1122_2tracks_wrongD_fail = new TH1F("h_stub4_bx1122_2tracks_wrongD_fail", "h_stub4_bx1122_2tracks_wrongD_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1122_2tracks_wrongD_fail->Sumw2();

   TH1F* h_stub4_bx1112 = new TH1F("h_stub4_bx1112", "h_stub4_bx1112", binnum_highpt, bins_highpt); h_stub4_bx1112->Sumw2();
   TH1F* h_stub4_bx1112_wrong = new TH1F("h_stub4_bx1112_wrong", "h_stub4_bx1112_wrong", binnum_highpt, bins_highpt); h_stub4_bx1112_wrong->Sumw2();
   TH1F* h_stub4_bx1112_wrongU = new TH1F("h_stub4_bx1112_wrongU", "h_stub4_bx1112_wrongU", binnum_highpt, bins_highpt); h_stub4_bx1112_wrongU->Sumw2();
   TH1F* h_stub4_bx1112_wrongD = new TH1F("h_stub4_bx1112_wrongD", "h_stub4_bx1112_wrongD", binnum_highpt, bins_highpt); h_stub4_bx1112_wrongD->Sumw2();
   TH1F* h_stub4_bx1112_fail = new TH1F("h_stub4_bx1112_fail", "h_stub4_bx1112_fail", binnum_highpt, bins_highpt); h_stub4_bx1112_fail->Sumw2();
   TH1F* h_stub4_bx1112_wrong_fail = new TH1F("h_stub4_bx1112_wrong_fail", "h_stub4_bx1112_wrong_fail", binnum_highpt, bins_highpt); h_stub4_bx1112_wrong_fail->Sumw2();
   TH1F* h_stub4_bx1112_wrongU_fail = new TH1F("h_stub4_bx1112_wrongU_fail", "h_stub4_bx1112_wrongU_fail", binnum_highpt, bins_highpt); h_stub4_bx1112_wrongU_fail->Sumw2();
   TH1F* h_stub4_bx1112_wrongD_fail = new TH1F("h_stub4_bx1112_wrongD_fail", "h_stub4_bx1112_wrongD_fail", binnum_highpt, bins_highpt); h_stub4_bx1112_wrongD_fail->Sumw2();

   TH1F* h_stub4_bx1112_2tracks = new TH1F("h_stub4_bx1112_2tracks", "h_stub4_bx1112_2tracks", binnum_mediumpt, bins_mediumpt); h_stub4_bx1112_2tracks->Sumw2();
   TH1F* h_stub4_bx1112_2tracks_wrong = new TH1F("h_stub4_bx1112_2tracks_wrong", "h_stub4_bx1112_2tracks_wrong", binnum_mediumpt, bins_mediumpt); h_stub4_bx1112_2tracks_wrong->Sumw2();
   TH1F* h_stub4_bx1112_2tracks_wrongU = new TH1F("h_stub4_bx1112_2tracks_wrongU", "h_stub4_bx1112_2tracks_wrongU", binnum_mediumpt, bins_mediumpt); h_stub4_bx1112_2tracks_wrongU->Sumw2();
   TH1F* h_stub4_bx1112_2tracks_wrongD = new TH1F("h_stub4_bx1112_2tracks_wrongD", "h_stub4_bx1112_2tracks_wrongD", binnum_mediumpt, bins_mediumpt); h_stub4_bx1112_2tracks_wrongD->Sumw2();
   TH1F* h_stub4_bx1112_2tracks_fail = new TH1F("h_stub4_bx1112_2tracks_fail", "h_stub4_bx1112_2tracks_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1112_2tracks_fail->Sumw2();
   TH1F* h_stub4_bx1112_2tracks_wrong_fail = new TH1F("h_stub4_bx1112_2tracks_wrong_fail", "h_stub4_bx1112_2tracks_wrong_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1112_2tracks_wrong_fail->Sumw2();
   TH1F* h_stub4_bx1112_2tracks_wrongU_fail = new TH1F("h_stub4_bx1112_2tracks_wrongU_fail", "h_stub4_bx1112_2tracks_wrongU_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1112_2tracks_wrongU_fail->Sumw2();
   TH1F* h_stub4_bx1112_2tracks_wrongD_fail = new TH1F("h_stub4_bx1112_2tracks_wrongD_fail", "h_stub4_bx1112_2tracks_wrongD_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1112_2tracks_wrongD_fail->Sumw2();

   TH1F* h_stub4_bx1222 = new TH1F("h_stub4_bx1222", "h_stub4_bx1222", binnum_highpt, bins_highpt); h_stub4_bx1222->Sumw2();
   TH1F* h_stub4_bx1222_wrong = new TH1F("h_stub4_bx1222_wrong", "h_stub4_bx1222_wrong", binnum_highpt, bins_highpt); h_stub4_bx1222_wrong->Sumw2();
   TH1F* h_stub4_bx1222_wrongU = new TH1F("h_stub4_bx1222_wrongU", "h_stub4_bx1222_wrongU", binnum_highpt, bins_highpt); h_stub4_bx1222_wrongU->Sumw2();
   TH1F* h_stub4_bx1222_wrongD = new TH1F("h_stub4_bx1222_wrongD", "h_stub4_bx1222_wrongD", binnum_highpt, bins_highpt); h_stub4_bx1222_wrongD->Sumw2();
   TH1F* h_stub4_bx1222_fail = new TH1F("h_stub4_bx1222_fail", "h_stub4_bx1222_fail", binnum_highpt, bins_highpt); h_stub4_bx1222_fail->Sumw2();
   TH1F* h_stub4_bx1222_wrong_fail = new TH1F("h_stub4_bx1222_wrong_fail", "h_stub4_bx1222_wrong_fail", binnum_highpt, bins_highpt); h_stub4_bx1222_wrong_fail->Sumw2();
   TH1F* h_stub4_bx1222_wrongU_fail = new TH1F("h_stub4_bx1222_wrongU_fail", "h_stub4_bx1222_wrongU_fail", binnum_highpt, bins_highpt); h_stub4_bx1222_wrongU_fail->Sumw2();
   TH1F* h_stub4_bx1222_wrongD_fail = new TH1F("h_stub4_bx1222_wrongD_fail", "h_stub4_bx1222_wrongD_fail", binnum_highpt, bins_highpt); h_stub4_bx1222_wrongD_fail->Sumw2();

   TH1F* h_stub4_bx1222_2tracks = new TH1F("h_stub4_bx1222_2tracks", "h_stub4_bx1222_2tracks", binnum_mediumpt, bins_mediumpt); h_stub4_bx1222_2tracks->Sumw2();
   TH1F* h_stub4_bx1222_2tracks_wrong = new TH1F("h_stub4_bx1222_2tracks_wrong", "h_stub4_bx1222_2tracks_wrong", binnum_mediumpt, bins_mediumpt); h_stub4_bx1222_2tracks_wrong->Sumw2();
   TH1F* h_stub4_bx1222_2tracks_wrongU = new TH1F("h_stub4_bx1222_2tracks_wrongU", "h_stub4_bx1222_2tracks_wrongU", binnum_mediumpt, bins_mediumpt); h_stub4_bx1222_2tracks_wrongU->Sumw2();
   TH1F* h_stub4_bx1222_2tracks_wrongD = new TH1F("h_stub4_bx1222_2tracks_wrongD", "h_stub4_bx1222_2tracks_wrongD", binnum_mediumpt, bins_mediumpt); h_stub4_bx1222_2tracks_wrongD->Sumw2();
   TH1F* h_stub4_bx1222_2tracks_fail = new TH1F("h_stub4_bx1222_2tracks_fail", "h_stub4_bx1222_2tracks_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1222_2tracks_fail->Sumw2();
   TH1F* h_stub4_bx1222_2tracks_wrong_fail = new TH1F("h_stub4_bx1222_2tracks_wrong_fail", "h_stub4_bx1222_2tracks_wrong_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1222_2tracks_wrong_fail->Sumw2();
   TH1F* h_stub4_bx1222_2tracks_wrongU_fail = new TH1F("h_stub4_bx1222_2tracks_wrongU_fail", "h_stub4_bx1222_2tracks_wrongU_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1222_2tracks_wrongU_fail->Sumw2();
   TH1F* h_stub4_bx1222_2tracks_wrongD_fail = new TH1F("h_stub4_bx1222_2tracks_wrongD_fail", "h_stub4_bx1222_2tracks_wrongD_fail", binnum_mediumpt, bins_mediumpt); h_stub4_bx1222_2tracks_wrongD_fail->Sumw2();

   TH1F* h_stub3_bx112_slow = new TH1F("h_stub3_bx112_slow", "h_stub3_bx112_slow", binnum_highpt, bins_highpt); h_stub3_bx112_slow->Sumw2();
   TH1F* h_stub3_bx112_slow_wrong = new TH1F("h_stub3_bx112_slow_wrong", "h_stub3_bx112_slow_wrong", binnum_highpt, bins_highpt); h_stub3_bx112_slow_wrong->Sumw2();
   TH1F* h_stub3_bx112_slow_wrongU = new TH1F("h_stub3_bx112_slow_wrongU", "h_stub3_bx112_slow_wrongU", binnum_highpt, bins_highpt); h_stub3_bx112_slow_wrongU->Sumw2();
   TH1F* h_stub3_bx112_slow_wrongD = new TH1F("h_stub3_bx112_slow_wrongD", "h_stub3_bx112_slow_wrongD", binnum_highpt, bins_highpt); h_stub3_bx112_slow_wrongD->Sumw2();
   TH1F* h_stub3_bx112_slow_fail = new TH1F("h_stub3_bx112_slow_fail", "h_stub3_bx112_slow_fail", binnum_highpt, bins_highpt); h_stub3_bx112_slow_fail->Sumw2();
   TH1F* h_stub3_bx112_slow_wrong_fail = new TH1F("h_stub3_bx112_slow_wrong_fail", "h_stub3_bx112_slow_wrong_fail", binnum_highpt, bins_highpt); h_stub3_bx112_slow_wrong_fail->Sumw2();
   TH1F* h_stub3_bx112_slow_wrongU_fail = new TH1F("h_stub3_bx112_slow_wrongU_fail", "h_stub3_bx112_slow_wrongU_fail", binnum_highpt, bins_highpt); h_stub3_bx112_slow_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx112_slow_wrongD_fail = new TH1F("h_stub3_bx112_slow_wrongD_fail", "h_stub3_bx112_slow_wrongD_fail", binnum_highpt, bins_highpt); h_stub3_bx112_slow_wrongD_fail->Sumw2();

   TH1F* h_stub3_bx112_slow_2tracks = new TH1F("h_stub3_bx112_slow_2tracks", "h_stub3_bx112_slow_2tracks", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_slow_2tracks->Sumw2();
   TH1F* h_stub3_bx112_slow_2tracks_wrong = new TH1F("h_stub3_bx112_slow_2tracks_wrong", "h_stub3_bx112_slow_2tracks_wrong", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_slow_2tracks_wrong->Sumw2();
   TH1F* h_stub3_bx112_slow_2tracks_wrongU = new TH1F("h_stub3_bx112_slow_2tracks_wrongU", "h_stub3_bx112_slow_2tracks_wrongU", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_slow_2tracks_wrongU->Sumw2();
   TH1F* h_stub3_bx112_slow_2tracks_wrongD = new TH1F("h_stub3_bx112_slow_2tracks_wrongD", "h_stub3_bx112_slow_2tracks_wrongD", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_slow_2tracks_wrongD->Sumw2();
   TH1F* h_stub3_bx112_slow_2tracks_fail = new TH1F("h_stub3_bx112_slow_2tracks_fail", "h_stub3_bx112_slow_2tracks_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_slow_2tracks_fail->Sumw2();
   TH1F* h_stub3_bx112_slow_2tracks_wrong_fail = new TH1F("h_stub3_bx112_slow_2tracks_wrong_fail", "h_stub3_bx112_slow_2tracks_wrong_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_slow_2tracks_wrong_fail->Sumw2();
   TH1F* h_stub3_bx112_slow_2tracks_wrongU_fail = new TH1F("h_stub3_bx112_slow_2tracks_wrongU_fail", "h_stub3_bx112_slow_2tracks_wrongU_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_slow_2tracks_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx112_slow_2tracks_wrongD_fail = new TH1F("h_stub3_bx112_slow_2tracks_wrongD_fail", "h_stub3_bx112_slow_2tracks_wrongD_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_slow_2tracks_wrongD_fail->Sumw2();

   TH1F* h_stub3_bx112_fast = new TH1F("h_stub3_bx112_fast", "h_stub3_bx112_fast", binnum_highpt, bins_highpt); h_stub3_bx112_fast->Sumw2();
   TH1F* h_stub3_bx112_fast_wrong = new TH1F("h_stub3_bx112_fast_wrong", "h_stub3_bx112_fast_wrong", binnum_highpt, bins_highpt); h_stub3_bx112_fast_wrong->Sumw2();
   TH1F* h_stub3_bx112_fast_wrongU = new TH1F("h_stub3_bx112_fast_wrongU", "h_stub3_bx112_fast_wrongU", binnum_highpt, bins_highpt); h_stub3_bx112_fast_wrongU->Sumw2();
   TH1F* h_stub3_bx112_fast_wrongD = new TH1F("h_stub3_bx112_fast_wrongD", "h_stub3_bx112_fast_wrongD", binnum_highpt, bins_highpt); h_stub3_bx112_fast_wrongD->Sumw2();
   TH1F* h_stub3_bx112_fast_fail = new TH1F("h_stub3_bx112_fast_fail", "h_stub3_bx112_fast_fail", binnum_highpt, bins_highpt); h_stub3_bx112_fast_fail->Sumw2();
   TH1F* h_stub3_bx112_fast_wrong_fail = new TH1F("h_stub3_bx112_fast_wrong_fail", "h_stub3_bx112_fast_wrong_fail", binnum_highpt, bins_highpt); h_stub3_bx112_fast_wrong_fail->Sumw2();
   TH1F* h_stub3_bx112_fast_wrongU_fail = new TH1F("h_stub3_bx112_fast_wrongU_fail", "h_stub3_bx112_fast_wrongU_fail", binnum_highpt, bins_highpt); h_stub3_bx112_fast_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx112_fast_wrongD_fail = new TH1F("h_stub3_bx112_fast_wrongD_fail", "h_stub3_bx112_fast_wrongD_fail", binnum_highpt, bins_highpt); h_stub3_bx112_fast_wrongD_fail->Sumw2();

   TH1F* h_stub3_bx112_fast_2tracks = new TH1F("h_stub3_bx112_fast_2tracks", "h_stub3_bx112_fast_2tracks", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_fast_2tracks->Sumw2();
   TH1F* h_stub3_bx112_fast_2tracks_wrong = new TH1F("h_stub3_bx112_fast_2tracks_wrong", "h_stub3_bx112_fast_2tracks_wrong", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_fast_2tracks_wrong->Sumw2();
   TH1F* h_stub3_bx112_fast_2tracks_wrongU = new TH1F("h_stub3_bx112_fast_2tracks_wrongU", "h_stub3_bx112_fast_2tracks_wrongU", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_fast_2tracks_wrongU->Sumw2();
   TH1F* h_stub3_bx112_fast_2tracks_wrongD = new TH1F("h_stub3_bx112_fast_2tracks_wrongD", "h_stub3_bx112_fast_2tracks_wrongD", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_fast_2tracks_wrongD->Sumw2();
   TH1F* h_stub3_bx112_fast_2tracks_fail = new TH1F("h_stub3_bx112_fast_2tracks_fail", "h_stub3_bx112_fast_2tracks_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_fast_2tracks_fail->Sumw2();
   TH1F* h_stub3_bx112_fast_2tracks_wrong_fail = new TH1F("h_stub3_bx112_fast_2tracks_wrong_fail", "h_stub3_bx112_fast_2tracks_wrong_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_fast_2tracks_wrong_fail->Sumw2();
   TH1F* h_stub3_bx112_fast_2tracks_wrongU_fail = new TH1F("h_stub3_bx112_fast_2tracks_wrongU_fail", "h_stub3_bx112_fast_2tracks_wrongU_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_fast_2tracks_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx112_fast_2tracks_wrongD_fail = new TH1F("h_stub3_bx112_fast_2tracks_wrongD_fail", "h_stub3_bx112_fast_2tracks_wrongD_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx112_fast_2tracks_wrongD_fail->Sumw2();


   TH1F* h_stub3_bx122_slow = new TH1F("h_stub3_bx122_slow", "h_stub3_bx122_slow", binnum_highpt, bins_highpt); h_stub3_bx122_slow->Sumw2();
   TH1F* h_stub3_bx122_slow_wrong = new TH1F("h_stub3_bx122_slow_wrong", "h_stub3_bx122_slow_wrong", binnum_highpt, bins_highpt); h_stub3_bx122_slow_wrong->Sumw2();
   TH1F* h_stub3_bx122_slow_wrongU = new TH1F("h_stub3_bx122_slow_wrongU", "h_stub3_bx122_slow_wrongU", binnum_highpt, bins_highpt); h_stub3_bx122_slow_wrongU->Sumw2();
   TH1F* h_stub3_bx122_slow_wrongD = new TH1F("h_stub3_bx122_slow_wrongD", "h_stub3_bx122_slow_wrongD", binnum_highpt, bins_highpt); h_stub3_bx122_slow_wrongD->Sumw2();
   TH1F* h_stub3_bx122_slow_fail = new TH1F("h_stub3_bx122_slow_fail", "h_stub3_bx122_slow_fail", binnum_highpt, bins_highpt); h_stub3_bx122_slow_fail->Sumw2();
   TH1F* h_stub3_bx122_slow_wrong_fail = new TH1F("h_stub3_bx122_slow_wrong_fail", "h_stub3_bx122_slow_wrong_fail", binnum_highpt, bins_highpt); h_stub3_bx122_slow_wrong_fail->Sumw2();
   TH1F* h_stub3_bx122_slow_wrongU_fail = new TH1F("h_stub3_bx122_slow_wrongU_fail", "h_stub3_bx122_slow_wrongU_fail", binnum_highpt, bins_highpt); h_stub3_bx122_slow_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx122_slow_wrongD_fail = new TH1F("h_stub3_bx122_slow_wrongD_fail", "h_stub3_bx122_slow_wrongD_fail", binnum_highpt, bins_highpt); h_stub3_bx122_slow_wrongD_fail->Sumw2();

   TH1F* h_stub3_bx122_slow_2tracks = new TH1F("h_stub3_bx122_slow_2tracks", "h_stub3_bx122_slow_2tracks", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_slow_2tracks->Sumw2();
   TH1F* h_stub3_bx122_slow_2tracks_wrong = new TH1F("h_stub3_bx122_slow_2tracks_wrong", "h_stub3_bx122_slow_2tracks_wrong", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_slow_2tracks_wrong->Sumw2();
   TH1F* h_stub3_bx122_slow_2tracks_wrongU = new TH1F("h_stub3_bx122_slow_2tracks_wrongU", "h_stub3_bx122_slow_2tracks_wrongU", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_slow_2tracks_wrongU->Sumw2();
   TH1F* h_stub3_bx122_slow_2tracks_wrongD = new TH1F("h_stub3_bx122_slow_2tracks_wrongD", "h_stub3_bx122_slow_2tracks_wrongD", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_slow_2tracks_wrongD->Sumw2();
   TH1F* h_stub3_bx122_slow_2tracks_fail = new TH1F("h_stub3_bx122_slow_2tracks_fail", "h_stub3_bx122_slow_2tracks_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_slow_2tracks_fail->Sumw2();
   TH1F* h_stub3_bx122_slow_2tracks_wrong_fail = new TH1F("h_stub3_bx122_slow_2tracks_wrong_fail", "h_stub3_bx122_slow_2tracks_wrong_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_slow_2tracks_wrong_fail->Sumw2();
   TH1F* h_stub3_bx122_slow_2tracks_wrongU_fail = new TH1F("h_stub3_bx122_slow_2tracks_wrongU_fail", "h_stub3_bx122_slow_2tracks_wrongU_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_slow_2tracks_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx122_slow_2tracks_wrongD_fail = new TH1F("h_stub3_bx122_slow_2tracks_wrongD_fail", "h_stub3_bx122_slow_2tracks_wrongD_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_slow_2tracks_wrongD_fail->Sumw2();

   TH1F* h_stub3_bx122_fast = new TH1F("h_stub3_bx122_fast", "h_stub3_bx122_fast", binnum_highpt, bins_highpt); h_stub3_bx122_fast->Sumw2();
   TH1F* h_stub3_bx122_fast_wrong = new TH1F("h_stub3_bx122_fast_wrong", "h_stub3_bx122_fast_wrong", binnum_highpt, bins_highpt); h_stub3_bx122_fast_wrong->Sumw2();
   TH1F* h_stub3_bx122_fast_wrongU = new TH1F("h_stub3_bx122_fast_wrongU", "h_stub3_bx122_fast_wrongU", binnum_highpt, bins_highpt); h_stub3_bx122_fast_wrongU->Sumw2();
   TH1F* h_stub3_bx122_fast_wrongD = new TH1F("h_stub3_bx122_fast_wrongD", "h_stub3_bx122_fast_wrongD", binnum_highpt, bins_highpt); h_stub3_bx122_fast_wrongD->Sumw2();
   TH1F* h_stub3_bx122_fast_fail = new TH1F("h_stub3_bx122_fast_fail", "h_stub3_bx122_fast_fail", binnum_highpt, bins_highpt); h_stub3_bx122_fast_fail->Sumw2();
   TH1F* h_stub3_bx122_fast_wrong_fail = new TH1F("h_stub3_bx122_fast_wrong_fail", "h_stub3_bx122_fast_wrong_fail", binnum_highpt, bins_highpt); h_stub3_bx122_fast_wrong_fail->Sumw2();
   TH1F* h_stub3_bx122_fast_wrongU_fail = new TH1F("h_stub3_bx122_fast_wrongU_fail", "h_stub3_bx122_fast_wrongU_fail", binnum_highpt, bins_highpt); h_stub3_bx122_fast_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx122_fast_wrongD_fail = new TH1F("h_stub3_bx122_fast_wrongD_fail", "h_stub3_bx122_fast_wrongD_fail", binnum_highpt, bins_highpt); h_stub3_bx122_fast_wrongD_fail->Sumw2();

   TH1F* h_stub3_bx122_fast_2tracks = new TH1F("h_stub3_bx122_fast_2tracks", "h_stub3_bx122_fast_2tracks", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_fast_2tracks->Sumw2();
   TH1F* h_stub3_bx122_fast_2tracks_wrong = new TH1F("h_stub3_bx122_fast_2tracks_wrong", "h_stub3_bx122_fast_2tracks_wrong", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_fast_2tracks_wrong->Sumw2();
   TH1F* h_stub3_bx122_fast_2tracks_wrongU = new TH1F("h_stub3_bx122_fast_2tracks_wrongU", "h_stub3_bx122_fast_2tracks_wrongU", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_fast_2tracks_wrongU->Sumw2();
   TH1F* h_stub3_bx122_fast_2tracks_wrongD = new TH1F("h_stub3_bx122_fast_2tracks_wrongD", "h_stub3_bx122_fast_2tracks_wrongD", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_fast_2tracks_wrongD->Sumw2();
   TH1F* h_stub3_bx122_fast_2tracks_fail = new TH1F("h_stub3_bx122_fast_2tracks_fail", "h_stub3_bx122_fast_2tracks_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_fast_2tracks_fail->Sumw2();
   TH1F* h_stub3_bx122_fast_2tracks_wrong_fail = new TH1F("h_stub3_bx122_fast_2tracks_wrong_fail", "h_stub3_bx122_fast_2tracks_wrong_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_fast_2tracks_wrong_fail->Sumw2();
   TH1F* h_stub3_bx122_fast_2tracks_wrongU_fail = new TH1F("h_stub3_bx122_fast_2tracks_wrongU_fail", "h_stub3_bx122_fast_2tracks_wrongU_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_fast_2tracks_wrongU_fail->Sumw2();
   TH1F* h_stub3_bx122_fast_2tracks_wrongD_fail = new TH1F("h_stub3_bx122_fast_2tracks_wrongD_fail", "h_stub3_bx122_fast_2tracks_wrongD_fail", binnum_mediumpt, bins_mediumpt); h_stub3_bx122_fast_2tracks_wrongD_fail->Sumw2();



   TH1F* h_qual_nstub2 = new TH1F("h_qual_nstub2","h_qual_nstub2",4,12,16); h_qual_nstub2->Sumw2();
   TH1F* h_qual_nstub3 = new TH1F("h_qual_nstub3","h_qual_nstub3",4,12,16); h_qual_nstub3->Sumw2();
   TH1F* h_qual_nstub4 = new TH1F("h_qual_nstub4","h_qual_nstub4",4,12,16); h_qual_nstub4->Sumw2();

   TH1F* h_dxy = new TH1F("h_dxy","h_dxy",80,0,4); h_dxy->Sumw2();
   TH1F* h_dxy3stubs = new TH1F("h_dxy3stubs","h_dxy3stubs",80,0,4); h_dxy3stubs->Sumw2();
   TH1F* h_dxy4stubs = new TH1F("h_dxy4stubs","h_dxy4stubs",80,0,4); h_dxy4stubs->Sumw2();
   TH1F* h_dxy3stubs_wrong = new TH1F("h_dxy3stubs_wrong","h_dxy3stubs_wrong",80,0,4); h_dxy3stubs_wrong->Sumw2();
   TH1F* h_dxy4stubs_wrong = new TH1F("h_dxy4stubs_wrong","h_dxy4stubs_wrong",80,0,4); h_dxy4stubs_wrong->Sumw2();
   TH1F* h_phi_2BX        = new TH1F("h_phi_2BX","h_phi_2BX",50,-3.14159,3.14159); h_phi_2BX->Sumw2();
   
   TH1F* h_phi_2BX_wrong  = new TH1F("h_phi_2BX_wrong","h_phi_2BX_wrong",50,-3.14159,3.14159); h_phi_2BX_wrong->Sumw2();
   TH1F* h_phi_gt2BX      = new TH1F("h_phi_gt2BX","h_phi_gt2BX",50,-3.14159,3.14159); h_phi_gt2BX->Sumw2();
   TH1F* h_phi_gt2BX_wrong= new TH1F("h_phi_gt2BX_wrong","h_phi_gt2BX_wrong",50,-3.14159,3.14159); h_phi_gt2BX_wrong->Sumw2();
   TH1F* h_1234ordering   = new TH1F("h_1234ordering","h_1234ordering",24,0,24); h_1234ordering->Sumw2();

   TH1F* h_nstub = new TH1F("h_nstub","h_nstub",3,2,5); h_nstub->Sumw2();


   const double PHI_SHIFT = 0.0;      
   const int    NK    = 4800;
   const double KMIN  = -600.0, KMAX = 600.0;  

   TH1F* h_pt = new TH1F("h_pt","h_pt",400,12.5,1300); h_pt->Sumw2();
   TH1F* h_charge = new TH1F("h_charge","h_charge", 3, -1.5, 1.5); h_charge->Sumw2();
   TH1F* h_dxy_prompt= new TH1F("h_dxy_prompt","d_{xy}, prompt",400,0,4); h_dxy_prompt->Sumw2();
   TH2F* h2_K_processor_plus_stub3  = new TH2F("h2_K_processor_plus_stub3", "K vs #phi, q>0, 3 stubs", 12, 0, 12, NK, KMIN, KMAX); h2_K_processor_plus_stub3->Sumw2();
   TH2F* h2_K_processor_minus_stub3 = new TH2F("h2_K_processor_minus_stub3","K vs #phi, q<0, 3 stubs", 12, 0, 12, NK, KMIN, KMAX); h2_K_processor_minus_stub3->Sumw2();
   TH2F* h2_K_processor_plus_stub4  = new TH2F("h2_K_processor_plus_stub4", "K vs #phi, q>0, 4 stubs", 12, 0, 12, NK, KMIN, KMAX); h2_K_processor_plus_stub4->Sumw2();
   TH2F* h2_K_processor_minus_stub4 = new TH2F("h2_K_processor_minus_stub4","K vs #phi, q<0, 4 stubs", 12, 0, 12, NK, KMIN, KMAX); h2_K_processor_minus_stub4->Sumw2();

   TH2F* h2_K_eta_plus  = new TH2F("h2_K_eta_plus", "K vs #eta, q>0", 30,-1.0,1.0, NK, KMIN, KMAX); h2_K_eta_plus->Sumw2(); 
   TH2F* h2_K_eta_minus = new TH2F("h2_K_eta_minus","K vs #eta, q<0", 30,-1.0,1.0, NK, KMIN, KMAX); h2_K_eta_minus->Sumw2();

   TH1F* h_K = new TH1F("h_K","hw curvature K", NK, KMIN, KMAX); h_K->Sumw2();
   TProfile* p_K_vs_phi   = new TProfile("p_K_vs_phi",  "", 12, -3.1415, 3.1415);
   TProfile* p_K_vs_eta   = new TProfile("p_K_vs_eta",  "", 20, -1.0, 1.0);
   TProfile* p_K_vs_nstub = new TProfile("p_K_vs_nstub","",  3, 2, 5);
   TProfile* p_K_vs_dxy   = new TProfile("p_K_vs_dxy",  "", 40, 0, 4);

   TH1F* h_sector_plus  = new TH1F("h_sector_plus", "sector, q>0", 12, 0, 12); h_sector_plus->Sumw2();
   TH1F* h_sector_minus = new TH1F("h_sector_minus","sector, q<0", 12, 0, 12); h_sector_minus->Sumw2();

   TProfile* p_dxy_vs_phi       = new TProfile("p_dxy_vs_phi",      "", 24, -3.14159, 3.14159);
   TProfile* p_dxy_vs_phi_plus  = new TProfile("p_dxy_vs_phi_plus", "", 24, -3.14159, 3.14159);
   TProfile* p_dxy_vs_phi_minus = new TProfile("p_dxy_vs_phi_minus","", 24, -3.14159, 3.14159);
   TProfile* p_dxy_vs_eta       = new TProfile("p_dxy_vs_eta",      "", 20, -1.0, 1.0);

   TH2F* h2_dxy_vs_phi = new TH2F("h2_dxy_vs_phi","d_{xy} vs #phi", 24,-3.14159,3.14159, 400,0,4); h2_dxy_vs_phi->Sumw2();

   TH2F* h2_K_phi_plus_tightdxy  = new TH2F("h2_K_phi_plus_tightdxy", "", 12,-3.14159,3.14159, NK,KMIN,KMAX); h2_K_phi_plus_tightdxy->Sumw2(); 
   TH2F* h2_K_phi_minus_tightdxy = new TH2F("h2_K_phi_minus_tightdxy","", 12,-3.14159,3.14159, NK,KMIN,KMAX); h2_K_phi_minus_tightdxy->Sumw2();


   TH1F* h_dK_hw = new TH1F("h_dK_hw","K_{reco} - K_{gen} [LSB]", 1025, -128.125, 128.125); h_dK_hw->Sumw2();
   TProfile* p_dK_hw_vs_eta   = new TProfile("p_dK_hw_vs_eta",  "", 20, -1.0, 1.0);
   TProfile* p_dK_hw_vs_phi   = new TProfile("p_dK_hw_vs_phi",  "", 12, -3.14159, 3.14159);
   TProfile* p_dK_hw_vs_beta  = new TProfile("p_dK_hw_vs_beta", "", 20, 0, 1);
   TProfile* p_dK_hw_vs_Kgen  = new TProfile("p_dK_hw_vs_Kgen", "", 40, -400, 400); 
   TProfile* p_dK_hw_vs_dxy   = new TProfile("p_dK_hw_vs_dxy",  "", 40, 0, 4);
   TH2F* h2_K_hw_reco_gen = new TH2F("h2_K_hw_reco_gen","K_{reco} vs K_{gen} [LSB]", 200,-400,400, 200,-400,400); h2_K_hw_reco_gen->Sumw2();

   TProfile* p_dK_hw_vs_phi_plus  = new TProfile("p_dK_hw_vs_phi_plus", "", 12, -3.14159, 3.14159);
   TProfile* p_dK_hw_vs_phi_minus = new TProfile("p_dK_hw_vs_phi_minus","", 12, -3.14159, 3.14159);


   //float lumiweight= (37770.0/15.046) + (5440.0/15.060) + (11470.0/15.048);//recorded
   //    float lumiweight= (40400.0/15.0) + (5820.0/15.0) + (12320.0/15.0); //delivered
   float lumiweight=5380.0;
   //float xsweight=1.000 * lumiweight; //20000 gen events
   //if (name=="data_obs") xsweight=1.0;
   static TRandom3 randGen(1234);


   float fidweight=Get_fidweight(name,input);
   cout<<fidweight<<endl;

   if (name.find("fid0p") != std::string::npos) fidweight = 0.5*fidweight; // sum 2 samples for better stat

   long ngenl=ngen;
   if (output=="out_2024G_part1.root") ngenl=2604637102;
   else if (output=="out_2024H.root") ngenl=584267654;
   else if (output=="out_2024I_part1.root") ngenl=895013722;
   else if (name=="data_obs") ngenl=0;

   Int_t nentries_wtn = (Int_t) arbre->GetEntries();
   for (Int_t i = 0; i < nentries_wtn; i++) {
        arbre->GetEntry(i);
        if (i % 100000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
        fflush(stdout);

      //if (name!="data_obs" or is_earlier_colliding) continue; // FIXME keep only noncolliding bunches

      if (region == "VR"){
         if (name=="data_obs" and is_earlier_colliding) continue; // FIXME keep only noncolliding bunches in data, dont touch mc
      }
      
      if (region == "SR"){
         if (name=="data_obs" and !is_earlier_colliding) continue; // FIXME keep only colliding bunches
      }

      float xsweight=1.0;
      if (name!="data_obs") xsweight=(lumiweight*fidweight)/ngen;
      //if (name.find("GluinoBall") != std::string::npos or name.find("ZPrimeTo2TauPrime") != std::string::npos or name.find("HSCP") != std::string::npos) xsweight = 0.1*xsweight; //FIXME rescaling to compute limits !!!!!!!!!


      float weight=1.0;

      //###########################################################
      //############## CHECK QUALITY CRITERIA #####################
      //###########################################################
      
      if (nstub1==2) h_qual_nstub2->Fill(qual1);
      if (nstub1==3) h_qual_nstub3->Fill(qual1);
      if (nstub1==4) h_qual_nstub4->Fill(qual1);
      h_nstub->Fill(nstub1);
      if (pt2>15) h_nstub->Fill(nstub2);
      h_dxy->Fill(dxy1);
      h_pt->Fill(pt1);
      double K = charge1/pt1;

      bool pass_quality_1=true;
      if (nstub1==2 and qual1<13) pass_quality_1=false;
      if (nstub1==3 and qual1<14) pass_quality_1=false;
      if (nstub1==4 and qual1<15) pass_quality_1=false;
      //if (dxy1>0) pass_quality_1=false;
      //if (dxy1>0) continue;

      bool has_2goodtracks=false;
         if (idx2<90 and pt2>100 and dxy2<1 and qual2>12){
            if (nstub2==2 and qual2>12) has_2goodtracks=true;
            if (nstub2==3 and qual2>13) has_2goodtracks=true;
            if (nstub2==4 and qual2>14) has_2goodtracks=true;
         }
         bool has_2tracks=false;
         if (idx2<90) has_2tracks=true;

            bool is_accepted=false;

      
      
            // Muon reconstruction efficiency
      float musf=1.0;
      float aweight=1.0;
      float weight_2tracks=1.0;
         aweight = aweight*musf;
      if (name=="data_obs") aweight=1.0;
      if (name=="data_obs") weight_2tracks=1.0;

      float w1 = xsweight*aweight;
      float w2 = xsweight*aweight*weight_2tracks;
      if (name=="data_obs"){
         w1=1.0;
         w2=1.0;
      }


      bool is_slow = (stationspread1==4320 or stationspread1==3210);
      bool is_tagged=false;

      if (pt1>15 and bxspread1>0 and nstub1>2 and pass_quality_1){
         if (nstub1==4 and pt1>50 and (bxspread1==3100 or bxspread1==3110 or bxspread1==3310 or bxspread1==3200 or bxspread1==3220 or bxspread1==3320)){ is_tagged=true; }//bx124
         if (nstub1==4 and pt1>50 and (bxspread1==2100 or bxspread1==2110 or bxspread1==2210)){  is_tagged=true;}//bx123
         if (nstub1==4 and pt1>15 and bxspread1==3210){ is_tagged=true; }//bx1234
         if (nstub1==4 and pt1>50 and bxspread1==1110){ is_tagged=true; }//bx1222
         if (nstub1==4 and pt1>50 and bxspread1==1100){ is_tagged=true; }//bx1122
         if (nstub1==4 and pt1>50 and bxspread1==1000){ is_tagged=true; }//bx1112

         if (nstub1==3 and pt1>50 and is_slow and (bxspread1==3200 or bxspread1==3100)){  is_tagged=true; }//bx124 slow
         if (nstub1==3 and pt1>50 and !is_slow and (bxspread1==3200 or bxspread1==3100)){ is_tagged=true; }//bx124 fast
         if (nstub1==3 and pt1>50 and is_slow and bxspread1==2100){  is_tagged=true; }//bx123 slow
         if (nstub1==3 and pt1>50 and !is_slow and bxspread1==2100){  is_tagged=true; }//bx123 fast
         if (nstub1==3 and pt1>50 and is_slow and bxspread1==1100){  is_tagged=true; }//bx122 slow
         if (nstub1==3 and pt1>50 and !is_slow and bxspread1==1100){  is_tagged=true; }//bx122 fast
         if (nstub1==3 and pt1>50 and is_slow and bxspread1==1000){is_tagged=true; }//bx112 slow
         if (nstub1==3 and pt1>50 and !is_slow and bxspread1==1000){  is_tagged=true; }//bx112 fast

         if (is_tagged){
            if (nstub1 == 3) h_dxy3stubs->Fill(dxy1);
            if (nstub1 == 4) h_dxy4stubs->Fill(dxy1);
            h_dxy_prompt->Fill(dxy1, w1);
         }

         //if (!is_tagged) cout<<"nstub pt bxspread "<<nstub1<<" "<<pt1<<" "<<bxspread1<<endl;
      }


      
if (pass_quality_1) {

   h_K->Fill(HwK1, w1);
   p_K_vs_nstub->Fill(nstub1, HwK1, w1);
   p_K_vs_dxy->Fill(dxy1, HwK1, w1);

   // ---- alignment sample: prompt, in-time, >=3 stubs ---------------------
   // No pt cut. The quantile ladder needs the whole spectrum; cutting on pt
   // cuts on |K|, which truncates the two charges at different true momenta
   // (they are offset by 2*Delta) and biases the estimator.
   bool is_prompt = (bxspread1 == 0) && (name == "data_obs" ? is_colliding : true);

   if (is_prompt && nstub1 >= 3) {

      const bool pos = (charge1 > 0);

      h_charge->Fill(charge1, w1);
      p_K_vs_phi->Fill(phi1, HwK1, w1);
      p_K_vs_eta->Fill(eta1, HwK1, w1);

      int isec = (int)((phi1 + 3.14159) / (3.14159/6.0));
      if (isec < 0) isec = 0;  if (isec > 11) isec = 11;
      (pos ? h_sector_plus : h_sector_minus)->Fill(isec, w1);

      if (nstub1 == 3) (pos ? h2_K_processor_plus_stub3 : h2_K_processor_minus_stub3)->Fill(processor1, HwK1, w1);
      if (nstub1 == 4) (pos ? h2_K_processor_plus_stub4 : h2_K_processor_minus_stub4)->Fill(processor1, HwK1, w1);

      (pos ? h2_K_eta_plus : h2_K_eta_minus)->Fill(eta1, HwK1, w1);

      // dxy: the independent probe of the same displacement
      h_dxy->Fill(dxy1, w1);
      p_dxy_vs_phi->Fill(phi1, dxy1, w1);
      (pos ? p_dxy_vs_phi_plus : p_dxy_vs_phi_minus)->Fill(phi1, dxy1, w1);
      p_dxy_vs_eta->Fill(eta1, dxy1, w1);
      h2_dxy_vs_phi->Fill(phi1, dxy1, w1);

      if (dxy1 < 0.5)
         (pos ? h2_K_phi_plus_tightdxy : h2_K_phi_minus_tightdxy)->Fill(phi1, HwK1, w1);
   }

   // ---- MC truth -------------------------------------------------------
   if (name != "data_obs" && genpt1 > 0) {
      double K_gen_hw = (charge1/genpt1) * 2.0 * 1.17 / 0.8569 / (1.25/8192.0);
      double dK_hw    = HwK1 - K_gen_hw;
      h_dK_hw->Fill(dK_hw, w1);
      p_dK_hw_vs_eta->Fill(geneta1, dK_hw, w1);
      p_dK_hw_vs_phi->Fill(phi1, dK_hw, w1);
      p_dK_hw_vs_beta->Fill(genbeta1, dK_hw, w1);
      p_dK_hw_vs_Kgen->Fill(K_gen_hw, dK_hw, w1);
      p_dK_hw_vs_dxy->Fill(dxy1, dK_hw, w1);
      h2_K_hw_reco_gen->Fill(K_gen_hw, HwK1, w1);
      (charge1 > 0 ? p_dK_hw_vs_phi_plus : p_dK_hw_vs_phi_minus)->Fill(phi1, dK_hw, w1);
   }
}

      //########################################################################
      //############################### ANALYSIS ###############################
      //########################################################################
      
      // ############### Across 3 or 4 BXs ################
      
      bool is_wrong = false;
      int cat = -1;
      int cat_wrong = -1;
      int wrong_dir = 0;
      
      bool is_1234_1 = (nstub1==4 and (bxspread1==3210 or bxspread1==3120 or bxspread1==3201 or bxspread1==3102 or bxspread1==3012 or bxspread1==3021 or bxspread1==1230 or bxspread1==1320 or bxspread1==2130 or bxspread1==2310));
      if (name=="fid0p45" or name=="fid0p50" or name=="fid0p55" or name=="fid0p60" or name=="fid0p65" or name=="fid0p70" or name=="fid0p75" or name=="fid0p80") is_1234_1 = false;
      if (is_1234_1){
         if (bxspread1==3210){
            if (qual1>=15) {
               h_stub4_bx1234->Fill(pt1,w1); is_accepted=true;
               cat = 3;
            }
            else if (qual1>=1) {h_stub4_bx1234_fail->Fill(pt1,w1); is_accepted=true;}
         }
         else{
            if (qual1>=15) {
               h_stub4_bx1234_wrong->Fill(pt1,w1); is_accepted=true;
               h_phi_gt2BX_wrong->Fill(phi1,w1);
               if (bxspread1==3120 or bxspread1==3201 or bxspread1==3102 or bxspread1==3012 or bxspread1==3021) {h_stub4_bx1234_wrongU->Fill(pt1,w1); is_wrong = true; cat_wrong = 3; wrong_dir = +1;}
               if (bxspread1==3120 or bxspread1==1230 or bxspread1==1320 or bxspread1==2130 or bxspread1==2310) {h_stub4_bx1234_wrongD->Fill(pt1,w1); is_wrong = true; cat_wrong = 3; wrong_dir = -1;}
            }
            else if (qual1>=1) {
               h_stub4_bx1234_wrong_fail->Fill(pt1,w1); is_accepted=true;
               if (bxspread1==3120 or bxspread1==3201 or bxspread1==3102 or bxspread1==3012 or bxspread1==3021) h_stub4_bx1234_wrongU_fail->Fill(pt1,w1);
               if (bxspread1==3120 or bxspread1==1230 or bxspread1==1320 or bxspread1==2130 or bxspread1==2310) h_stub4_bx1234_wrongD_fail->Fill(pt1,w1);
            }
         // Save individual alternative shapes without quality criteria to compare them together
            if (bxspread1==3120) shape_stub4_bx1234_3120->Fill(pt1,w1);
            if (bxspread1==3201) shape_stub4_bx1234_3201->Fill(pt1,w1);
            if (bxspread1==3102) shape_stub4_bx1234_3102->Fill(pt1,w1);
            if (bxspread1==3012) shape_stub4_bx1234_3012->Fill(pt1,w1);
            if (bxspread1==3021) shape_stub4_bx1234_3021->Fill(pt1,w1);
            if (bxspread1==1230) shape_stub4_bx1234_1230->Fill(pt1,w1);
            if (bxspread1==1320) shape_stub4_bx1234_1320->Fill(pt1,w1);
            if (bxspread1==2130) shape_stub4_bx1234_2130->Fill(pt1,w1);
            if (bxspread1==2310) shape_stub4_bx1234_2310->Fill(pt1,w1);
         }
      }

      bool is_123_3stubs_fast_1 = (nstub1==3 and (bxspread1==2100 or bxspread1==2010 or bxspread1==1200) and (stationspread1==4310 or stationspread1==4210));
      if (name=="fid0p45" or name=="fid0p50" or name=="fid0p55" or name=="fid0p60" or name=="fid0p65" or name=="fid0p70" or name=="fid0p75" or name=="fid0p80") is_123_3stubs_fast_1 = false;
      if (is_123_3stubs_fast_1){
         if (bxspread1==2100){
            if (qual1>=14) {
               h_stub3_bx123_fast->Fill(pt1,w1); is_accepted=true;
               cat = 10;
            }
            else {h_stub3_bx123_fast_fail->Fill(pt1,w1); is_accepted=true;}
         }
         else{
            if (qual1>=14) {
               h_stub3_bx123_fast_wrong->Fill(pt1,w1); is_accepted=true;
               h_phi_gt2BX_wrong->Fill(phi1,w1);
               if (bxspread1==2010) {h_stub3_bx123_fast_wrongU->Fill(pt1,w1); is_wrong = true; cat_wrong = 10; wrong_dir = +1;}
               if (bxspread1==1200) {h_stub3_bx123_fast_wrongD->Fill(pt1,w1); is_wrong = true; cat_wrong = 10; wrong_dir = -1;}
            }
            else {
               h_stub3_bx123_fast_wrong_fail->Fill(pt1,w1); is_accepted=true;
               if (bxspread1==2010) h_stub3_bx123_fast_wrongU_fail->Fill(pt1,w1);
               if (bxspread1==1200) h_stub3_bx123_fast_wrongD_fail->Fill(pt1,w1);
            }
         }
      }

      bool is_123_3stubs_slow_1 = (nstub1==3 and (bxspread1==2100 or bxspread1==2010 or bxspread1==1200) and (stationspread1==4320 or stationspread1==3210));
      if (name=="fid0p45" or name=="fid0p50" or name=="fid0p55" or name=="fid0p60" or name=="fid0p65" or name=="fid0p70" or name=="fid0p75" or name=="fid0p80") is_123_3stubs_slow_1 = false;
         if (is_123_3stubs_slow_1){
            if (bxspread1==2100){
               if (qual1>=14) {
                  h_stub3_bx123_slow->Fill(pt1,w1); is_accepted=true;
                  cat = 9;
               }
               else {h_stub3_bx123_slow_fail->Fill(pt1,w1); is_accepted=true;}
            }
            else{
               if (qual1>=14) {
                  h_stub3_bx123_slow_wrong->Fill(pt1,w1); is_accepted=true;
                  if (bxspread1==2010) {h_stub3_bx123_slow_wrongU->Fill(pt1,w1); is_wrong = true; cat_wrong = 9; wrong_dir = +1;}
                  if (bxspread1==1200){ h_stub3_bx123_slow_wrongD->Fill(pt1,w1); is_wrong = true;cat_wrong = 9; wrong_dir = -1;}
               }
               else {
                  h_stub3_bx123_slow_wrong_fail->Fill(pt1,w1); is_accepted=true;
                  if (bxspread1==2010) h_stub3_bx123_slow_wrongU_fail->Fill(pt1,w1);
                  if (bxspread1==1200) h_stub3_bx123_slow_wrongD_fail->Fill(pt1,w1);
               }
            }
         }

      bool is_123_4stubs_1 = (nstub1==4 and (bxspread1==2100 or bxspread1==2110 or bxspread1==2210 or bxspread1==2120 or bxspread1==2102 or bxspread1==2101 or bxspread1==2201 or bxspread1==2001 or bxspread1==2010 or bxspread1==2012 or bxspread1==2021 or bxspread1==2011 or (bxspread1==120 or bxspread1==210 or bxspread1==1120 or bxspread1==1210 or bxspread1==1220 or bxspread1==1020 or bxspread1==1200)));
      if (name=="fid0p45" or name=="fid0p50" or name=="fid0p55" or name=="fid0p60" or name=="fid0p65" or name=="fid0p70" or name=="fid0p75" or name=="fid0p80") is_123_4stubs_1 = false;
         if (is_123_4stubs_1){
            if (bxspread1==2100 or bxspread1==2110 or bxspread1==2210){
               if (qual1>=15) {
                  h_stub4_bx123->Fill(pt1,w1); is_accepted=true;
                  cat = 2;
               }
               else if (qual1>=1) {h_stub4_bx123_fail->Fill(pt1,w1); is_accepted=true;}
            }
            else{
               if (qual1>=15) {
                  h_stub4_bx123_wrong->Fill(pt1,w1); is_accepted=true;
                  h_phi_gt2BX_wrong->Fill(phi1,w1);
                  if (bxspread1==2120 or bxspread1==2102 or bxspread1==2101 or bxspread1==2201 or bxspread1==2001 or bxspread1==2010 or bxspread1==2012 or bxspread1==2021 or bxspread1==2011){ h_stub4_bx123_wrongU->Fill(pt1,w1); is_wrong = true; cat_wrong = 2; wrong_dir = +1;}
                  if (bxspread1==2120 or bxspread1==2010 or bxspread1==120 or bxspread1==210 or bxspread1==1120 or bxspread1==1210 or bxspread1==1220 or bxspread1==1020 or bxspread1==1200) {h_stub4_bx123_wrongD->Fill(pt1,w1); is_wrong = true; cat_wrong = 2; wrong_dir = -1;}
               }
               else if (qual1>=1) {
                  h_stub4_bx123_wrong_fail->Fill(pt1,w1); is_accepted=true;
                  if (bxspread1==2120 or bxspread1==2102 or bxspread1==2101 or bxspread1==2201 or bxspread1==2001 or bxspread1==2010 or bxspread1==2012 or bxspread1==2021 or bxspread1==2011) h_stub4_bx123_wrongU_fail->Fill(pt1,w1);
                  if (bxspread1==2120 or bxspread1==2010 or bxspread1==120 or bxspread1==210 or bxspread1==1120 or bxspread1==1210 or bxspread1==1220 or bxspread1==1020 or bxspread1==1200) h_stub4_bx123_wrongD_fail->Fill(pt1,w1); 
               }
            }
         }
      
      if (nstub1==4 and (bxspread1==112 or bxspread1==122 or bxspread1==12) and qual1>=15) h_stub4_bx123_neutron->Fill(pt1,w1);
      bool is_124_4stubs_1 = (nstub1==4 and ((bxspread1==3100 or bxspread1==3110 or bxspread1==3310 or bxspread1==3011 or bxspread1==3101 or bxspread1==3001 or bxspread1==3010 or bxspread1==3031 or bxspread1==3013 or bxspread1==3130 or bxspread1==3103 or bxspread1==3031 or bxspread1==3301) or (bxspread1==3200 or bxspread1==3220 or bxspread1==3320 or bxspread1==3022 or bxspread1==3202 or bxspread1==3002 or bxspread1==3020 or bxspread1==3032 or bxspread1==3023 or bxspread1==3230 or bxspread1==3203 or bxspread1==3032 or bxspread1==3302) or (bxspread1==1330 or bxspread1==1130 or bxspread1==1030 or bxspread1==1300 or bxspread1==2330 or bxspread1==2230 or bxspread1==2030 or bxspread1==2300 or bxspread1==130 or bxspread1==310 or bxspread1==230 or bxspread1==320)));
      if (name=="fid0p45" or name=="fid0p50" or name=="fid0p55" or name=="fid0p60" or name=="fid0p65" or name=="fid0p70" or name=="fid0p75" or name=="fid0p80") is_124_4stubs_1 = false;
         if (is_124_4stubs_1){
            if (bxspread1==3100 or bxspread1==3110 or bxspread1==3310 or bxspread1==3200 or bxspread1==3220 or bxspread1==3320){
               if (qual1>=15) {
                  h_stub4_bx124->Fill(pt1,w1); is_accepted=true;
                  cat = 1;
               }
               else if (qual1>=1) {h_stub4_bx124_fail->Fill(pt1,w1); is_accepted=true;}
            }
            else{
               if (qual1>=15) {
                  h_stub4_bx124_wrong->Fill(pt1,w1); is_accepted=true;
                  if ( bxspread1==3011 or bxspread1==3101 or bxspread1==3001 or bxspread1==3010 or bxspread1==3031 or bxspread1==3013 or bxspread1==3130 or bxspread1==3103 or bxspread1==3031 or bxspread1==3301 or bxspread1==3022 or bxspread1==3202 or bxspread1==3002 or bxspread1==3020 or bxspread1==3032 or bxspread1==3023 or bxspread1==3230 or bxspread1==3203 or bxspread1==3032 or bxspread1==3302) {h_stub4_bx124_wrongU->Fill(pt1,w1); is_wrong = true; cat_wrong = 1; wrong_dir = +1;}
                  if (bxspread1==1330 or bxspread1==1130 or bxspread1==1030 or bxspread1==1300 or bxspread1==2330 or bxspread1==2230 or bxspread1==2030 or bxspread1==2300 or bxspread1==130 or bxspread1==310 or bxspread1==230 or bxspread1==320) {h_stub4_bx124_wrongD->Fill(pt1,w1); is_wrong = true; cat_wrong = 1; wrong_dir = -1;}
               }
               else if (qual1>=1) {
                  h_stub4_bx124_wrong_fail->Fill(pt1,w1); is_accepted=true;
                  if ( bxspread1==3011 or bxspread1==3101 or bxspread1==3001 or bxspread1==3010 or bxspread1==3031 or bxspread1==3013 or bxspread1==3130 or bxspread1==3103 or bxspread1==3031 or bxspread1==3301 or bxspread1==3022 or bxspread1==3202 or bxspread1==3002 or bxspread1==3020 or bxspread1==3032 or bxspread1==3023 or bxspread1==3230 or bxspread1==3203 or bxspread1==3032 or bxspread1==3302) h_stub4_bx124_wrongU_fail->Fill(pt1,w1);
                  if (bxspread1==1330 or bxspread1==1130 or bxspread1==1030 or bxspread1==1300 or bxspread1==2330 or bxspread1==2230 or bxspread1==2030 or bxspread1==2300 or bxspread1==130 or bxspread1==310 or bxspread1==230 or bxspread1==320) h_stub4_bx124_wrongD_fail->Fill(pt1,w1);
               }
            }
         }

      bool is_124_3stubs_slow_1 = (nstub1==3 and (bxspread1==3200 or bxspread1==3100 or bxspread1==3020 or bxspread1==3010 or (bxspread1==2300 or bxspread1==1300)) and (stationspread1==4320 or stationspread1==3210));
         if (is_124_3stubs_slow_1){
            if (bxspread1==3200 or bxspread1==3100){
               if (qual1>=14) {
                  h_stub3_bx124_slow->Fill(pt1,w1); is_accepted=true;
                  cat = 7;
               }
               else {h_stub3_bx124_slow_fail->Fill(pt1,w1); is_accepted=true;}
            }
            else{
               if (qual1>=14) {
                  h_stub3_bx124_slow_wrong->Fill(pt1,w1); is_accepted=true;
                  if (bxspread1==3020 or bxspread1==3010) {h_stub3_bx124_slow_wrongU->Fill(pt1,w1); is_wrong = true; cat_wrong = 7; wrong_dir = +1;}
                  if (bxspread1==2300 or bxspread1==1300) {h_stub3_bx124_slow_wrongD->Fill(pt1,w1); is_wrong = true; cat_wrong = 7; wrong_dir = -1;}
               }
               else {
                  h_stub3_bx124_slow_wrong_fail->Fill(pt1,w1); is_accepted=true;
                  if (bxspread1==3020 or bxspread1==3010) h_stub3_bx124_slow_wrongU_fail->Fill(pt1,w1);
                  if (bxspread1==2300 or bxspread1==1300) h_stub3_bx124_slow_wrongD_fail->Fill(pt1,w1);
               }
            }
         }

      bool is_124_3stubs_fast_1 = (nstub1==3 and (bxspread1==3200 or bxspread1==3100 or bxspread1==3020 or bxspread1==3010 or (bxspread1==2300 or bxspread1==1300)) and (stationspread1==4310 or stationspread1==4210));
      if (name=="fid0p45" or name=="fid0p50" or name=="fid0p55" or name=="fid0p60" or name=="fid0p65" or name=="fid0p70" or name=="fid0p75" or name=="fid0p80") is_124_3stubs_fast_1 = false;
         if (is_124_3stubs_fast_1){
            if (bxspread1==3200 or bxspread1==3100){
               if (qual1>=14) {
                  h_stub3_bx124_fast->Fill(pt1,w1); is_accepted=true;
                  cat = 8;
               }
               else {h_stub3_bx124_fast_fail->Fill(pt1,w1); is_accepted=true;}
            }
            else{
               if (qual1>=14) {
                  h_stub3_bx124_fast_wrong->Fill(pt1,w1); is_accepted=true;
                  if ( bxspread1==3020 or bxspread1==3010) {h_stub3_bx124_fast_wrongU->Fill(pt1,w1); is_wrong = true; cat_wrong = 8; wrong_dir = +1;}
                  if (bxspread1==2300 or bxspread1==1300) {h_stub3_bx124_fast_wrongD->Fill(pt1,w1); is_wrong = true; cat_wrong = 8; wrong_dir = -1;}
               }
               else {
                  h_stub3_bx124_fast_wrong_fail->Fill(pt1,w1); is_accepted=true;
                  if ( bxspread1==3020 or bxspread1==3010) h_stub3_bx124_fast_wrongU_fail->Fill(pt1,w1);
                  if (bxspread1==2300 or bxspread1==1300) h_stub3_bx124_fast_wrongD_fail->Fill(pt1,w1);
               }
            }
         }

      //################ Across 2 BX ################

      bool is_1122_4stubs_1 = (nstub1==4 and (bxspread1==1100 or bxspread1==1001 or bxspread1==1010 or (bxspread1==110)));
         if (is_1122_4stubs_1){
            if (bxspread1==1100){
               if (qual1>=15) {
                  h_stub4_bx1122->Fill(pt1,w2); is_accepted=true;
               }
               else if (qual1<15) {h_stub4_bx1122_fail->Fill(pt1,w2); is_accepted=true;}
            }
            else{
               if (qual1>=15) {
                  h_stub4_bx1122_wrong->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==1001 or bxspread1==1010) {h_stub4_bx1122_wrongU->Fill(pt1,w2); is_wrong = true;} //FIXME 1010 used twice
                  if (bxspread1==1010 or bxspread1==110){ h_stub4_bx1122_wrongD->Fill(pt1,w2); is_wrong = true;}
               }
               else if (qual1<15) {
                  h_stub4_bx1122_wrong_fail->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==1001 or bxspread1==1010) h_stub4_bx1122_wrongU_fail->Fill(pt1,w2);
                  if (bxspread1==1010 or bxspread1==110) h_stub4_bx1122_wrongD_fail->Fill(pt1,w2);
               }
            }
         }

      bool is_1122_2tracks_4stubs_1 = (has_2tracks and nstub1==4 and (bxspread1==1100 or bxspread1==1001 or bxspread1==1010 or (bxspread1==110)));
      if (is_1122_2tracks_4stubs_1){
         if (bxspread1==1100){
            if (qual1>=15 and has_2goodtracks) {h_stub4_bx1122_2tracks->Fill(pt1,w2); is_accepted=true; cat = 5;}
            else if (qual1<15) {h_stub4_bx1122_2tracks_fail->Fill(pt1,w2); is_accepted=true;}
         }
         else{
            if (qual1>=15 and has_2goodtracks) {
               h_stub4_bx1122_2tracks_wrong->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==1001 or bxspread1==1010) {h_stub4_bx1122_2tracks_wrongU->Fill(pt1,w2); is_wrong = true; cat_wrong = 5; wrong_dir = +1;}
               if (bxspread1==1010 or bxspread1==110){ h_stub4_bx1122_2tracks_wrongD->Fill(pt1,w2); is_wrong = true; cat_wrong = 5; wrong_dir = -1;}
            }
            else if (qual1<15) {
               h_stub4_bx1122_2tracks_wrong_fail->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==1001 or bxspread1==1010) h_stub4_bx1122_2tracks_wrongU_fail->Fill(pt1,w2);
               if (bxspread1==1010 or bxspread1==110) h_stub4_bx1122_2tracks_wrongD_fail->Fill(pt1,w2);
            }
         }
      }
      ////////

      bool is_1112_4stubs_1 = (nstub1==4 and (bxspread1==1000 or bxspread1==100 or bxspread1==10)); 
      if (is_1112_4stubs_1){
         if (bxspread1==1000){
            if (qual1>=15) {
               h_stub4_bx1112->Fill(pt1,w2); is_accepted=true;
            }
            else if (qual1<15) {h_stub4_bx1112_fail->Fill(pt1,w2); is_accepted=true;}
         }
         else{
            if (qual1>=15) {
               h_stub4_bx1112_wrong->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==100) {h_stub4_bx1112_wrongU->Fill(pt1,w2); is_wrong = true;}//FIXME not same U/D convention
               if (bxspread1==10) {h_stub4_bx1112_wrongD->Fill(pt1,w2); is_wrong = true;}
            }
            else if (qual1<15) {
               h_stub4_bx1112_wrong_fail->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==100) h_stub4_bx1112_wrongU_fail->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==10) h_stub4_bx1112_wrongD_fail->Fill(pt1,w2); is_accepted=true;
            }
         }
      }

      bool is_1112_2tracks_4stubs_1 = (has_2tracks and nstub1==4 and (bxspread1==1000 or bxspread1==100 or bxspread1==10));
      if (is_1112_2tracks_4stubs_1){
         if (bxspread1==1000){
            if (qual1>=15 and has_2goodtracks) { h_stub4_bx1112_2tracks->Fill(pt1,w2); is_accepted=true; cat = 6;}
            else if (qual1<15) {h_stub4_bx1112_2tracks_fail->Fill(pt1,w2); is_accepted=true;}
         }
         else{
            if (qual1>=15 and has_2goodtracks) {
               h_stub4_bx1112_2tracks_wrong->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==100){ h_stub4_bx1112_2tracks_wrongU->Fill(pt1,w2); is_wrong = true; cat_wrong = 6; wrong_dir = +1;}
               if (bxspread1==10){ h_stub4_bx1112_2tracks_wrongD->Fill(pt1,w2); is_wrong = true;cat_wrong = 6; wrong_dir = -1;}
            }
            else if (qual1<15) {
               h_stub4_bx1112_2tracks_wrong_fail->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==100) h_stub4_bx1112_2tracks_wrongU_fail->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==10) h_stub4_bx1112_2tracks_wrongD_fail->Fill(pt1,w2); is_accepted=true;
            }
         }
      }
      //////////

      bool is_1222_4stubs_1 = (nstub1==4 and (bxspread1==1110 or bxspread1==1101 or bxspread1==1011)); 
         if (is_1222_4stubs_1){
            if (bxspread1==1110){
               if (qual1>=15) {
               h_stub4_bx1222->Fill(pt1,w2); is_accepted=true;
            }
               else if (qual1<15) {h_stub4_bx1222_fail->Fill(pt1,w2); is_accepted=true;} 
            }
            else{
               if (qual1>=15) {
               h_stub4_bx1222_wrong->Fill(pt1,w2); is_accepted=true;
         if (bxspread1==1101) {h_stub4_bx1222_wrongU->Fill(pt1,w2); is_wrong = true;}//FIXME not same U/D convention
                  if (bxspread1==1011) {h_stub4_bx1222_wrongD->Fill(pt1,w2); is_wrong = true;}
            }
               else if (qual1<15) { 
               h_stub4_bx1222_wrong_fail->Fill(pt1,w2); is_accepted=true;
         if (bxspread1==1101) h_stub4_bx1222_wrongU_fail->Fill(pt1,w2);
         if (bxspread1==1011) h_stub4_bx1222_wrongD_fail->Fill(pt1,w2);
            } 
            }
         }

      bool is_1222_2tracks_4stubs_1 = (has_2tracks and nstub1==4 and (bxspread1==1110 or bxspread1==1101 or bxspread1==1011));
         if (is_1222_2tracks_4stubs_1){
            if (bxspread1==1110){
               if (qual1>=15 and has_2goodtracks) {h_stub4_bx1222_2tracks->Fill(pt1,w2); is_accepted=true; cat = 4;}
               else if (qual1<15) {h_stub4_bx1222_2tracks_fail->Fill(pt1,w2); is_accepted=true;}
            }
            else{
               if (qual1>=15 and has_2goodtracks) {
                  h_stub4_bx1222_2tracks_wrong->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==1101) {h_stub4_bx1222_2tracks_wrongU->Fill(pt1,w2); is_wrong = true; cat_wrong = 4; wrong_dir = +1;} 
                  if (bxspread1==1011) {h_stub4_bx1222_2tracks_wrongD->Fill(pt1,w2); is_wrong = true; cat_wrong = 4; wrong_dir = -1;}
               }
               else if (qual1<15) {
                  h_stub4_bx1222_2tracks_wrong_fail->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==1101) h_stub4_bx1222_2tracks_wrongU_fail->Fill(pt1,w2);
                  if (bxspread1==1011) h_stub4_bx1222_2tracks_wrongD_fail->Fill(pt1,w2);
               }
            }
         }
      ///////////

      bool is_112_3stubs_fast_1 = (nstub1==3 and (bxspread1==1000 or bxspread1==100) and (stationspread1==4310 or stationspread1==4210));
         if (is_112_3stubs_fast_1){
            if (bxspread1==1000){
               if (qual1>=14) {
                  h_stub3_bx112_fast->Fill(pt1,w2); is_accepted=true;
               }
               else if (qual1<14) {h_stub3_bx112_fast_fail->Fill(pt1,w2); is_accepted=true;}
            }
            else{
               if (qual1>=14) {
                  h_stub3_bx112_fast_wrong->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==100) {h_stub3_bx112_fast_wrongU->Fill(pt1,w2); is_wrong = true;}
                  if (bxspread1==100) {h_stub3_bx112_fast_wrongD->Fill(pt1,w2); is_wrong = true;} //FIXME up and down same
               }
               else if (qual1<14){
                  h_stub3_bx112_fast_wrong_fail->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==100) h_stub3_bx112_fast_wrongU_fail->Fill(pt1,w2);
                  if (bxspread1==100) h_stub3_bx112_fast_wrongD_fail->Fill(pt1,w2);
               }
            }
         }

      bool is_112_3stubs_fast_2tracks_1 = (has_2tracks and nstub1==3 and (bxspread1==1000 or bxspread1==100) and (stationspread1==4310 or stationspread1==4210));
         if (is_112_3stubs_fast_2tracks_1){
            if (bxspread1==1000){
               if (qual1>=14 and has_2goodtracks) { h_stub3_bx112_fast_2tracks->Fill(pt1,w2); is_accepted=true; cat = 14;}
               else if (qual1<14) {h_stub3_bx112_fast_2tracks_fail->Fill(pt1,w2); is_accepted=true;}
            }
            else{
               if (qual1>=14 and has_2goodtracks) {
                  h_stub3_bx112_fast_2tracks_wrong->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==100){ h_stub3_bx112_fast_2tracks_wrongU->Fill(pt1,w2); is_wrong = true; cat_wrong = 14; wrong_dir = +1;}
                  if (bxspread1==100){ h_stub3_bx112_fast_2tracks_wrongD->Fill(pt1,w2);  is_wrong = true; cat_wrong = 14; wrong_dir = -1;}
               }
               else if (qual1<14){
                  h_stub3_bx112_fast_2tracks_wrong_fail->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==100) h_stub3_bx112_fast_2tracks_wrongU_fail->Fill(pt1,w2);
                  if (bxspread1==100) h_stub3_bx112_fast_2tracks_wrongD_fail->Fill(pt1,w2);
               }
            }
         }
      ///////////

      bool is_112_3stubs_slow_1 = (nstub1==3 and (bxspread1==1000 or bxspread1==100) and (stationspread1==3210 or stationspread1==4320));
         if (is_112_3stubs_slow_1){
            if (bxspread1==1000){
               if (qual1>=14) {
               h_stub3_bx112_slow->Fill(pt1,w2); is_accepted=true;
            }
               else if (qual1<14) {h_stub3_bx112_slow_fail->Fill(pt1,w2); is_accepted=true;}
            }
            else{
               if (qual1>=14) {
               h_stub3_bx112_slow_wrong->Fill(pt1,w2); is_accepted=true;
         if (bxspread1==100){ h_stub3_bx112_slow_wrongU->Fill(pt1,w2); is_wrong = true;}
         if (bxspread1==100){ h_stub3_bx112_slow_wrongD->Fill(pt1,w2);  is_wrong = true;}//FIXME up and down same
            }
               else if (qual1<14){
               h_stub3_bx112_slow_wrong_fail->Fill(pt1,w2); is_accepted=true;
         if (bxspread1==100) h_stub3_bx112_slow_wrongU_fail->Fill(pt1,w2);
                  if (bxspread1==100) h_stub3_bx112_slow_wrongD_fail->Fill(pt1,w2);
            }
            }
         }

      bool is_112_3stubs_slow_2tracks_1 = (has_2tracks and nstub1==3 and (bxspread1==1000 or bxspread1==100) and (stationspread1==3210 or stationspread1==4320));
         if (is_112_3stubs_slow_2tracks_1){
            if (bxspread1==1000){
               if (qual1>=14 and has_2goodtracks) { h_stub3_bx112_slow_2tracks->Fill(pt1,w2); is_accepted=true; cat = 13;}
               else if (qual1<14) {h_stub3_bx112_slow_2tracks_fail->Fill(pt1,w2); is_accepted=true;}
            }
            else{
               if (qual1>=14 and has_2goodtracks) {
                  h_stub3_bx112_slow_2tracks_wrong->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==100) {h_stub3_bx112_slow_2tracks_wrongU->Fill(pt1,w2); is_wrong = true; cat_wrong = 13; wrong_dir = +1;}
                  if (bxspread1==100){ h_stub3_bx112_slow_2tracks_wrongD->Fill(pt1,w2);  is_wrong = true; cat_wrong = 13; wrong_dir = -1;}
               }
               else if (qual1<14){
                  h_stub3_bx112_slow_2tracks_wrong_fail->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==100) h_stub3_bx112_slow_2tracks_wrongU_fail->Fill(pt1,w2);
                  if (bxspread1==100) h_stub3_bx112_slow_2tracks_wrongD_fail->Fill(pt1,w2);
               }
            }
         }
      /////////

      bool is_122_3stubs_fast_1 = (nstub1==3 and (bxspread1==1100 or bxspread1==1010) and (stationspread1==4310 or stationspread1==4210));
         if (is_122_3stubs_fast_1){
            if (bxspread1==1100){
               if (qual1>=14) {
                  h_stub3_bx122_fast->Fill(pt1,w2); is_accepted=true;
               }
               else if (qual1<14) {h_stub3_bx122_fast_fail->Fill(pt1,w2); is_accepted=true;}
            }
            else{
               if (qual1>=14) {
                  h_stub3_bx122_fast_wrong->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==1010) {h_stub3_bx122_fast_wrongU->Fill(pt1,w2); is_wrong = true;}//FIXME same up and down
                  if (bxspread1==1010) {h_stub3_bx122_fast_wrongD->Fill(pt1,w2);is_wrong = true;}
               }
               else if (qual1<14){
                  h_stub3_bx122_fast_wrong_fail->Fill(pt1,w2); is_accepted=true;
                  if (bxspread1==1010) h_stub3_bx122_fast_wrongU_fail->Fill(pt1,w2);
                  if (bxspread1==1010) h_stub3_bx122_fast_wrongD_fail->Fill(pt1,w2);
               }
            }
         }

      bool is_122_3stubs_fast_2tracks_1 = (has_2tracks and nstub1==3 and (bxspread1==1100 or bxspread1==1010) and (stationspread1==4310 or stationspread1==4210));
      if (is_122_3stubs_fast_2tracks_1){
         if (bxspread1==1100){
            if (qual1>=14 and has_2goodtracks) {h_stub3_bx122_fast_2tracks->Fill(pt1,w2); is_accepted=true; cat = 12;}
            else if (qual1<14) {h_stub3_bx122_fast_2tracks_fail->Fill(pt1,w2); is_accepted=true;}
         }
         else{
            if (qual1>=14 and has_2goodtracks) {
               h_stub3_bx122_fast_2tracks_wrong->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==1010) {h_stub3_bx122_fast_2tracks_wrongU->Fill(pt1,w2); is_wrong = true; cat_wrong = 12; wrong_dir = +1;}
               if (bxspread1==1010) {h_stub3_bx122_fast_2tracks_wrongD->Fill(pt1,w2);is_wrong = true; cat_wrong = 12; wrong_dir = -1;}
            }
            else if (qual1<14){
               h_stub3_bx122_fast_2tracks_wrong_fail->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==1010) h_stub3_bx122_fast_2tracks_wrongU_fail->Fill(pt1,w2);
               if (bxspread1==1010) h_stub3_bx122_fast_2tracks_wrongD_fail->Fill(pt1,w2);
            }
         }
      }
      ////////////

      bool is_122_3stubs_slow_1 = (nstub1==3 and (bxspread1==1100 or bxspread1==1010) and (stationspread1==3210 or stationspread1==4320));
      if (is_122_3stubs_slow_1){
         if (bxspread1==1100){
            if (qual1>=14) {
               h_stub3_bx122_slow->Fill(pt1,w2); is_accepted=true;
            }
            else if (qual1<14) {h_stub3_bx122_slow_fail->Fill(pt1,w2); is_accepted=true;}
         }
         else{
            if (qual1>=14) {
               h_stub3_bx122_slow_wrong->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==1010) {h_stub3_bx122_slow_wrongU->Fill(pt1,w2); is_wrong = true;}//FIXME same up and down
               if (bxspread1==1010) {h_stub3_bx122_slow_wrongD->Fill(pt1,w2);is_wrong = true;}
            }
            else if (qual1<14){
               h_stub3_bx122_slow_wrong_fail->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==1010) h_stub3_bx122_slow_wrongU_fail->Fill(pt1,w2);
               if (bxspread1==1010) h_stub3_bx122_slow_wrongD_fail->Fill(pt1,w2);
            }
         }
      }

      bool is_122_3stubs_slow_2tracks_1 = (has_2tracks and nstub1==3 and (bxspread1==1100 or bxspread1==1010) and (stationspread1==3210 or stationspread1==4320));
      if (is_122_3stubs_slow_2tracks_1){
         if (bxspread1==1100){
            if (qual1>=14 and has_2goodtracks) {h_stub3_bx122_slow_2tracks->Fill(pt1,w2); is_accepted=true; cat = 11;}
            else if (qual1<14) {h_stub3_bx122_slow_2tracks_fail->Fill(pt1,w2); is_accepted=true;}
         }
         else{
            if (qual1>=14 and has_2goodtracks) {
               h_stub3_bx122_slow_2tracks_wrong->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==1010) {h_stub3_bx122_slow_2tracks_wrongU->Fill(pt1,w2); is_wrong = true; cat_wrong = 11; wrong_dir = +1;}
               if (bxspread1==1010){ h_stub3_bx122_slow_2tracks_wrongD->Fill(pt1,w2);is_wrong = true; cat_wrong = 11; wrong_dir = -1;}
            }
            else if (qual1<14){
               h_stub3_bx122_slow_2tracks_wrong_fail->Fill(pt1,w2); is_accepted=true;
               if (bxspread1==1010) h_stub3_bx122_slow_2tracks_wrongU_fail->Fill(pt1,w2);
               if (bxspread1==1010) h_stub3_bx122_slow_2tracks_wrongD_fail->Fill(pt1,w2);
            }
         }
      }

      if(is_wrong){
         if(nstub1 == 3) h_dxy3stubs_wrong->Fill(dxy1);
         if(nstub1 == 4) h_dxy4stubs_wrong->Fill(dxy1);
         //h_qpt_wrong->Fill(charge1/pt1);
      }

      // if (pt1 > 0) {
      //    double K = charge1/pt1;
      //    if (cat >= 0)                        p_qpt_vs_cat->Fill(cat + 0.5, K, w1);
      //    if (cat_wrong >= 0 && wrong_dir > 0) p_qpt_vs_cat_wrongU->Fill(cat_wrong + 0.5, K, w1);
      //    if (cat_wrong >= 0 && wrong_dir < 0) p_qpt_vs_cat_wrongD->Fill(cat_wrong + 0.5, K, w1);
      // }


   } // end of loop over events

   TFile *fout = TFile::Open(output.c_str(), "RECREATE");
    fout->cd();



   //  // ---- MET ----
   //  h_metbx0_mu0->Write();  h_metbxm1_mu0->Write(); h_metbxm2_mu0->Write();
   //  h_metbxm3_mu0->Write(); h_metbxm4_mu0->Write(); h_metbxm5_mu0->Write();
   //  h_metbx0_mu0fail->Write();  h_metbxm1_mu0fail->Write(); h_metbxm2_mu0fail->Write();
   //  h_metbxm3_mu0fail->Write(); h_metbxm4_mu0fail->Write(); h_metbxm5_mu0fail->Write();

    // ---- basic distributions ----
    h_qual_nstub2->Write();
    h_qual_nstub3->Write();
    h_qual_nstub4->Write();
    h_nstub->Write();
    h_pt->Write();
    h_charge->Write();
    h_dxy->Write();
    h_dxy3stubs->Write();
    h_dxy4stubs->Write();
    h_dxy3stubs_wrong->Write();
    h_dxy4stubs_wrong->Write();
   //  h_phi_2BX->Write();
   //  h_phi_2BX_wrong->Write();
   //  h_phi_gt2BX->Write();
   //  h_phi_gt2BX_wrong->Write();
   //  h_1234ordering->Write();

    // ---- alignment: K marginals (these four are the measurement) ----
    h2_K_processor_plus_stub3->Write();
    h2_K_processor_minus_stub3->Write();
    h2_K_processor_plus_stub4->Write();
    h2_K_processor_minus_stub4->Write();

    h2_K_eta_plus->Write();
    h2_K_eta_minus->Write();

    h_K->Write();
    p_K_vs_phi->Write();
    p_K_vs_eta->Write();
    p_K_vs_nstub->Write();
    p_K_vs_dxy->Write();
    h_sector_plus->Write();
    h_sector_minus->Write();

    // ---- alignment: dxy ----
    h_dxy_prompt->Write();
    p_dxy_vs_phi->Write();
    p_dxy_vs_phi_plus->Write();
    p_dxy_vs_phi_minus->Write();
    p_dxy_vs_eta->Write();
    h2_dxy_vs_phi->Write();
    h2_K_phi_plus_tightdxy->Write();
    h2_K_phi_minus_tightdxy->Write();

    // ---- MC closure (empty for data_obs) ----
    h_dK_hw->Write();
    p_dK_hw_vs_eta->Write();
    p_dK_hw_vs_phi->Write();
    p_dK_hw_vs_beta->Write();
    p_dK_hw_vs_Kgen->Write();
    p_dK_hw_vs_dxy->Write();
    p_dK_hw_vs_phi_plus->Write();
    p_dK_hw_vs_phi_minus->Write();
    h2_K_hw_reco_gen->Write();

    fout->Close();

   //  WriteHistToFile(fout,h_phi_2BX,name,"phi_2BX");
   //  WriteHistToFile(fout,h_phi_2BX_wrong,name,"phi_2BX_wrong");
   //  WriteHistToFile(fout,h_phi_gt2BX,name,"phi_gt2BX");
   //  WriteHistToFile(fout,h_phi_gt2BX_wrong,name,"phi_gt2BX_wrong");

   //  WriteHistToFile(fout,h_stub4_bx1234,name,"stub4_bx1234");
   //  WriteHistToFile(fout,h_stub4_bx1234_wrong,name,"stub4_bx1234_wrong");
   //  WriteHistToFile(fout,h_stub4_bx1234_wrongU,name,"stub4_bx1234_wrongU");
   //  WriteHistToFile(fout,h_stub4_bx1234_wrongD,name,"stub4_bx1234_wrongD");
   //  WriteHistToFile(fout,h_stub4_bx1234_fail,name,"stub4_bx1234_fail");
   //  WriteHistToFile(fout,h_stub4_bx1234_wrong_fail,name,"stub4_bx1234_wrong_fail");
   //  WriteHistToFile(fout,h_stub4_bx1234_wrongU_fail,name,"stub4_bx1234_wrongU_fail");
   //  WriteHistToFile(fout,h_stub4_bx1234_wrongD_fail,name,"stub4_bx1234_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx123_slow,name,"stub3_bx123_slow");
   //  WriteHistToFile(fout,h_stub3_bx123_slow_wrong,name,"stub3_bx123_slow_wrong");
   //  WriteHistToFile(fout,h_stub3_bx123_slow_wrongU,name,"stub3_bx123_slow_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx123_slow_wrongD,name,"stub3_bx123_slow_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx123_slow_fail,name,"stub3_bx123_slow_fail");
   //  WriteHistToFile(fout,h_stub3_bx123_slow_wrong_fail,name,"stub3_bx123_slow_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx123_slow_wrongU_fail,name,"stub3_bx123_slow_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx123_slow_wrongD_fail,name,"stub3_bx123_slow_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx123_fast,name,"stub3_bx123_fast");
   //  WriteHistToFile(fout,h_stub3_bx123_fast_wrong,name,"stub3_bx123_fast_wrong");
   //  WriteHistToFile(fout,h_stub3_bx123_fast_wrongU,name,"stub3_bx123_fast_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx123_fast_wrongD,name,"stub3_bx123_fast_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx123_fast_fail,name,"stub3_bx123_fast_fail");
   //  WriteHistToFile(fout,h_stub3_bx123_fast_wrong_fail,name,"stub3_bx123_fast_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx123_fast_wrongU_fail,name,"stub3_bx123_fast_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx123_fast_wrongD_fail,name,"stub3_bx123_fast_wrongD_fail");

   //  WriteHistToFile(fout,h_stub4_bx123,name,"stub4_bx123");
   //  WriteHistToFile(fout,h_stub4_bx123_neutron,name,"stub4_bx123_neutron");
   //  WriteHistToFile(fout,h_stub4_bx123_wrong,name,"stub4_bx123_wrong");
   //  WriteHistToFile(fout,h_stub4_bx123_wrongU,name,"stub4_bx123_wrongU");
   //  WriteHistToFile(fout,h_stub4_bx123_wrongD,name,"stub4_bx123_wrongD");
   //  WriteHistToFile(fout,h_stub4_bx123_fail,name,"stub4_bx123_fail");
   //  WriteHistToFile(fout,h_stub4_bx123_wrong_fail,name,"stub4_bx123_wrong_fail");
   //  WriteHistToFile(fout,h_stub4_bx123_wrongU_fail,name,"stub4_bx123_wrongU_fail");
   //  WriteHistToFile(fout,h_stub4_bx123_wrongD_fail,name,"stub4_bx123_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx124_slow,name,"stub3_bx124_slow");
   //  WriteHistToFile(fout,h_stub3_bx124_slow_wrong,name,"stub3_bx124_slow_wrong");
   //  WriteHistToFile(fout,h_stub3_bx124_slow_wrongU,name,"stub3_bx124_slow_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx124_slow_wrongD,name,"stub3_bx124_slow_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx124_slow_fail,name,"stub3_bx124_slow_fail");
   //  WriteHistToFile(fout,h_stub3_bx124_slow_wrong_fail,name,"stub3_bx124_slow_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx124_slow_wrongU_fail,name,"stub3_bx124_slow_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx124_slow_wrongD_fail,name,"stub3_bx124_slow_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx124_fast,name,"stub3_bx124_fast");
   //  WriteHistToFile(fout,h_stub3_bx124_fast_wrong,name,"stub3_bx124_fast_wrong");
   //  WriteHistToFile(fout,h_stub3_bx124_fast_wrongU,name,"stub3_bx124_fast_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx124_fast_wrongD,name,"stub3_bx124_fast_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx124_fast_fail,name,"stub3_bx124_fast_fail");
   //  WriteHistToFile(fout,h_stub3_bx124_fast_wrong_fail,name,"stub3_bx124_fast_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx124_fast_wrongU_fail,name,"stub3_bx124_fast_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx124_fast_wrongD_fail,name,"stub3_bx124_fast_wrongD_fail");

   //  WriteHistToFile(fout,h_stub4_bx124,name,"stub4_bx124");
   //  WriteHistToFile(fout,h_stub4_bx124_wrong,name,"stub4_bx124_wrong");
   //  WriteHistToFile(fout,h_stub4_bx124_wrongU,name,"stub4_bx124_wrongU");
   //  WriteHistToFile(fout,h_stub4_bx124_wrongD,name,"stub4_bx124_wrongD");
   //  WriteHistToFile(fout,h_stub4_bx124_fail,name,"stub4_bx124_fail");
   //  WriteHistToFile(fout,h_stub4_bx124_wrong_fail,name,"stub4_bx124_wrong_fail");
   //  WriteHistToFile(fout,h_stub4_bx124_wrongU_fail,name,"stub4_bx124_wrongU_fail");
   //  WriteHistToFile(fout,h_stub4_bx124_wrongD_fail,name,"stub4_bx124_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx112_fast,name,"stub3_bx112_fast");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_wrong,name,"stub3_bx112_fast_wrong");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_wrongU,name,"stub3_bx112_fast_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_wrongD,name,"stub3_bx112_fast_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_fail,name,"stub3_bx112_fast_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_wrong_fail,name,"stub3_bx112_fast_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_wrongU_fail,name,"stub3_bx112_fast_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_wrongD_fail,name,"stub3_bx112_fast_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx112_fast_2tracks,name,"stub3_bx112_fast_2tracks");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_2tracks_wrong,name,"stub3_bx112_fast_2tracks_wrong");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_2tracks_wrongU,name,"stub3_bx112_fast_2tracks_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_2tracks_wrongD,name,"stub3_bx112_fast_2tracks_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_2tracks_fail,name,"stub3_bx112_fast_2tracks_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_2tracks_wrong_fail,name,"stub3_bx112_fast_2tracks_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_2tracks_wrongU_fail,name,"stub3_bx112_fast_2tracks_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_fast_2tracks_wrongD_fail,name,"stub3_bx112_fast_2tracks_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx112_slow,name,"stub3_bx112_slow");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_wrong,name,"stub3_bx112_slow_wrong");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_wrongU,name,"stub3_bx112_slow_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_wrongD,name,"stub3_bx112_slow_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_fail,name,"stub3_bx112_slow_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_wrong_fail,name,"stub3_bx112_slow_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_wrongU_fail,name,"stub3_bx112_slow_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_wrongD_fail,name,"stub3_bx112_slow_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx112_slow_2tracks,name,"stub3_bx112_slow_2tracks");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_2tracks_wrong,name,"stub3_bx112_slow_2tracks_wrong");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_2tracks_wrongU,name,"stub3_bx112_slow_2tracks_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_2tracks_wrongD,name,"stub3_bx112_slow_2tracks_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_2tracks_fail,name,"stub3_bx112_slow_2tracks_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_2tracks_wrong_fail,name,"stub3_bx112_slow_2tracks_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_2tracks_wrongU_fail,name,"stub3_bx112_slow_2tracks_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx112_slow_2tracks_wrongD_fail,name,"stub3_bx112_slow_2tracks_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx122_fast,name,"stub3_bx122_fast");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_wrong,name,"stub3_bx122_fast_wrong");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_wrongU,name,"stub3_bx122_fast_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_wrongD,name,"stub3_bx122_fast_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_fail,name,"stub3_bx122_fast_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_wrong_fail,name,"stub3_bx122_fast_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_wrongU_fail,name,"stub3_bx122_fast_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_wrongD_fail,name,"stub3_bx122_fast_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx122_fast_2tracks,name,"stub3_bx122_fast_2tracks");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_2tracks_wrong,name,"stub3_bx122_fast_2tracks_wrong");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_2tracks_wrongU,name,"stub3_bx122_fast_2tracks_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_2tracks_wrongD,name,"stub3_bx122_fast_2tracks_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_2tracks_fail,name,"stub3_bx122_fast_2tracks_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_2tracks_wrong_fail,name,"stub3_bx122_fast_2tracks_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_2tracks_wrongU_fail,name,"stub3_bx122_fast_2tracks_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_fast_2tracks_wrongD_fail,name,"stub3_bx122_fast_2tracks_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx122_slow,name,"stub3_bx122_slow");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_wrong,name,"stub3_bx122_slow_wrong");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_wrongU,name,"stub3_bx122_slow_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_wrongD,name,"stub3_bx122_slow_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_fail,name,"stub3_bx122_slow_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_wrong_fail,name,"stub3_bx122_slow_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_wrongU_fail,name,"stub3_bx122_slow_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_wrongD_fail,name,"stub3_bx122_slow_wrongD_fail");

   //  WriteHistToFile(fout,h_stub3_bx122_slow_2tracks,name,"stub3_bx122_slow_2tracks");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_2tracks_wrong,name,"stub3_bx122_slow_2tracks_wrong");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_2tracks_wrongU,name,"stub3_bx122_slow_2tracks_wrongU");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_2tracks_wrongD,name,"stub3_bx122_slow_2tracks_wrongD");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_2tracks_fail,name,"stub3_bx122_slow_2tracks_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_2tracks_wrong_fail,name,"stub3_bx122_slow_2tracks_wrong_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_2tracks_wrongU_fail,name,"stub3_bx122_slow_2tracks_wrongU_fail");
   //  WriteHistToFile(fout,h_stub3_bx122_slow_2tracks_wrongD_fail,name,"stub3_bx122_slow_2tracks_wrongD_fail");

   //  WriteHistToFile(fout,h_stub4_bx1122,name,"stub4_bx1122");
   //  WriteHistToFile(fout,h_stub4_bx1122_wrong,name,"stub4_bx1122_wrong");
   //  WriteHistToFile(fout,h_stub4_bx1122_wrongU,name,"stub4_bx1122_wrongU");
   //  WriteHistToFile(fout,h_stub4_bx1122_wrongD,name,"stub4_bx1122_wrongD");
   //  WriteHistToFile(fout,h_stub4_bx1122_fail,name,"stub4_bx1122_fail");
   //  WriteHistToFile(fout,h_stub4_bx1122_wrong_fail,name,"stub4_bx1122_wrong_fail");
   //  WriteHistToFile(fout,h_stub4_bx1122_wrongU_fail,name,"stub4_bx1122_wrongU_fail");
   //  WriteHistToFile(fout,h_stub4_bx1122_wrongD_fail,name,"stub4_bx1122_wrongD_fail");

   //  WriteHistToFile(fout,h_stub4_bx1122_2tracks,name,"stub4_bx1122_2tracks");
   //  WriteHistToFile(fout,h_stub4_bx1122_2tracks_wrong,name,"stub4_bx1122_2tracks_wrong");
   //  WriteHistToFile(fout,h_stub4_bx1122_2tracks_wrongU,name,"stub4_bx1122_2tracks_wrongU");
   //  WriteHistToFile(fout,h_stub4_bx1122_2tracks_wrongD,name,"stub4_bx1122_2tracks_wrongD");
   //  WriteHistToFile(fout,h_stub4_bx1122_2tracks_fail,name,"stub4_bx1122_2tracks_fail");
   //  WriteHistToFile(fout,h_stub4_bx1122_2tracks_wrong_fail,name,"stub4_bx1122_2tracks_wrong_fail");
   //  WriteHistToFile(fout,h_stub4_bx1122_2tracks_wrongU_fail,name,"stub4_bx1122_2tracks_wrongU_fail");
   //  WriteHistToFile(fout,h_stub4_bx1122_2tracks_wrongD_fail,name,"stub4_bx1122_2tracks_wrongD_fail");

   //  WriteHistToFile(fout,h_stub4_bx1112,name,"stub4_bx1112");
   //  WriteHistToFile(fout,h_stub4_bx1112_wrong,name,"stub4_bx1112_wrong");
   //  WriteHistToFile(fout,h_stub4_bx1112_wrongU,name,"stub4_bx1112_wrongU");
   //  WriteHistToFile(fout,h_stub4_bx1112_wrongD,name,"stub4_bx1112_wrongD");
   //  WriteHistToFile(fout,h_stub4_bx1112_fail,name,"stub4_bx1112_fail");
   //  WriteHistToFile(fout,h_stub4_bx1112_wrong_fail,name,"stub4_bx1112_wrong_fail");
   //  WriteHistToFile(fout,h_stub4_bx1112_wrongU_fail,name,"stub4_bx1112_wrongU_fail");
   //  WriteHistToFile(fout,h_stub4_bx1112_wrongD_fail,name,"stub4_bx1112_wrongD_fail");

   //  WriteHistToFile(fout,h_stub4_bx1112_2tracks,name,"stub4_bx1112_2tracks");
   //  WriteHistToFile(fout,h_stub4_bx1112_2tracks_wrong,name,"stub4_bx1112_2tracks_wrong");
   //  WriteHistToFile(fout,h_stub4_bx1112_2tracks_wrongU,name,"stub4_bx1112_2tracks_wrongU");
   //  WriteHistToFile(fout,h_stub4_bx1112_2tracks_wrongD,name,"stub4_bx1112_2tracks_wrongD");
   //  WriteHistToFile(fout,h_stub4_bx1112_2tracks_fail,name,"stub4_bx1112_2tracks_fail");
   //  WriteHistToFile(fout,h_stub4_bx1112_2tracks_wrong_fail,name,"stub4_bx1112_2tracks_wrong_fail");
   //  WriteHistToFile(fout,h_stub4_bx1112_2tracks_wrongU_fail,name,"stub4_bx1112_2tracks_wrongU_fail");
   //  WriteHistToFile(fout,h_stub4_bx1112_2tracks_wrongD_fail,name,"stub4_bx1112_2tracks_wrongD_fail");

   //  WriteHistToFile(fout,h_stub4_bx1222,name,"stub4_bx1222");
   //  WriteHistToFile(fout,h_stub4_bx1222_wrong,name,"stub4_bx1222_wrong");
   //  WriteHistToFile(fout,h_stub4_bx1222_wrongU,name,"stub4_bx1222_wrongU");
   //  WriteHistToFile(fout,h_stub4_bx1222_wrongD,name,"stub4_bx1222_wrongD");
   //  WriteHistToFile(fout,h_stub4_bx1222_fail,name,"stub4_bx1222_fail");
   //  WriteHistToFile(fout,h_stub4_bx1222_wrong_fail,name,"stub4_bx1222_wrong_fail");
   //  WriteHistToFile(fout,h_stub4_bx1222_wrongU_fail,name,"stub4_bx1222_wrongU_fail");
   //  WriteHistToFile(fout,h_stub4_bx1222_wrongD_fail,name,"stub4_bx1222_wrongD_fail");

   //  WriteHistToFile(fout,h_stub4_bx1222_2tracks,name,"stub4_bx1222_2tracks");
   //  WriteHistToFile(fout,h_stub4_bx1222_2tracks_wrong,name,"stub4_bx1222_2tracks_wrong");
   //  WriteHistToFile(fout,h_stub4_bx1222_2tracks_wrongU,name,"stub4_bx1222_2tracks_wrongU");
   //  WriteHistToFile(fout,h_stub4_bx1222_2tracks_wrongD,name,"stub4_bx1222_2tracks_wrongD");
   //  WriteHistToFile(fout,h_stub4_bx1222_2tracks_fail,name,"stub4_bx1222_2tracks_fail");
   //  WriteHistToFile(fout,h_stub4_bx1222_2tracks_wrong_fail,name,"stub4_bx1222_2tracks_wrong_fail");
   //  WriteHistToFile(fout,h_stub4_bx1222_2tracks_wrongU_fail,name,"stub4_bx1222_2tracks_wrongU_fail");
   //  WriteHistToFile(fout,h_stub4_bx1222_2tracks_wrongD_fail,name,"stub4_bx1222_2tracks_wrongD_fail");

    //cout<<h_stub4_bx1222_2tracks->Integral()<<endl;

    fout->Close();
}