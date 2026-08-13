// build:
//   g++ -O3 TreeClassification.cc -o TreeClassification.exe $(root-config --cflags --glibs)
// run:
//   ./TreeClassification.exe /eos/user/f/flarover/DATA/AnalysisHSCP/modKBMTF/DupRemove/395343_394_0/output_1-1.root /eos/user/f/flarover/DATA/AnalysisHSCP/modKBMTF/test.root data_obs
//
// Unbinned counterpart of ModROC_SlowAnalysis_classification.cc.
// Only SIGNAL-LIKE tracks are written out, i.e. exactly the tracks that fill
// the nominal h_<category> histograms of the ROC script:
//   * data: earlier BX must be colliding (non-colliding bunches are dropped)
//   * only the "right" (correctly ordered) BX patterns -- no _wrong sideband
//   * quality above the per-category threshold  -- no _fail
//   * the pT threshold of the category (see table below)

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TNamed.h>
#include <TRandom3.h>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <TChain.h>
#include <fstream>

#include "tr_Tree.h"
#include "myHelper.h"
#include "fiducial_weight.h"

using namespace std;

// ---------------------------------------------------------------------------
//  Category table
// ---------------------------------------------------------------------------
//  right     : BX patterns treated as correctly-ordered (signal-like).
//              The _wrong sidebands are deliberately NOT listed: they are not
//              written to the tree at all.
//  stations  : allowed stationspread values; empty = no requirement
//  qualMin   : quality threshold of the category (qual1 >= qualMin);
//              tracks below it go to the _fail histograms and are dropped here
//  ptMin/Max : in the ROC script the pT selection is implicit in the histogram
//              range -- entries below the first bin edge land in the underflow
//              and never enter the datacard. Those lower edges are the real
//              analysis thresholds and ARE applied here (applyPtMin).
//              ptMax is the last bin edge; the ROC script never adds the
//              overflow either, but that is a binning artefact rather than a
//              physics cut, so it is off by default (applyPtMax).
//  need2gt   : the *_2tracks signal regions, which require a second good track
//              and in exchange use a looser pT threshold (50 instead of 150).
//
//  NOTE: as in the ROC script, a 2-BX track that has a second good track fills
//  BOTH the single-track and the _2tracks region, so it produces two entries
//  here (same track, different cat). Select on `cat` when counting.

struct Category {
   string        name;
   int           nstub;
   set<int>      right;
   set<int>      stations;   // empty => any
   int           qualMin;
   float         ptMin;
   float         ptMax;
   bool          need2gt;
};

static const bool applyPtMin = true;
static const bool applyPtMax = false;   // true => byte-identical to the histogram yields

static const vector<Category> CATS = {

  // ---------------- 4 stubs, >2 BX ----------------
  { "stub4_bx1234", 4, {3210},                               {},           15,  15.f,  50.f, false },
  { "stub4_bx123",  4, {2100,2110,2210},                     {},           15, 100.f, 500.f, false },
  { "stub4_bx124",  4, {3100,3110,3310,3200,3220,3320},      {},           15,  50.f, 500.f, false },

  // ---------------- 3 stubs, >2 BX ----------------
  { "stub3_bx123_slow", 3, {2100},      {4320,3210}, 14, 100.f, 500.f, false },
  { "stub3_bx123_fast", 3, {2100},      {4310,4210}, 14,  50.f, 500.f, false },
  { "stub3_bx124_slow", 3, {3200,3100}, {4320,3210}, 14, 100.f, 500.f, false },
  { "stub3_bx124_fast", 3, {3200,3100}, {4310,4210}, 14,  50.f, 500.f, false },

  // ---------------- 4 stubs, 2 BX ----------------
  { "stub4_bx1122", 4, {1100}, {}, 15, 150.f, 750.f, false },
  { "stub4_bx1112", 4, {1000}, {}, 15, 150.f, 750.f, false },
  { "stub4_bx1222", 4, {1110}, {}, 15, 150.f, 750.f, false },

  // ---------------- 3 stubs, 2 BX ----------------
  { "stub3_bx112_slow", 3, {1000}, {3210,4320}, 14, 150.f, 750.f, false },
  { "stub3_bx112_fast", 3, {1000}, {4310,4210}, 14, 150.f, 750.f, false },
  { "stub3_bx122_slow", 3, {1100}, {3210,4320}, 14, 150.f, 750.f, false },
  { "stub3_bx122_fast", 3, {1100}, {4310,4210}, 14, 150.f, 750.f, false },

  // ---------------- 2 BX, second good track ----------------
  { "stub4_bx1122_2tracks", 4, {1100}, {},          15, 50.f, 500.f, true },
  { "stub4_bx1112_2tracks", 4, {1000}, {},          15, 50.f, 500.f, true },
  { "stub4_bx1222_2tracks", 4, {1110}, {},          15, 50.f, 500.f, true },

  { "stub3_bx112_slow_2tracks", 3, {1000}, {3210,4320}, 14, 50.f, 500.f, true },
  { "stub3_bx112_fast_2tracks", 3, {1000}, {4310,4210}, 14, 50.f, 500.f, true },
  { "stub3_bx122_slow_2tracks", 3, {1100}, {3210,4320}, 14, 50.f, 500.f, true },
  { "stub3_bx122_fast_2tracks", 3, {1100}, {4310,4210}, 14, 50.f, 500.f, true },
};


int main(int argc, char** argv) {

   if (argc < 4) { cerr << "usage: exe <in> <out> <name>\n"; return 1; }
   const string input  = argv[1];
   const string output = argv[2];
   const string name   = argv[3];
   const bool   isData = (name == "data_obs");

   TFile *fin = TFile::Open(input.c_str(), "READ");
   TTree *arbre = (TTree*) fin->Get("Events");

   // ---------------- input branches ----------------
   arbre->SetBranchAddress("bunchCrossing", &bunchCrossing);
   arbre->SetBranchAddress("nL1KBMTFSkimmed", &nL1KBMTFSkimmed);
   if (!isData) arbre->SetBranchAddress("ngen", &ngen);
   if (isData) {
      arbre->SetBranchAddress("is_colliding", &is_colliding);
      arbre->SetBranchAddress("is_earlier_colliding", &is_earlier_colliding);
   }
   arbre->SetBranchAddress("nstub1", &nstub1);
   arbre->SetBranchAddress("bxspread1", &bxspread1);
   arbre->SetBranchAddress("stationspread1", &stationspread1);
   arbre->SetBranchAddress("firstbx1", &firstbx1);
   arbre->SetBranchAddress("qual1", &qual1);
   arbre->SetBranchAddress("pt1", &pt1);
   arbre->SetBranchAddress("eta1", &eta1);
   arbre->SetBranchAddress("phi1", &phi1);
   arbre->SetBranchAddress("dxy1", &dxy1);
   arbre->SetBranchAddress("charge1", &charge1);
   arbre->SetBranchAddress("recobeta1", &recobeta1);
   arbre->SetBranchAddress("genpt1", &genpt1);
   arbre->SetBranchAddress("geneta1", &geneta1);
   arbre->SetBranchAddress("genbeta1", &genbeta1);
   // second track -- nstub2 was missing before, which silently left it at 0 and
   // made has_2goodtracks always false.
   arbre->SetBranchAddress("idx2", &idx2);
   arbre->SetBranchAddress("nstub2", &nstub2);
   arbre->SetBranchAddress("pt2", &pt2);
   arbre->SetBranchAddress("dxy2", &dxy2);
   arbre->SetBranchAddress("qual2", &qual2);

   // ---------------- output ----------------
   TFile *fout = TFile::Open(output.c_str(), "RECREATE");

   Int_t    o_cat, o_nstub, o_qual, o_bxspread, o_stationspread, o_firstbx, o_charge;
   UInt_t   o_bx;
   Double_t o_pt, o_eta, o_phi, o_dxy, o_recobeta;
   Double_t o_genpt, o_geneta, o_genbeta;
   Double_t o_weight;
   Int_t    o_nstub2, o_qual2;
   Double_t o_pt2, o_dxy2;
   Bool_t   o_has_2tracks, o_has_2goodtracks;

   TTree *tout = new TTree("Tracks", "signal-like classified tracks");

   // NOTE: leaflist letter MUST match the C++ type. Declaring /F while writing
   // a Double_t silently corrupts the following branch -- the exact bug that
   // hit dxy1/dxy2 before. Everything above is Double_t => /D.
   tout->Branch("cat",            &o_cat,            "cat/I");
   tout->Branch("bunchCrossing",  &o_bx,             "bunchCrossing/i");
   tout->Branch("nstub",          &o_nstub,          "nstub/I");
   tout->Branch("qual",           &o_qual,           "qual/I");
   tout->Branch("bxspread",       &o_bxspread,       "bxspread/I");
   tout->Branch("stationspread",  &o_stationspread,  "stationspread/I");
   tout->Branch("firstbx",        &o_firstbx,        "firstbx/I");
   tout->Branch("charge",         &o_charge,         "charge/I");

   tout->Branch("pt",             &o_pt,             "pt/D");
   tout->Branch("eta",            &o_eta,            "eta/D");
   tout->Branch("phi",            &o_phi,            "phi/D");
   tout->Branch("dxy",            &o_dxy,            "dxy/D");
   tout->Branch("recobeta",       &o_recobeta,       "recobeta/D");

   tout->Branch("genpt",          &o_genpt,          "genpt/D");
   tout->Branch("geneta",         &o_geneta,         "geneta/D");
   tout->Branch("genbeta",        &o_genbeta,        "genbeta/D");

   tout->Branch("weight",         &o_weight,         "weight/D");

   tout->Branch("nstub2",         &o_nstub2,         "nstub2/I");
   tout->Branch("qual2",          &o_qual2,          "qual2/I");
   tout->Branch("pt2",            &o_pt2,            "pt2/D");
   tout->Branch("dxy2",           &o_dxy2,           "dxy2/D");
   tout->Branch("has_2tracks",    &o_has_2tracks,    "has_2tracks/O");
   tout->Branch("has_2goodtracks",&o_has_2goodtracks,"has_2goodtracks/O");

   // ---------------- normalization ----------------
   // same xsweight as the ROC script: (lumi * fiducial) / ngen, 1 for data
   const float lumiweight = 5380.0;
   float fidweight = Get_fidweight(name, input);
   if (name.find("fid0p") != std::string::npos) fidweight = 0.5*fidweight; // sum 2 samples for better stat
   static TRandom3 randGen(1234);

   Long64_t nentries = arbre->GetEntries();

   //  event loop
   for (Long64_t i = 0; i < nentries; ++i) {
      arbre->GetEntry(i);
      if (i % 100000 == 0) { fprintf(stdout,"\r  %lld / %lld ", i, nentries); fflush(stdout); }

      // ---- non-colliding bunches are not signal-like: drop them (data only) ----
      if (isData && !is_earlier_colliding) continue;

      // ---- muon momentum smearing (MC only) ----
      if (!isData) {
         pt1 *= randGen.Gaus(1.0, 0.1148);
         pt2 *= randGen.Gaus(1.0, 0.1148);
      }

      const double xsweight = (isData || ngen <= 0) ? 1.0 : (lumiweight*fidweight)/ngen;

      // ---- second-track flags ----
      const bool has_2tracks = (idx2 < 90);
      bool has_2goodtracks = false;
      if (idx2 < 90 && pt2 > 100 && dxy2 < 1 && qual2 > 12) {
         if (nstub2 == 2 && qual2 > 12) has_2goodtracks = true;
         if (nstub2 == 3 && qual2 > 13) has_2goodtracks = true;
         if (nstub2 == 4 && qual2 > 14) has_2goodtracks = true;
      }

      // ---- classify against the table ----
      // No `break`: the 2-BX regions come in a single-track and a _2tracks
      // flavour, and the ROC script fills both for the same track.
      for (size_t c = 0; c < CATS.size(); ++c) {
         const Category &C = CATS[c];

         if (nstub1 != C.nstub) continue;
         if (!C.right.count(bxspread1)) continue;                        // signal-like BX ordering only
         if (!C.stations.empty() && !C.stations.count(stationspread1)) continue;
         if (qual1 < C.qualMin) continue;                                // drop the _fail categories
         if (C.need2gt && !has_2goodtracks) continue;
         if (applyPtMin && pt1 < C.ptMin) continue;
         if (applyPtMax && pt1 >= C.ptMax) continue;

         o_cat = (Int_t) c;
         o_bx  = bunchCrossing;

         o_nstub = nstub1; o_qual = qual1; o_charge = charge1;
         o_bxspread = bxspread1; o_stationspread = stationspread1;
         o_firstbx = firstbx1;

         o_pt = pt1; o_eta = eta1; o_phi = phi1; o_dxy = dxy1;
         o_recobeta = recobeta1;

         o_genpt   = isData ? 0. : genpt1;
         o_geneta  = isData ? 0. : geneta1;
         o_genbeta = isData ? 0. : genbeta1;

         o_weight  = xsweight;

         o_nstub2 = nstub2; o_qual2 = qual2; o_pt2 = pt2; o_dxy2 = dxy2;
         o_has_2tracks     = has_2tracks;
         o_has_2goodtracks = has_2goodtracks;

         tout->Fill();
      }
   }
   cout << endl;

   // category id -> name, stored alongside the tree
  //  string joined;
  //  for (size_t c = 0; c < CATS.size(); ++c) {
  //     if (c) joined += ",";
  //     joined += CATS[c].name;
  //  }
  //  TNamed catmap("catmap", joined.c_str());

   fout->cd();
   tout->Write();
   //catmap.Write();
   fout->Close();
   fin->Close();
   return 0;
   }