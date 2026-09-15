// =============================================================================
//  MisalignmentFactor.cc   --   DATA RETRIEVAL ONLY
//
//  Fills and writes the histograms needed to measure the charge-odd curvature
//  offset of the KBMTF. It deliberately does NO analysis: no percentiles, no
//  fits, no derived numbers, no tables. All of that lives in
//  MisalignmentFactor.ipynb, which reads this file's output.
//
//    compile:  g++ -O3 MisalignmentFactor.cc -o MisalignmentFactor.exe $(root-config --cflags --glibs)
//    run:      ./MisalignmentFactor.exe <in.root> <out.root>
//
//  Runs on the Zmumu data skims and on the Zmumu/DY simulation: branch names
//  match and the scalar types are read off the file rather than assumed.


#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TNamed.h>
#include <TString.h>

#include <cstdio>
#include <cmath>

// ---------------------------------------------------------------- ptLUT scale
static const double LSB = 1.25 / 8192.0;

// |K| -> pT with the misalignment term left out; it is what we are measuring.
// Used only for the pT cut and the pT histograms, never for the measurement.
static double PtFromK(double absK) {
   double k = absK * LSB;
   double F = 0.8569 * k / (1.0 + 0.1144 * k);
   return (F > 1e-12) ? 1.17 / F : 1e9;
}

// ---------------------------------------------------------------- binning
static const int    NSEC  = 12;              // BMTF sectors, 30 deg each
static const int    NETA  = 6;               // signed: a twist is odd in eta
static const double ETAMAX = 0.83;

static const int    NK    = 2400;            // 0.25 hw units
static const double KMAX  = 600.0;
static const double PTMIN = 15.0;

// sector centred on phi = s*30deg + rotation
static int sectorOf(double phi, double rotation) {
   double w = TMath::Pi() / 6.0;
   int s = (int) floor((phi - rotation + w / 2.0) / w);
   while (s < 0)     s += NSEC;
   while (s >= NSEC) s -= NSEC;
   return s;
}

static int etaBinOf(double eta) {
   int b = (int) ((eta + ETAMAX) / (2 * ETAMAX) * NETA);
   return (b < 0) ? 0 : (b >= NETA ? NETA - 1 : b);
}

// ---------------------------------------------------------------- branch read
// charge / qual / nstub / dxy / bxspread are not the same scalar type in every
// ntuple version (int16 in the data skims, int32 or double elsewhere). Reading
// through the leaf adapts to whatever is on file and cannot return garbage.
struct Br {
   TLeaf* leaf = 0;
   void set(TTree* t, const char* name, const char* alt = 0) {
      leaf = t->GetLeaf(name);
      if (!leaf && alt) leaf = t->GetLeaf(alt);
   }
   bool ok() const { return leaf != 0; }
   double operator()() const { return leaf ? leaf->GetValue(0) : 0.0; }
};

// one muon slot; the ntuples carry two, <name>1 and <name>2
struct Muon {
   Br K, eta, phi, charge, qual, nstub, dxy, bxspread;
   void bind(TTree* t, int i) {
      char b[32], c[32];
      sprintf(b, "HwK%d", i); sprintf(c, "hwK%d", i); K.set(t, b, c);
      sprintf(b, "eta%d", i);      eta.set(t, b);
      sprintf(b, "phi%d", i);      phi.set(t, b);
      sprintf(b, "charge%d", i);   charge.set(t, b);
      sprintf(b, "qual%d", i);     qual.set(t, b);
      sprintf(b, "nstub%d", i);    nstub.set(t, b);
      sprintf(b, "dxy%d", i);      dxy.set(t, b);
      sprintf(b, "bxspread%d", i); bxspread.set(t, b);
   }
   bool present() const { return K.ok() && eta.ok() && phi.ok() && charge.ok(); }
};

// ---------------------------------------------------------------- booking
static TH1D* bookK(const char* name, const char* title) {
   TH1D* h = new TH1D(name, Form("%s;|K| [hw units];tracks", title), NK, 0, KMAX);
   h->Sumw2();
   return h;
}
static TH1D* book1(const char* name, const char* title, int n, double lo, double hi) {
   TH1D* h = new TH1D(name, title, n, lo, hi);
   h->Sumw2();
   return h;
}

// =============================================================================
int main(int argc, char** argv) {

   if (argc < 3) { printf("usage: %s <in.root> <out.root>\n", argv[0]); return 1; }

   TFile* fin = TFile::Open(argv[1]);
   if (!fin || fin->IsZombie()) { printf("cannot open %s\n", argv[1]); return 1; }
   TTree* tree = (TTree*) fin->Get("Events");
   if (!tree) { printf("no tree 'Events' in %s\n", argv[1]); return 1; }

   Muon mu[2];
   mu[0].bind(tree, 1);
   mu[1].bind(tree, 2);
   if (!mu[0].present()) { printf("muon-1 branches not found\n"); return 1; }
   const int nslot = mu[1].present() ? 2 : 1;

   const bool haveBx = mu[0].bxspread.ok();
   if (!haveBx)
      printf("\n  WARNING: no bxspread branch on this file. The one-BX requirement\n"
             "           cannot be applied and every track is kept. Any result from\n"
             "           this output includes out-of-time background.\n\n");

   printf("=================================================================\n");
   printf(" in     %s\n", argv[1]);
   printf(" tree   Events, %lld entries, %d muon slot(s)\n", tree->GetEntries(), nslot);
   printf(" cuts   bxspread == 0%s,  pT > %.0f,  |eta| < %.2f,  nstub >= 3,\n",
          haveBx ? "" : " (NOT APPLIED)", PTMIN, ETAMAX);
   printf("        dxy < 1,  qual >= 12 (>=13 nstub 3, >=14 nstub 4)\n");
   printf("=================================================================\n");

   // ------------------------------------------------------------- histograms
   TH1D *hKallP = bookK("hK_plus_all",  "|K|, q>0, all sectors");
   TH1D *hKallM = bookK("hK_minus_all", "|K|, q<0, all sectors");

   TH1D *hSecP[NSEC], *hSecM[NSEC], *hRotP[NSEC], *hRotM[NSEC];
   for (int s = 0; s < NSEC; ++s) {
      hSecP[s] = bookK(Form("hK_plus_sec%02d",  s), Form("|K|, q>0, sector %d", s));
      hSecM[s] = bookK(Form("hK_minus_sec%02d", s), Form("|K|, q<0, sector %d", s));
      hRotP[s] = bookK(Form("hK_plus_rot%02d",  s), Form("|K|, q>0, rotated sector %d", s));
      hRotM[s] = bookK(Form("hK_minus_rot%02d", s), Form("|K|, q<0, rotated sector %d", s));
   }

   TH1D *hEtaP[NETA], *hEtaM[NETA];
   for (int e = 0; e < NETA; ++e) {
      hEtaP[e] = bookK(Form("hK_plus_eta%d",  e), Form("|K|, q>0, eta band %d", e));
      hEtaM[e] = bookK(Form("hK_minus_eta%d", e), Form("|K|, q<0, eta band %d", e));
   }

   // sector x eta half: the material for deciding if one number per sector is
   // enough, or if the offset also moves with the wheel
   TH1D *hSEP[NSEC][2], *hSEM[NSEC][2];
   for (int s = 0; s < NSEC; ++s)
      for (int e = 0; e < 2; ++e) {
         hSEP[s][e] = bookK(Form("hK_plus_s%02de%d",  s, e), Form("|K|, q>0, sector %d, eta %s", s, e ? ">0" : "<0"));
         hSEM[s][e] = bookK(Form("hK_minus_s%02de%d", s, e), Form("|K|, q<0, sector %d, eta %s", s, e ? ">0" : "<0"));
      }

   TH1D *hNsP[5], *hNsM[5];        // only 3 and 4 are filled
   for (int n = 3; n <= 4; ++n) {
      hNsP[n] = bookK(Form("hK_plus_nstub%d",  n), Form("|K|, q>0, nstub %d", n));
      hNsM[n] = bookK(Form("hK_minus_nstub%d", n), Form("|K|, q<0, nstub %d", n));
   }

   // phi vs |K|, so the notebook can build any phi binning it likes
   TH2D* h2P = new TH2D("h2_phiK_plus",  "q>0;#phi [rad];|K| [hw units]", 144, -TMath::Pi(), TMath::Pi(), 1600, 0, KMAX);
   TH2D* h2M = new TH2D("h2_phiK_minus", "q<0;#phi [rad];|K| [hw units]", 144, -TMath::Pi(), TMath::Pi(), 1600, 0, KMAX);
   TH2D* h2EP = new TH2D("h2_etaK_plus",  "q>0;#eta;|K| [hw units]", 84, -1.05, 1.05, 1600, 0, KMAX);
   TH2D* h2EM = new TH2D("h2_etaK_minus", "q<0;#eta;|K| [hw units]", 84, -1.05, 1.05, 1600, 0, KMAX);
   h2P->Sumw2(); h2M->Sumw2(); h2EP->Sumw2(); h2EM->Sumw2();

   // plain distributions, for sanity and for the notebook's first look
   TH1D* hSignedP = book1("hKsigned_plus",  "signed K, q>0;K [hw units];tracks", 1600, -KMAX, KMAX);
   TH1D* hSignedM = book1("hKsigned_minus", "signed K, q<0;K [hw units];tracks", 1600, -KMAX, KMAX);
   TH1D* hPtP     = book1("hPt_plus",   "q>0;p_{T} [GeV];tracks", 250, 0, 1000);
   TH1D* hPtM     = book1("hPt_minus",  "q<0;p_{T} [GeV];tracks", 250, 0, 1000);
   TH1D* hPhiP    = book1("hPhi_plus",  "q>0;#phi [rad];tracks", 144, -TMath::Pi(), TMath::Pi());
   TH1D* hPhiM    = book1("hPhi_minus", "q<0;#phi [rad];tracks", 144, -TMath::Pi(), TMath::Pi());
   TH1D* hEtaHP   = book1("hEta_plus",  "q>0;#eta;tracks", 100, -1.0, 1.0);
   TH1D* hEtaHM   = book1("hEta_minus", "q<0;#eta;tracks", 100, -1.0, 1.0);
   TH1D* hNstubH  = book1("hNstub",     "selected tracks;n stubs;tracks", 5, -0.5, 4.5);
   TH1D* hQualH   = book1("hQual",      "selected tracks;hw quality;tracks", 17, -0.5, 16.5);
   TH1D* hDxyH    = book1("hDxy",       "selected tracks;hw dxy;tracks", 100, 0, 4);
   TH1D* hBxAll   = book1("hBxspread_beforecut",
                          "all tracks reaching the BX cut;bxspread;tracks", 12, -0.5, 11.5);

   TH1D* hCut = book1("cutflow", "cutflow;;tracks", 9, 0.5, 9.5);
   const char* cutName[9] = {"read", "finite", "|eta|<0.83", "bxspread==0",
                             "nstub>=3", "dxy<1", "quality", "sign(K)==q", "pT>15"};
   for (int i = 0; i < 9; ++i) hCut->GetXaxis()->SetBinLabel(i + 1, cutName[i]);

   // ------------------------------------------------------------- event loop
   const double ROT = TMath::Pi() / 12.0;      // half a sector
   Long64_t nentries = tree->GetEntries();

   for (Long64_t i = 0; i < nentries; ++i) {
      tree->GetEntry(i);
      if (i % 10000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries);
      fflush(stdout);

      for (int m = 0; m < nslot; ++m) {
         hCut->Fill(1);

         double K = mu[m].K();
         double eta = mu[m].eta(), phi = mu[m].phi(), q = mu[m].charge();
         if (!std::isfinite(K) || !std::isfinite(eta) || !std::isfinite(phi)) continue;
         if (K == 0 || fabs(q) != 1) continue;
         hCut->Fill(2);

         if (fabs(eta) > ETAMAX) continue;
         hCut->Fill(3);

         int bx = mu[m].bxspread.ok() ? (int) mu[m].bxspread() : 0;
         hBxAll->Fill(bx);
         if (haveBx && bx != 0) continue;           // <-- the one-BX requirement
         hCut->Fill(4);

         int nstub = mu[m].nstub.ok() ? (int) mu[m].nstub() : 4;
         if (nstub < 3) continue;
         hCut->Fill(5);

         double dxy = mu[m].dxy.ok() ? mu[m].dxy() : 0;
         if (dxy >= 1) continue;
         hCut->Fill(6);

         int qual = mu[m].qual.ok() ? (int) mu[m].qual() : 15;
         if (qual < 12) continue;
         if (nstub == 3 && qual < 13) continue;
         if (nstub == 4 && qual < 14) continue;
         hCut->Fill(7);

         if ((K > 0) != (q > 0)) continue;
         hCut->Fill(8);

         double absK = fabs(K);
         if (absK >= KMAX || PtFromK(absK) < PTMIN) continue;
         hCut->Fill(9);

         // ---- fill ----
         const int s  = sectorOf(phi, 0.0);
         const int sr = sectorOf(phi, ROT);
         const int e  = etaBinOf(eta);
         const int eh = (eta > 0) ? 1 : 0;
         const double pt = PtFromK(absK);
         const bool pos = (q > 0);

         if (pos) {
            hKallP->Fill(absK); hSecP[s]->Fill(absK); hRotP[sr]->Fill(absK);
            hEtaP[e]->Fill(absK); hSEP[s][eh]->Fill(absK);
            if (nstub == 3 || nstub == 4) hNsP[nstub]->Fill(absK);
            h2P->Fill(phi, absK); h2EP->Fill(eta, absK);
            hSignedP->Fill(K); hPtP->Fill(pt); hPhiP->Fill(phi); hEtaHP->Fill(eta);
         } else {
            hKallM->Fill(absK); hSecM[s]->Fill(absK); hRotM[sr]->Fill(absK);
            hEtaM[e]->Fill(absK); hSEM[s][eh]->Fill(absK);
            if (nstub == 3 || nstub == 4) hNsM[nstub]->Fill(absK);
            h2M->Fill(phi, absK); h2EM->Fill(eta, absK);
            hSignedM->Fill(K); hPtM->Fill(pt); hPhiM->Fill(phi); hEtaHM->Fill(eta);
         }
         hNstubH->Fill(nstub); hQualH->Fill(qual); hDxyH->Fill(dxy);
      }
   }

   // ------------------------------------------------------------- write
   TFile* fout = TFile::Open(argv[2], "RECREATE");

   fout->cd();
   hKallP->Write(); hKallM->Write();

   TDirectory* dSec = fout->mkdir("sectors");    dSec->cd();
   for (int s = 0; s < NSEC; ++s) { hSecP[s]->Write(); hSecM[s]->Write(); }

   TDirectory* dRot = fout->mkdir("sectors_rot"); dRot->cd();
   for (int s = 0; s < NSEC; ++s) { hRotP[s]->Write(); hRotM[s]->Write(); }

   TDirectory* dEta = fout->mkdir("eta");         dEta->cd();
   for (int e = 0; e < NETA; ++e) { hEtaP[e]->Write(); hEtaM[e]->Write(); }

   TDirectory* dSE = fout->mkdir("sector_eta");   dSE->cd();
   for (int s = 0; s < NSEC; ++s)
      for (int e = 0; e < 2; ++e) { hSEP[s][e]->Write(); hSEM[s][e]->Write(); }

   TDirectory* dNs = fout->mkdir("nstub");        dNs->cd();
   for (int n = 3; n <= 4; ++n) { hNsP[n]->Write(); hNsM[n]->Write(); }

   TDirectory* dMap = fout->mkdir("maps");        dMap->cd();
   h2P->Write(); h2M->Write(); h2EP->Write(); h2EM->Write();

   TDirectory* dBas = fout->mkdir("basic");       dBas->cd();
   hSignedP->Write(); hSignedM->Write();
   hPtP->Write(); hPtM->Write();
   hPhiP->Write(); hPhiM->Write();
   hEtaHP->Write(); hEtaHM->Write();
   hNstubH->Write(); hQualH->Write(); hDxyH->Write();
   hBxAll->Write(); hCut->Write();

   // provenance, so a plot made weeks from now can be traced back
   fout->cd();
   TNamed prov("provenance",
               Form("input=%s;entries=%lld;slots=%d;bxcut=%d;ptmin=%.1f;etamax=%.2f;"
                    "lsb=%.8e;nsec=%d;kmax=%.0f;nkbins=%d",
                    argv[1], nentries, nslot, (int) haveBx, PTMIN, ETAMAX,
                    LSB, NSEC, KMAX, NK));
   prov.Write();

   // ------------------------------------------------------------- report
   printf("\n cutflow (tracks)\n");
   for (int i = 1; i <= 9; ++i)
      printf("   %-14s %10.0f\n", cutName[i - 1], hCut->GetBinContent(i));
   printf("\n selected   q>0 %.0f     q<0 %.0f\n", hKallP->Integral(), hKallM->Integral());
   printf(" wrote      %s\n", argv[2]);
   printf(" next       open MisalignmentFactor.ipynb and point it at that file\n");

   fout->Close();
   return 0;
}
