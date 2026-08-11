#ifndef basic_sel_H
#define basic_sel_H

#include <ROOT/RDataFrame.hxx>
#include "TLorentzVector.h"

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include <cmath>
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
using namespace std;
using namespace ROOT;
using namespace ROOT::VecOps;

int GetIndex_nostub(int rank, int ncand, ROOT::VecOps::RVec<Float_t> &LepCand_pt, ROOT::VecOps::RVec<Float_t> &LepCand_eta, ROOT::VecOps::RVec<Float_t> &LepCand_phi);

int GetIndex(int rank, int ncand, ROOT::VecOps::RVec<Float_t> &LepCand_pt, ROOT::VecOps::RVec<Float_t> &LepCand_eta, ROOT::VecOps::RVec<Float_t> &LepCand_phi, ROOT::VecOps::RVec<Short_t> &nstub, ROOT::VecOps::RVec<Short_t> &bx1, ROOT::VecOps::RVec<Short_t> &bx2, ROOT::VecOps::RVec<Short_t> &bx3, ROOT::VecOps::RVec<Short_t> &bx4);//, ROOT::VecOps::RVec<Float_t> &LepCand_hwK);

int GetIndex_nostub_hwK(int rank, int ncand, ROOT::VecOps::RVec<Float_t> &LepCand_pt, ROOT::VecOps::RVec<Float_t> &LepCand_eta, ROOT::VecOps::RVec<Float_t> &LepCand_phi);//, ROOT::VecOps::RVec<Float_t> &LepCand_hwK);

TLorentzVector GetLepVector(int index, ROOT::VecOps::RVec<Float_t> &LepCand_pt, ROOT::VecOps::RVec<Float_t> &LepCand_eta,ROOT::VecOps::RVec<Float_t> &LepCand_phi);

TLorentzVector GetLepVector_hwK(int index, ROOT::VecOps::RVec<Float_t> &LepCand_eta,ROOT::VecOps::RVec<Float_t> &LepCand_phi, ROOT::VecOps::RVec<Float_t> &LepCand_pt);

int GetBxSpread(int ncand, int index, ROOT::VecOps::RVec<Short_t> &nstub, ROOT::VecOps::RVec<Short_t> &bx1, ROOT::VecOps::RVec<Short_t> &bx2, ROOT::VecOps::RVec<Short_t> &bx3, ROOT::VecOps::RVec<Short_t> &bx4);

int GetFirstBx(int ncand, int index, ROOT::VecOps::RVec<Short_t> &nstub, ROOT::VecOps::RVec<Short_t> &bx1, ROOT::VecOps::RVec<Short_t> &bx2, ROOT::VecOps::RVec<Short_t> &bx3, ROOT::VecOps::RVec<Short_t> &bx4);

int GetStationSpread(int ncand, int index, ROOT::VecOps::RVec<Short_t> &nstub, ROOT::VecOps::RVec<Short_t> &st1, ROOT::VecOps::RVec<Short_t> &st2, ROOT::VecOps::RVec<Short_t> &st3, ROOT::VecOps::RVec<Short_t> &st4);

int GetNstub(int ncand, int index, ROOT::VecOps::RVec<Short_t> &nstub);

bool IsL1MuMatched(int ncand, int index, ROOT::VecOps::RVec<Float_t> &pt, ROOT::VecOps::RVec<Float_t> &eta,ROOT::VecOps::RVec<Float_t> &phi, int nl1, ROOT::VecOps::RVec<Float_t> &l1pt, ROOT::VecOps::RVec<Float_t> &l1eta,ROOT::VecOps::RVec<Float_t> &l1phi);

float GetMET(int diff, int firstbx, int track_bx, float met0, float metm1, float metm2, float metm3, float metm4, float metm5, float metm6, float metm7, float metm8, float metm9);

bool IsColliding(int run, int bx);

bool IsEarlierColliding(int run, int bx, int interval, bool is_colliding);

float Get_newpt(int oldK);

//float Get_genbeta(float eta, float phi, int ngen, ROOT::VecOps::RVec<Float_t> &Gen_eta, ROOT::VecOps::RVec<Float_t> &Gen_phi, ROOT::VecOps::RVec<Short_t> &Gen_pdgid, ROOT::VecOps::RVec<Float_t> &Gen_beta);
template <typename T>
float Get_genbeta(float eta, float phi, int ngen,
                  const ROOT::VecOps::RVec<Float_t> &Gen_eta,
                  const ROOT::VecOps::RVec<Float_t> &Gen_phi,
                  const ROOT::VecOps::RVec<Short_t> &Gen_pdgid,
                  const ROOT::VecOps::RVec<T> &Gen_val) {
    float out = -1.0f;
    TLorentzVector my_reco;
    my_reco.SetPtEtaPhiM(100.0, eta, phi, 1.0);
    for (int i = 0; i < ngen; ++i) {
        const long id = std::abs(static_cast<long>(Gen_pdgid[i]));
        if (id == 13 || id == 17 || id == 16975 || id == 1000015 ||
            (id >= 1000993 && id <= 1093334)) {
            TLorentzVector tmp_gen;
            tmp_gen.SetPtEtaPhiM(100.0, Gen_eta[i], Gen_phi[i], 1.0);
            if (tmp_gen.DeltaR(my_reco) < 0.3) out = static_cast<float>(Gen_val[i]);
        }
    }
    return out;
}

float Get_genpt(float eta, float phi, int ngen, ROOT::VecOps::RVec<Float_t> &Gen_eta, ROOT::VecOps::RVec<Float_t> &Gen_phi, ROOT::VecOps::RVec<Short_t> &Gen_pdgid, ROOT::VecOps::RVec<Float_t> &Gen_pt);

float Get_geneta(float eta, float phi, int ngen, ROOT::VecOps::RVec<Float_t> &Gen_eta, ROOT::VecOps::RVec<Float_t> &Gen_phi, ROOT::VecOps::RVec<Short_t> &Gen_pdgid);

#endif
