//! Apply a LUT for the corrections on the curvature hwK

#ifndef SIM_KCORRECTION_H
#define SIM_KCORRECTION_H



#include <cmath>
#include <algorithm>
#include "TMath.h"

namespace sim_kcorr {

const int    NSEC   = 12;
const int    NETA   = 5;
const double ETAMAX = 0.83;
const double LSB    = 1.25 / double(1 << 13);

// DT/BMTF sector: bin k = sector k+1, centred at k*30 deg, edges at +-15, 45, ... deg
inline int sectorBin(double phi){
   const double w = 2*TMath::Pi()/NSEC;
   double x = std::fmod(phi + 0.5*w, 2*TMath::Pi());
   if (x < 0) x += 2*TMath::Pi();
   return std::min(std::max(int(x / w), 0), NSEC-1);
}

inline int etaBin(double eta){
   int b = int((eta + ETAMAX) / (2*ETAMAX) * NETA);
   return std::min(std::max(b, 0), NETA-1);
}

//! #########################################----FOR MC----#########################################
const double SIM_DELTA_PHI[12] = {
    -0.5808543299800815,  -0.4482559634934374,  -0.39470292489693337,
    -0.32173030915407824, -0.29121521935873984,  0.9848557834538392,
     0.5273392636005425,   0.6048884361602552,  -0.07746586458272894,
     0.822882273931358,    0.7732244600819353,  -0.805566929640018
};


const double SIM_DELTA_ETA[NETA] = {
   -1.0229336713155375, -0.4689045276261802, -0.08774475563269833,
    0.7737097341785608,   0.889865916802925
};

const double SIM_DELTA_NSTUB[3] = {  -0.438, 0.434, -0.214 };  


inline double deltaNStub(int nstub){
   int k = nstub - 2;
   return (k >= 0 && k < 3) ? SIM_DELTA_NSTUB[k] : 0.;
}

enum Level { kNone = 0, kPhi = 1, kPhiEta = 2, kFull = 3 };

// signed hwK in, signed corrected hwK out
inline double correctK(double K, double phi, double eta, int nstub, int level = kFull){
   double k = K;
   if (level >= kPhi)    k -= SIM_DELTA_PHI[sectorBin(phi)];
   if (level >= kPhiEta) k -= SIM_DELTA_ETA[etaBin(eta)];
   if (level >= kFull)   k -= deltaNStub(nstub);
   return k;
}

inline double applyMaterialMap(double K){

   double FK = std::fabs(K);
   float lsb = 1.25 / float(1 << 13);
   FK = FK * lsb;
   FK = .8569 * FK / (1.0 + 0.1144 * FK);
   return FK / lsb;
}

inline int chargeFromK(double K){ return (K >= 0) ? +1 : -1; }

// original KBMTF LUT: material/B-field, global misalignment, BMTF scale
inline double ptLUT_orig(double K){
   int charge = chargeFromK(K);
   double FK = std::fabs(K);
   if (FK > 2047) FK = 2047.;
   if (FK < 9)    FK = 9.;
   FK = FK * LSB;
   FK = .8569 * FK / (1.0 + 0.1144 * FK);
   FK = FK - charge * 1.23e-03;
   //FK = FK / 1.17;
   double pt = 0;
   if (FK != 0) pt = 1 / FK;
   if (pt < 4) pt = 4;
   if (pt > 1000) pt = 1000;
   return pt;
}


float Get_newpt(int oldK){
  float K = oldK-9;
  if (K==0) K=1;
  float lsb = 1.25 / float(1 << 13);
  float FK = abs(K);

  if (FK > 2047)
    FK = 2047.;

  FK = FK * lsb;

  //step 1 -material and B-field
  FK = .8569 * FK / (1.0 + 0.1144 * FK);

  float pt = 0;
  if (FK != 0)
    pt = float(2.0 / FK);

  if (pt > 2000)
    pt = 2000;

  if (pt < 8)
    pt = 8;

  return pt/2;
}

// new LUT: Kcorr must already come from correctK()
// applyMaterial = the material/B-field step of the original LUT (physical scale, as in Get_pTfromK)
inline double ptLUT_corr(double Kcorr, bool applyMaterial = true){
   double FK = std::fabs(Kcorr);
   if (FK > 2047) FK = 2047.;
   if (FK < 3)    FK = 3.;
   FK = FK * LSB;
   if (applyMaterial) FK = .8569 * FK / (1.0 + 0.1144 * FK);
   double pt = 0;
   if (FK != 0) pt = 1 / FK;
   if (pt < 4) pt = 4;
   if (pt > 1000) pt = 1000;
   return pt;
}

} // namespace kcorr

#endif