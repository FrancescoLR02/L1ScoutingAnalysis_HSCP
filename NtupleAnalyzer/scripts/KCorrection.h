//! Apply a LUT for the corrections on the curvature hwK

#ifndef KCORRECTION_H
#define KCORRECTION_H



#include <cmath>
#include <algorithm>
#include "TMath.h"

namespace kcorr {

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

// step 1: measured on raw hwK
const double DELTA_PHI[NSEC] = {
    -4.517551982392351,  15.357480990344975,  16.890352083459184,
     4.568931798089349,  -8.445488610821378,  -3.2565913333840792,
     6.753763570952405,  17.763404130940582,  14.926302197139007,
     0.5210699539132264, -14.19312748508228,  -18.101331219109678
};

// step 2: measured on phi-corrected hwK
const double DELTA_ETA[NETA] = {
   -5.031174684878439, -1.7392444643104927, 0.9836903473143631,
    2.591548000863815,  2.8348744112647695
};

// step 3: measured on phi+eta corrected hwK, nStub = 2, 3, 4
// NOTE: re-extract after fixing the muon-2 fill (hwK2_lut -> hwK2_lut2) in the Zmumu script
const double DELTA_NSTUB[3] = { 0.905, -0.490, 0.127 };

inline double deltaNStub(int nstub){
   int k = nstub - 2;
   return (k >= 0 && k < 3) ? DELTA_NSTUB[k] : 0.;
}

enum Level { kNone = 0, kPhi = 1, kPhiEta = 2, kFull = 3 };

// signed hwK in, signed corrected hwK out
inline double correctK(double K, double phi, double eta, int nstub, int level = kFull){
   double k = K;
   if (level >= kPhi)    k -= DELTA_PHI[sectorBin(phi)];
   if (level >= kPhiEta) k -= DELTA_ETA[etaBin(eta)];
   if (level >= kFull)   k -= deltaNStub(nstub);
   return k;
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
   return pt;
}


float Get_newpt(double oldK){
  float K = static_cast<int>(oldK)-9;
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
inline double ptLUT_corr(double Kcorr, bool applyMaterial = false){
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