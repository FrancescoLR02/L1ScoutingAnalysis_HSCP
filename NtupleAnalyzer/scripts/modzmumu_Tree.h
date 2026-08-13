#ifndef MODZMUMU_TREE_H
#define MODZMUMU_TREE_H

#include <iostream>
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TString.h"

// Some branches do not have the same scalar type in every ntuple version:
// dxy1/dxy2 are Int_t in the 2024 skims, Double_t in simulation and Float_t in data
// (same story for charge and qual, which are Int_t in the 2024 skims and Short_t now).
// Binding them with a fixed type reads garbage, so they go through this small wrapper
// which looks up the type on file and converts to double.
class ScalarBranch {
   public:
      void connect(TTree* tree, const char* name) {
         kind = kNone;
         buffer.d = 0;
         TBranch* branch = tree->GetBranch(name);
         if (!branch) {
            std::cout<<"WARNING: branch "<<name<<" not found, using 0"<<std::endl;
            return;
         }
         TLeaf* leaf = branch->GetLeaf(name);
         TString type = leaf ? leaf->GetTypeName() : "";
         if      (type=="Double_t") kind = kDouble;
         else if (type=="Float_t")  kind = kFloat;
         else if (type=="Int_t")    kind = kInt;
         else if (type=="UInt_t")   kind = kUInt;
         else if (type=="Short_t")  kind = kShort;
         else if (type=="UShort_t") kind = kUShort;
         else if (type=="Bool_t")   kind = kBool;
         else {
            std::cout<<"WARNING: branch "<<name<<" has unsupported type "<<type<<", using 0"<<std::endl;
            return;
         }
         tree->SetBranchAddress(name, (void*) &buffer);
      }

      double value() const {
         switch (kind) {
            case kDouble: return buffer.d;
            case kFloat:  return buffer.f;
            case kInt:    return buffer.i;
            case kUInt:   return buffer.u;
            case kShort:  return buffer.s;
            case kUShort: return buffer.us;
            case kBool:   return buffer.b;
            default:      return 0;
         }
      }

      operator double() const { return value(); }

   private:
      enum Kind { kNone, kDouble, kFloat, kInt, kUInt, kShort, kUShort, kBool };
      Kind kind = kNone;
      union Buffer {
         Double_t d;
         Float_t  f;
         Int_t    i;
         UInt_t   u;
         Short_t  s;
         UShort_t us;
         Bool_t   b;
      } buffer;
};

   unsigned int          run;
   unsigned int          luminosityBlock;
   unsigned int          bunchCrossing;
   unsigned int          orbitNumber;
   Double_t        mmumu;
   Bool_t          isOS;
   Double_t        DRmumu;
   Double_t        xsweight;
   Float_t         met;
   int           bxspread1;
   int           bxspread2;
   Bool_t          isL1MuMatched1;
   Bool_t          isL1MuMatched2;
   int           nstub1;
   int           nstub2;
   short stub1Station1;
   short stub2Station1;
   short stub3Station1;
   short stub4Station1;
   short stub1Station2;
   short stub2Station2;
   short stub3Station2;
   short stub4Station2;
   short stub1Bx1;
   short stub2Bx1;
   short stub3Bx1;
   short stub4Bx1;
   short stub1Bx2;
   short stub2Bx2;
   short stub3Bx2;
   short stub4Bx2;
   Float_t        genpt1;
   Float_t        genpt2;
   Double_t        pt1;
   Double_t        eta1;
   Double_t        phi1;
   ScalarBranch    charge1;
   ScalarBranch    qual1;
   ScalarBranch    dxy1;
   Double_t        pt2;
   Double_t        eta2;
   Double_t        phi2;
   ScalarBranch    charge2;
   ScalarBranch    qual2;
   ScalarBranch    dxy2;

#endif