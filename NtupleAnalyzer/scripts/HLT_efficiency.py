#!/usr/bin/env python3
# HLT efficiency vs mass hypothesis
#
#   python3 hlt_eff_vs_mass.py [base_dir] [output.root]
#
# Expects one folder per mass point, e.g. HSCP_tauPrime_1000/*.root

import glob
import os
import re
from tqdm import tqdm
import sys

import ROOT

TREE    = "Events"          # tree name inside the root files
BRANCH  = "HLT_IsoMu24"     # 0 = not reconstructed, 1 = reconstructed
PATTERN = "HSCP_*"

base = sys.argv[1] if len(sys.argv) > 1 else "."
out  = sys.argv[2] if len(sys.argv) > 2 else "hlt_eff_vs_mass.root"

ROOT.gROOT.SetBatch(True)


def mass_of(folder):
    # last number in the folder name, in whatever unit the name uses
    return float(re.findall(r"[\d.]+", os.path.basename(folder))[-1])


points = []

for folder in tqdm(sorted(glob.glob(os.path.join(base, PATTERN)), key=mass_of)):
   files = glob.glob(os.path.join(folder, "*.root"))
   if not files:
      print("  no root files in %s, skipped" % folder)
      continue

   chain = ROOT.TChain(TREE)
   for f in files[:50]:
      chain.Add(f)

   tot = chain.GetEntries()
   if tot == 0:
      print("  empty tree in %s, skipped" % folder)
      continue

   passed = chain.GetEntries("%s == 1" % BRANCH)
   mass   = mass_of(folder)
   points.append((mass, passed, tot))

   print("  %-30s  %8d / %8d  eff = %.4f" %
         (os.path.basename(folder), passed, tot, float(passed) / tot))

if not points:
   sys.exit("no mass points found under '%s'" % base)

g = ROOT.TGraphAsymmErrors(len(points))
g.SetName(f"{PATTERN}")
g.SetTitle(";mass [GeV];HLT_IsoMu24 efficiency")

for i, (mass, passed, tot) in enumerate(points):
   eff = float(passed) / tot
   lo  = eff - ROOT.TEfficiency.ClopperPearson(tot, passed, 0.683, False)
   hi  = ROOT.TEfficiency.ClopperPearson(tot, passed, 0.683, True) - eff
   g.SetPoint(i, mass, eff)
   g.SetPointError(i, 0., 0., lo, hi)

fout = ROOT.TFile.Open(out, "RECREATE")
g.Write()
fout.Close()

print("\nwrote %d points to %s" % (len(points), out))