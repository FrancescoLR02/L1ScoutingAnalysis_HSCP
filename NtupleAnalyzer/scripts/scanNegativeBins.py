import ROOT
import sys


f = ROOT.TFile('/eos/user/f/flarover/DATA/AnalysisHSCP/reKBMTF/Combine/CombineClassification/Fake.root', "READ") 

print(f"Scanning slow.root for negative bins...\n")
print(f"{'Category (Folder)':<25} | {'Histogram Name':<45} | {'Bin':<5} | {'Value'}")
print("-" * 95)

found_negative = False

# Loop over the category folders (e.g., stub3_bx112_slow)
for dir_key in f.GetListOfKeys():
    if not dir_key.IsFolder(): continue
    d = dir_key.ReadObj()
    
    # Loop over the histograms inside the folder
    for hist_key in d.GetListOfKeys():
        h = hist_key.ReadObj()
        if h.InheritsFrom("TH1"):
            # Check every bin (from 1 to the max number of bins)
            for i in range(1, h.GetNbinsX() + 1):
                content = h.GetBinContent(i)
                if content < 0:
                    print(f"{d.GetName():<25} | {h.GetName():<45} | {i:<5} | {content:.6f}")
                    found_negative = True

f.Close()

if not found_negative:
    print("\nNo negative bins found. The file is clean!")
else:
    print("\nScan complete. These are the bins causing Combine to crash.")