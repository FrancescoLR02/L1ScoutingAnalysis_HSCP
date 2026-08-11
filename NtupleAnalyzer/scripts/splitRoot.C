// splitRoot.C
// Usage:
//   root -l -b -q 'splitRoot.C+("/eos/user/f/flarover/DATA/AnalysisHSCP/reKBMTF/DupRemove", "/eos/user/f/flarover/DATA/AnalysisHSCP/reKBMTF/SplitDupRemove", "Events")'

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TString.h>
#include <cmath>

void SplitOneFile(const char* fullpath, const char* outdir,
                   const char* treename = "tree",
                   Long64_t maxBytes = 2000000000LL /* 2GB */) {

    Long64_t fsize = 0;
    Long_t id, flags, modtime;
    if (gSystem->GetPathInfo(fullpath, &id, &fsize, &flags, &modtime) != 0) {
        printf("Cannot stat %s\n", fullpath);
        return;
    }

    if (fsize < maxBytes) {
        printf("%s (%.2f GB) below threshold, skipping\n", fullpath, fsize/1e9);
        return;
    }

    TFile *fin = TFile::Open(fullpath, "READ");   // read-only, original untouched
    if (!fin || fin->IsZombie()) { printf("Cannot open %s\n", fullpath); return; }

    TTree *tree = (TTree*)fin->Get(treename);
    if (!tree) { printf("No tree '%s' in %s\n", treename, fullpath); fin->Close(); return; }

    // Large read-ahead cache + cache all branches: batches EOS/xrootd reads
    // into big chunks instead of many small round-trips.
    tree->SetCacheSize(200*1024*1024);      // 200 MB
    tree->AddBranchToCache("*", kTRUE);

    Long64_t nentries = tree->GetEntries();
    int N = (int)std::ceil((double)fsize / (double)maxBytes);
    Long64_t chunk = nentries / N;

    TString base = gSystem->BaseName(fullpath);
    base.ReplaceAll(".root", "");

    printf("Splitting %s into %d parts (%lld entries) -> %s\n", fullpath, N, nentries, outdir);

    for (int i = 0; i < N; ++i) {
        Long64_t first = i * chunk;
        Long64_t nToCopy = (i == N-1) ? (nentries - first) : chunk;

        TString outname = Form("%s/%s_%d.root", outdir, base.Data(), i);
        TFile fout(outname, "RECREATE", "", fin->GetCompressionSettings());

        // "fast" -> basket-level copy, no decompress/recompress per entry
        TTree *newtree = tree->CopyTree("", "fast", nToCopy, first);
        newtree->Write("", TObject::kOverwrite);   // avoid extra autosave cycles
        fout.Close();
        printf("  -> %s (%lld entries)\n", outname.Data(), nToCopy);
    }

    fin->Close();   // original file was only ever opened READ, never written to
}

void splitRoot(const char* folder, const char* outdir, const char* treename = "tree") {
    if (gSystem->AccessPathName(outdir)) {
        gSystem->mkdir(outdir, kTRUE);
        printf("Created output folder %s\n", outdir);
    }

    TSystemDirectory dir("dir", folder);
    TList *files = dir.GetListOfFiles();
    if (!files) return;

    for (auto obj : *files) {
        TString fname = obj->GetName();
        if (!fname.EndsWith(".root")) continue;
        TString fullpath = Form("%s/%s", folder, fname.Data());
        SplitOneFile(fullpath, outdir, treename);
    }
}