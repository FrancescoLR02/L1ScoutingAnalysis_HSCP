#include <iostream>
#include <string>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"

using namespace std;
namespace fs = std::filesystem;

// Compile:
// g++ -O3 -std=c++17 count_entries.cc -o count_entries.exe $(root-config --cflags --glibs)
//
// Run:
// ./count_entries.exe /path/to/directory

int main(int argc, char** argv) {

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " /path/to/directory" << endl;
        return 1;
    }

    string directory = argv[1];

    Long64_t total_entries = 0;
    Long64_t nfiles = 0;
    Long64_t badfiles = 0;

    for (const auto& entry : fs::directory_iterator(directory)) {

        if (!entry.is_regular_file())
            continue;

        string filename = entry.path().string();

        // Only ROOT files
        if (entry.path().extension() != ".root")
            continue;

        nfiles++;

        TFile *file = TFile::Open(filename.c_str(), "READ");

        if (!file || file->IsZombie()) {
            cerr << "ERROR opening: " << filename << endl;
            badfiles++;
            continue;
        }

        TTree *tree = dynamic_cast<TTree*>(file->Get("Events"));

        if (!tree) {
            cerr << "ERROR: Events tree not found in: "
                 << filename << endl;
            badfiles++;
            file->Close();
            delete file;
            continue;
        }

        Long64_t nentries = tree->GetEntries();

        total_entries += nentries;

        // Print progress rather than thousands of lines
        if (nfiles % 100 == 0) {
            cout << "Processed " << nfiles
                 << " files, current total = "
                 << total_entries << endl;
        }

        file->Close();
        delete file;
    }

    cout << endl;
    cout << "========================================" << endl;
    cout << "Files processed : " << nfiles << endl;
    cout << "Bad files       : " << badfiles << endl;
    cout << "Total entries   : " << total_entries << endl;
    cout << "========================================" << endl;

    return 0;
}