set -e

INDIR=/eos/user/f/flarover/DATA/Drell-Yan
OUTDIR=output_ZmumuKBMTF

mkdir -p $OUTDIR
mkdir -p plots

./Make.sh ZmumuAnalysis.cc
./ZmumuAnalysis.exe $INDIR/modZmumuSim.root  $OUTDIR/DY.root            DY       KBMTF
./ZmumuAnalysis.exe $INDIR/modZmumuSim.root  $OUTDIR/DYall.root         DYall    KBMTF
./ZmumuAnalysis.exe $INDIR/modZmumuTEST.root $OUTDIR/Scouting_2025.root data_obs KBMTF
python3 Create_fake_Zmumu.py --selection="ZmumuKBMTF" --data="Scouting_2025"
hadd -f datacard_ZmumuKBMTF.root $OUTDIR/DY.root $OUTDIR/DYall.root $OUTDIR/Scouting_2025.root $OUTDIR/Fake.root
python3 Draw_Zmumu.py --selection="ZmumuKBMTF"
