set -e

# INDIR=/eos/user/f/flarover/DATA/AnalysisHSCP/modKBMTF/Zmumu
# OUTDIR=/eos/user/f/flarover/DATA/AnalysisHSCP/modKBMTF/Zmumu/output_modZmumu

INDIR=/eos/user/f/flarover/DATA/AnalysisHSCP/reKBMTF/Zmumu
OUTDIR=/eos/user/f/flarover/DATA/AnalysisHSCP/reKBMTF/Zmumu/output_modZmumu

mkdir -p $OUTDIR
mkdir -p $OUTDIR/plots

./Make.sh ZmumuAnalysis.cc
# ./ZmumuAnalysis.exe $INDIR/DY_Skim_SIM_all.root  $OUTDIR/DY_SIM.root DY KBMTF
# ./ZmumuAnalysis.exe $INDIR/DY_Skim.root $OUTDIR/DY.root data_obs KBMTF

./ZmumuAnalysis.exe $INDIR/reDY_SIM_SkimAll.root  $OUTDIR/DY_SIM.root DY KBMTF
./ZmumuAnalysis.exe $INDIR/reZmumu.root $OUTDIR/DY.root data_obs KBMTF


python3 Create_fake_Zmumu.py --selection="ZmumuKBMTF" --folder="re"
hadd -f $OUTDIR/datacard_ZmumuKBMTF.root $OUTDIR/DY.root $OUTDIR/DY_SIM.root $OUTDIR/Fake.root
python3 Draw_Zmumu.py --selection="ZmumuKBMTF" --folder="re"
