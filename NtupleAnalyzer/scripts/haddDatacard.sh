

DATA_DIR="/eos/user/f/flarover/DATA/AnalysisHSCP/reKBMTF/Classification/AddedTracks_Classification"
ALL_DIR="/eos/user/f/flarover/DATA/AnalysisHSCP/reKBMTF/Combine/AddedTracks_Classification"
SHAPE_DIR="/eos/user/f/flarover/HSCP_2025/COMBINE/CMSSW_15_0_10/src/auxiliaries/shapes"

FILE_NAME="AddedTrackSlow"

# hadd -f ${ALL_DIR}/data_obs.root ${DATA_DIR}/*.root

# python3 Create_fakeFra.py --input ${ALL_DIR}/data_obs.root --output ${ALL_DIR}/Fake.root

hadd -f ${SHAPE_DIR}/${FILE_NAME}.root ${ALL_DIR}/*.root

