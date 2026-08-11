import os

# Define input and output directories
input_dir = '/eos/user/f/flarover/CRAB_DATA/L1ScoutingSelection/crab_ScoutingSelection_2025G/260421_081217/0000/'
output_dir = '/eos/user/f/flarover/DATA/Temp'
exe_script = "/eos/user/f/flarover/HSCP_2025/HSCPanalysis/CMSSW_15_0_10/src/L1ScoutingAnalysisRDataFrame/NtupleAnalyzer/scripts/FinalSelection_Zmumu.py"

xrootd_prefix = "root://eosuser.cern.ch/"

os.makedirs("logs_KBMTF", exist_ok=True)
os.makedirs("wrappers_KBMTF", exist_ok=True)


submit_template = """\
universe   = vanilla
executable = {wrapper}
output     = logs_KBMTF/{job_name}.out
error      = logs_KBMTF/{job_name}.err
log        = logs_KBMTF/{job_name}.log
+JobFlavour = "workday"
+AccountingGroup = "group_u_CMST3.all"
queue
"""

for input_file in os.listdir(input_dir):
    if input_file.endswith(".root"):
        job_name = os.path.splitext(input_file)[0]
        
        input_path_xrootd = f"{xrootd_prefix}{os.path.join(input_dir, input_file)}"
        output_path_xrootd = f"{xrootd_prefix}{os.path.join(output_dir, f'output_{input_file}')}"

        # Create the Wrapper Script (Saved in your AFS home)
        wrapper_script = f"wrappers_KBMTF/wrapper_{job_name}.sh"
        with open(wrapper_script, "w") as wf:
            wf.write(f"""#!/bin/bash
echo "Setting up CMS environment from EOS..."
source /cvmfs/cms.cern.ch/cmsset_default.sh

cd /eos/user/f/flarover/HSCP_2025/HSCPanalysis/CMSSW_15_0_10/src
eval `scramv1 runtime -sh`
cd $_CONDOR_SCRATCH_DIR

echo "Running analysis..."
python3 {exe_script} {input_path_xrootd} temp_output.root

echo "Transferring file back to EOS..."
xrdcp -f temp_output.root {output_path_xrootd}

rm temp_output.root
echo "Done!"
""")
        os.chmod(wrapper_script, 0o755)

        submit_file = f"submit_KBMTF_{job_name}.sub"
        with open(submit_file, "w") as f:
            f.write(submit_template.format(
                wrapper=wrapper_script,
                job_name=job_name
            ))

        os.system(f"condor_submit {submit_file}")

print("All jobs successfully submitted from AFS!")