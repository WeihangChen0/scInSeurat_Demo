#!/bin/bash
R pipeline/01_prepareInput.R
bash pipeline/02_run_genie3.sh
bash pipeline/03_runSCENIC.sh
