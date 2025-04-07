# Human Task Analysis Pipeline: MD-dlPFC Attentional Control Task

This repository contains the analysis pipeline for the human version of the translational attentional control task described in our STAR Protocols submission. The task assesses executive function under sensory uncertainty and congruency, with data collected from both healthy control participants and individuals with schizophrenia.

The pipeline is designed to analyze behavioral performance, detect high-performance "bouts," fit psychometric functions, and compare group-level metrics across ambiguity and congruency conditions.

## Repository Contents

- `CueUncertaintyTask.7z`  
  Contains the PsychoPy task code and related materials for running the human behavioral experiment.

- `human_amb_data_Anna_Huang.7z`  
  Contains MATLAB scripts and helper functions for processing and analyzing participant data. You will need to extract this archive and provide the correct path to your behavioral data for analysis.

- `example_data/`  
  Includes an example of simulated data for testing and demonstration purposes (real patient data not publicly available due to privacy concerns).

## Requirements

- **MATLAB R2021b** or later
- Required MATLAB Toolboxes:
  - Statistics and Machine Learning Toolbox
- Required external MATLAB scripts (available via Add-ons or GitHub):
  - `shadedErrorBar` for plotting performance curves with confidence intervals
  - `chi2test` for significance testing

Make sure these dependencies are installed and added to your MATLAB path.

## Setup Instructions

**Edit Paths**
   - Update the path variables near the top of the script to point to your extracted data folders:
     ```matlab
     dataPath = 'C:/path/to/your/data/';
     savePath = 'C:/path/to/save/outputs/';
     ```

## Contact

For questions about the analysis pipeline, please contact:

- **Sahil Suresh** – [sahil.suresh@tufts.edu](mailto:sahil.suresh@tufts.edu)

## Citation

If you use this code, please cite the STAR Protocol associated with this analysis (reference will be added upon publication).

---

**Disclaimer:** This code is provided as-is for academic research purposes only.
