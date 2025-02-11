# Overview
In this Coldbox we took data with three different acquisition system:
1. CAEN digitiser
2. DAPHNE in standalone mode
3. DAPHNE through the DAQ system

Hardware setup and runs log [here](https://docs.google.com/spreadsheets/d/1N9xcb2VVlzzDcNfBjlj_buhH9LiBTdG8-cnisb-orsI/edit?gid=631129411#gid=631129411).

Results [spreadsheet](https://docs.google.com/spreadsheets/d/1UbbC-N2yJ7k_QW4HT-eM1Flrzsd2dG26-lBIFAHH8Ew/edit?gid=1663548084#gid=1663548084).


## CAEN digitiser

## DAPHNE standalone

## DAPHNE DAQ
The following should be run in sequence.


### Noise_VGain_scan
It contains the single macro to analyze all the noise runs for the membrane modules (M1->M4). A bit of hard code here,
but, in principle, just need to modify `output_ana_folder` and `LED` variables in the macro.
The output of this macro is NECESSARY as input of the `Bias_and_VGain_scan` analysis. The run numbers are computed in a `for` loop.
The macro produces root files containing the FFTs and csv files with the RMS of the baseline for each channel and vgain.
```bash
root NoiseVGain_ana.cpp
```

### FineBias_scan
We performed a bias scan to estimate the V breakdown of the modules and to observe the dependence of few variables (e.g. gain,
cross-talk..) on the overvoltages. The analysis is performed with the `FineBiasScan_ana.cpp` macro in `FineBiasScan` folder.
In the `config` subfolder you can find the general configutation file `ana_config.json` and the specific configuration files
for each module `config_module_<module>.json`. In each config-file you'll see two bolcks of variables: one is to configure class members
and the other to configure variables proper of the specific analysis.
Repeat the analysis by running the following command:
```bash
source launch_ana.sh <module>
```
where `<module>` is the module number (1, 2, 3 or 4). You may need to re-analyze some runs in case the automatic fits fail. Check whether they
failed by looking at the output root files.
Then, condense the results in a root file with the `Gain_vs_VBias.cpp` macro (a bit hardcoded).

### Bias_and_VGain_scan
!! Configuration files depend on the output of the previous steps !!

For membrane modules (M1->M4) we scanned over 5 bias voltages and 30 VGain values for each bias.
To perform the analysis just run the following command:
```bash
source launch_ana.sh <module>
```
where `<module>` is the module number (1, 2, 3 or 4).
This will launch the `VGainScans_ana.cpp` macro. The configuration files are in the `Bias_and_VGain_scan/config` folder.
In  particular `ana_config.json` contains the general configuration for the analysis (e.g. inputs, outputs and the settings for the 
class that are shared for all the modules), while `config_module_<module>.json` contains the specific settings for the module `<module>`.

After that, use the VGainScan_ResultAnalyzer.cpp macro to create the root files with containing the plots of the
results. A bit of hard code here.
