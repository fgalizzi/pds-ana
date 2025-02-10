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
### Noise_VGain_scan
It contains the single macro to analyze all the noise runs for the membrane modules (M1->M4). A bit of hard code here,
but, in princliple, just need to modify `outoul_ana_folder` and `LED` variables in the macro.
The output of this macro is NECESSARY as input of the `Bias_and_VGain_scan` analysis.

### Bias_and_VGain_scan
For membrane modules (M1->M4) we scanned over 5 bias voltages and 30 VGain values for each bias.
To perform the analysis just run the following command:
```bash
./Bias_and_VGain_scan/launch_ana.sh <module>
```
where `<module>` is the module number (1, 2, 3 or 4).
This will launch the VGainScans_ana.cpp macro. The configuration files are in the `Bias_and_VGain_scan/config` folder.
In  particular `ana_config.json` contains the general configuration for the analysis (e.g. inputs, outputs and the settings for the 
class that are shared for all the modules), while `config_module_<module>.json` contains the specific settings for the module `<module>`.

After that, use the VGainSCan_ResultAnalyzer.cpp macro to create the root files with containing the plots of the
results. A bit of hard code here.
