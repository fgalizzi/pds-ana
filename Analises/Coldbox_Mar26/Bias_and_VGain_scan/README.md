You can compile `VGainScans_ana.cpp`  and then quickly run analyses in parallel using GNU Parallel with the following command:

```bash
parallel -j 10 ./launch_ana.sh ::: config/config_module_*
```

# General
I assume the breakdown voltages as in the spreadsheet (tab: SiPM) which were estimated in LN2. To these values,
I add 0.2 V (0.1 V) for HPK (FBK) modules to correct for the LAr-LN2 difference.

Waveform selection is performed in the usual way: I discard waveforms which exceed the +-3$`\sigma`$ threshold
in the pre-pulse region. The $`\sigma`$ is computed as the RMS of the waveforms in the VGain-Noise scan runs.

## VD
### Module 1
The gaussians of the charge histogram look asymmetric with a weird tail on the left; especially at the highest overvoltage (+5 OV).
In addition, at low vgain $`\sigma_{cel}`$, the incremental width parameter, is larger than $`\sigma_{0}`$, the baseline width parameter, which is unexpected.

### Module 2
The gaussians of the charge histogram look asymmetric with a weird tail on the left; especially at the highest overvoltage (+5 OV).
In addition, at low vgain $`\sigma_{cel}`$, the incremental width parameter, is larger than $`\sigma_{0}`$, the baseline width parameter, which is unexpected.

## HD
Weird behaviour: the electronics noise is lower for the +7 OV (bias 874 DAC) runs of the bias scan. The noise components
that are reduced are the high-grequency bursts we observe in the waveforms. We cannot be 100% sure whether this is
actually due to the bias or it was an environmental effect. Actually, the data-taking was interrupted in the middle of
this scan... To investigate more looking at the data in the original folders.
#### About the reduction in noise
In this folder we put the good data originally sparse in different folders:
- "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/2026March/20260318_daphne-13_led_vgain_bias_scan_remote_rerun/"
    Bias 784 VGain [500, 2500], Bias 810 VGain [500, 2300]
- "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/2026March/20260318_daphne-13_led_vgain_bias_scan_remote_rerun_resume1/"
    Bias 810 VGain [2300, 2500], Bias 874 VGain [2300, 2500]
- "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/2026March/20260318_daphne-13_led_vgain_bias_scan_remote_rerun_resume2/"
    Bias 874 VGain [500, 2200]

We have vgain scans (500->2500, step of 100) at three bias voltage [784, 810, 874] in DAC unit.

We noted a reduction noise between the runs coming from the first folder (acquired at "2026-03-18T13:27:43.018788+00:00") and
the second (acquired at "2026-03-18T13:20:31.398904+00:00"). This is quite in coincidence with the message
"FEMBs have been turned back ON" in the slack channel #testing-crp-ehn1-daq, March 18th 2026

### Module 3
Quite poor SNR, we should probably go for more aggressive cuts given the high background rate.

### Module 4
Quite poor SNR, we should probably go for more aggressive cuts given the high background rate.

## PoF
For PoF modules I'm assuming that the overvoltage is 5 V.
As noise RMS I'm using the one of HD channel 12.
