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

### Module 3
Quite poor SNR, we should probably go for more aggressive cuts given the high background rate.

### Module 4
Quite poor SNR, we should probably go for more aggressive cuts given the high background rate.

## PoF
For PoF modules I'm assuming that the overvoltage is 5 V.
As noise RMS I'm using the one of HD channel 12.
