You can compile `VGainScans_ana.cpp`  and then quickly run analyses in parallel using GNU Parallel with the following command:

```bash
parallel -j 10 ./launch_ana.sh ::: config/config_module_*
```

General
The wafeform selection is based on the baseline RMS of channel 1 because we have a good vgain scan
in noise runs only for few channels.
In the config files, very high vgain values were excluded because they lead to very poor SNR.

About M6
Default (reference) setting channels: 8 (1), 15 (2).
Neither channel 1 nor 2 produce good results when read by DAPHNE channel 15.
