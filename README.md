# dune-pd-ana

A C++ analysis framework to estimate the performance of the [DUNE](https://www.dunescience.org) Photon Detection System

## In this repo

Here you can find macro (file.cpp in folder [Macro](Macro)) useful to plot the
persistance of a set of waveforms, compute the Signal-to-Noise ratio from
a calibration run, compute the FFTs of the electronic noise... In
[Macro_description.txt](Macro_description.txt) you can find a brief overview of what the macros do.
"Header" contains the .hpp files where all the functions are delcared,
divided in categories. The flc.hpp file is the only one you should change to
adapt the codes to your dataset: it contains all the constant values needed for
your analysis.

To run a macro, you can use the command "root macro.cpp"

### Author
Federico Galizzi - f.galizzi1@campus.unimib.it
