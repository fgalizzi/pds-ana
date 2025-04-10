#ifndef all_the_includes
  #define all_the_includes
  // C++ general includes
  #include <iostream>
  #include "vector"
  #include <stdio.h>
  #include <fstream>
  #include <sstream>
  #include <filesystem>
  #include <ostream>
  #include <string>
  #include <vector>

  #include <TChain.h>
  #include <TString.h>
  #include <TF1.h>
  #include <TH1.h>
  #include <TH2.h>
  #include <TROOT.h>
  #include <TChain.h>
  #include <TFile.h>
  #include <TGraphErrors.h>
  #include <TFitResult.h>



#endif // !all_the_includes









#ifndef hdf5torootclass_cxx
  #define hdf5torootclass_cxx
  #include "../ProtoduneHD/hdf5torootclass.h"
  #include "../ProtoduneHD/wffunctions2.h"
#endif // !hdf5torootclass_cxx

#ifndef my_headers_hpp
  #define my_headers_hpp
  #include "../../Header/G_Func.hpp"
  #include "../../Header/G_Read.hpp"
  #include "../../Header/G_WF.hpp"
  #include "../../Header/G_Utility.hpp"
#endif // !my_headers_hpp
//
//
#ifndef default
  #include "default.hpp"
#endif // !default
//

#include "../classe.hpp"
#include "../private_methods.hpp"

#include "AverageWF.cpp"
#include "Avg_Muon.cpp"
#include "Avg_Alpha.cpp"
#include "DCR.cpp"
// #include "DUNEStyle.h"
// #include "FFT_Comp.cpp"
#include "Filt_Analysis.cpp"
#include "Full_Resolution.cpp"
#include "Jitter.cpp"
#include "LED_Analysis.cpp"
#include "LoadFitParameters.cpp"
#include "Loops.cpp"
#include "Optimize_Integration_Window.cpp"
// #include "MAX_Analysis.cpp"
// #include "MuonDeco.cpp"
#include "Muon_PDHD.cpp"
#include "Convolution.cpp"
#include "Noise_PSD.cpp"
// #include "Pdhd_FFT.cpp"
#include "Persistence.cpp"
#include "ProtoDUNE_Calibration.cpp"
#include "SPE.cpp"
#include "Build_Template.cpp"
#include "Saturation.cpp"
#include "Self_Trigger.cpp"
#include "ST_Analysis.cpp"
#include "TooAnnoying.cpp"
#include "configDCR.cpp"
// #include "map_channel_timestamp.cpp"
#include "read.cpp"
#include "set.cpp"
// #include "setup.cpp"
#include "update.cpp"
