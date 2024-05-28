//
//Set all the constant values I need from cost.hpp
//

void cla::set(){
  #include "const.hpp"
  WF_FILE = wf_file;
  DATA_FORMAT = data_format;
  INVERT = invert;

  MEMORYDEPTH = memorydepth;
  N_WF = n_wf;           
  RES = res;           
  PREPULSE_TICKS = prepulse_ticks;
  TICK_LEN = tick_len;  

  SAT_UP = sat_up;
  SAT_LOW = sat_low;
  BSL = bsl;

  //Calibration
  INT_LOW = int_low;     
  INT_UP = int_up;     
  NBINS = nbins;     
  NMAXPEAKS = nmaxpeaks;
  MU0_LOW = mu0_low; 
  MU0_UP = mu0_up; 
  SPE_LOW = spe_low;
  SPE_UP = spe_up;
  S0_LOW = s0_low;
  S0_UP = s0_up;
  SC_LOW = sc_low;
  SC_UP = sc_up;
  FIT_LOW = fit_low;
  FIT_UP = fit_up;
  HMIN = hmin;
  HMAX = hmax;
  //Deconvolution
  //INT_PROMPT;
  //F_PROMPT;

}

