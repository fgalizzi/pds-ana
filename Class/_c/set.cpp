//
//Set all the constant values I need from cost.hpp
//

void cla::set(){
  wf_file        = WF_FILE;
  oldwf_file     = WF_FILE;
  data_format    = DATA_FORMAT;
  invert         = INVERT;
  templ_f        = TEMPL_F;
  noise_f        = NOISE_F;
  muon_f         = MUON_F;
  trg_f          = TRG_F;

  memorydepth       = MEMORYDEPTH;
  n_wf              = N_WF;
  res               = RES;
  prepulse_ticks    = PREPULSE_TICKS;
  oldprepulse_ticks = PREPULSE_TICKS;
  pretrg            = PRETRG;
  afttrg            = AFTTRG;
  tick_len          = TICK_LEN;

  sat_up         = SAT_UP;
  sat_low        = SAT_LOW;
  bsl            = BSL;
  rms            = BSL/3;

  //Calibration
  int_low   = INT_LOW;
  int_up    = INT_UP;
  nbins     = NBINS;
  nmaxpeaks = NMAXPEAKS;
  mu0_low   = MU0_LOW;
  mu0_up    = MU0_UP;
  spe_low   = SPE_LOW;
  spe_up    = SPE_UP;
  s0_low    = S0_LOW;
  s0_up     = S0_UP;
  sc_low    = SC_LOW;
  sc_up     = SC_UP;
  fit_low   = FIT_LOW;
  fit_up    = FIT_UP;
  hmin      = HMIN;
  hmax      = HMAX;
  spe_ampl  = SPE_AMPL; 
  spe_charge= SPE_CHARGE;  
  //Deconvolution
  int_prompt = INT_PROMPT;
  f_prompt   = F_PROMPT;
  amp_low    = AMP_LOW;
  amp_up     = AMP_UP;
  n2_        = N2_;
  deco_sm    = DECO_SM;
  fit_l      = FIT_L;
  fit_u      = FIT_U;
  a_fast     = A_FAST;
  tau_fast   = TAU_FAST; 
  a_slow     = A_SLOW;   
  tau_slow   = TAU_SLOW; 
  sigma      = SIGMA;
  t_0        = T_0;      

  //ProtoDUNE
  channel = CHANNEL;


  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 
}

