std::string wf_file = "0_wave0_45V00_0ADC_Ch0.dat";
// daphne caen csv 
std::string data_format  = "caen";
bool invert              = 0;
const int memorydepth    = 5000;
const int n_wf           = 50;
const int res            = 14;
const int prepulse_ticks = 1450;
const double tick_len    = 0.002;
 
const double sat_up      = 110;
const double sat_low     = -30;
const double bsl         = 15;
 
 
//calibration
const int int_low        = 1520;
const int int_up         = 1600;
const int nbins          = 500;
const int nmaxpeaks      = 6;
const double mu0_low     = -120;
const double mu0_up      = 63;
const double spe_low     = 449;
const double spe_up      = 572;
const double s0_low      = 20;
const double s0_up       = 180;
const double sc_low      = 5;
const double sc_up       = 80;
const double fit_low     = -200;
const double fit_up      = 1800;
const double hmin        = -200;
const double hmax        = 3000;
//deconvolution
const int int_prompt     = 0;
const double f_prompt    = 0;
