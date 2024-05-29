std::string WF_FILE = "channel_1.dat";
std::string TEMPL_F = "Template.dat";
std::string NOISE_F = "Noise.dat";
std::string MUON_F  = "Muon.dat";
// daphne caen csv 
std::string DATA_FORMAT  = "esteban";
bool INVERT              = 1;
const int MEMORYDEPTH    = 4096;
const int N_WF           = 3000;
const int RES            = 14;
const int PREPULSE_TICKS = 345;
const double TICK_LEN    = 0.016;
 
const double SAT_UP      = 800;
const double SAT_LOW     = -300;
const double BSL         = 100;
 
 
//calibration
const int INT_LOW        = 353;
const int INT_UP         = 393;
const int NBINS          = 300;
const int NMAXPEAKS      = 6;
const double MU0_LOW     = -169;
const double MU0_UP      = 74;
const double SPE_LOW     = 429;
const double SPE_UP      = 591;
const double S0_LOW      = 20;
const double S0_UP       = 180;
const double SC_LOW      = 5;
const double SC_UP       = 80;
const double FIT_LOW     = -200;
const double FIT_UP      = 1800;
const double HMIN        = -1000;
const double HMAX        = 9000;
const double SPE_AMPL    = 10;
const double SPE_CHARGE  = 10;
const double PEDESTAL    = 10;
//deconvolution
const int INT_PROMPT     = 70;
const double F_PROMPT    = 0.32;
const double AMP_LOW     = -8000;
const double AMP_UP      = 6000;
const double DECO_SM     = 0.141051;
const double N2_         = 1e-07;
const double FIT_L       = 0.23;
const double FIT_U       = 2.98;
const double A_FAST      = 5993.6;
const double TAU_FAST    = 0.004;
const double A_SLOW      = 25;
const double TAU_SLOW    = 1.12;
const double SIGMA       = 0.05;
const double T_0         = 0.8;
