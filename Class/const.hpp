std::string WF_FILE = "";
std::string TEMPL_F = "Template.dat";
std::string NOISE_F = "Noise.dat";
std::string MUON_F  = "Muon.dat";
// daphne caen csv 
std::string DATA_FORMAT  = "csv";
bool INVERT              = 1;
const int MEMORYDEPTH    = 1000;
const int N_WF           = 10000;
const int RES            = 14;
const int PREPULSE_TICKS = 1450;
const double TICK_LEN    = 0.002;
 
const double SAT_UP      = 110;
const double SAT_LOW     = -30;
const double BSL         = 15;
 
 
//calibration
const int INT_LOW        = 1520;
const int INT_UP         = 1600;
const int NBINS          = 500;
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
const double HMIN        = -200;
const double HMAX        = 3000;
const double SPE_AMPL    = 10;
//deconvolution
const int INT_PROMPT     = 1550;
const double F_PROMPT    = 0.8;
const double AMP_LOW     = 400; 
const double AMP_UP      = 4000; 
const double N2_         = 1.e-7;
const double DECO_SM     = 0.08;
