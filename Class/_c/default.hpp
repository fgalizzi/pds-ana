std::string WF_FILE = "channel_1.dat";
std::string TEMPL_F = "";
std::string NOISE_F = "";
std::string MUON_F  = "";
std::string TRG_F   = "";
// daphne caen csv 
std::string DATA_FORMAT  = "daphne";
bool INVERT        = 1;
int MEMORYDEPTH    = 1024;
int N_WF           = 1000;
int RES            = 14;
int PREPULSE_TICKS = 20;
int PRETRG         = 10;
int AFTTRG         = 1000;
double TICK_LEN    = 0.016;
 
double SAT_UP      = 20000;
double SAT_LOW     = -20000;
double BSL         = 16000;
double RMS         = 4000;
 
 
//calibration
int INT_LOW        = 100;
int INT_UP         = 151;
int NBINS          = 300;
int NMAXPEAKS      = 6;
double MU0_LOW     = -56;
double MU0_UP      = 41;
double SPE_LOW     = 19;
double SPE_UP      = 84;
double S0_LOW      = 20;
double S0_UP       = 180;
double SC_LOW      = 5;
double SC_UP       = 80;
double FIT_LOW     = -200;
double FIT_UP      = 1800;
double HMIN        = -50;
double HMAX        = 700;
double SPE_AMPL    = 3.11725;
double SPE_CHARGE  = 137.102;
double PEDESTAL    = 0;
//deconvolution
int INT_PROMPT     = 1500;
int ROLL           = 70;
double AMP         = 70;
double F_PROMPT    = 0.25;
double AMP_LOW     = 550;
double AMP_UP      = 2000;
double DECO_SM     = 0.141051;
double N2_         = 5e-06;
double FIT_L       = 6.15;
double FIT_U       = 18;
double A_FAST      = 1643.64;
double TAU_FAST    = 0.004;
double A_SLOW      = 10.125;
double TAU_SLOW    = 1.12;
double SIGMA       = 0.03645;
double T_0         = 0.042;
double YERR        = 100;
//ProtoDUNE
size_t CHANNEL     = 10900;
