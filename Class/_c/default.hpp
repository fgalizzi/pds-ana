std::string WF_FILE = "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/Daphne_DAQ/binaries/run_33753/channel_0.dat";
std::string TEMPL_F = "Template.dat";
std::string NOISE_F = "";
std::string MUON_F  = "FFT_20241204_M4_AFE0_Ch3_Vgain1000_clean.root";
std::string TRG_F  = "NoiseTd_20241204_M4_AFE0_Ch3_Vgain1000.dat";
// daphne caen csv 
std::string DATA_FORMAT  = "esteban";
bool INVERT              = 1;
int MEMORYDEPTH    = 1024;
int N_WF           ;
int RES            = 14;
int PREPULSE_TICKS = 124;
int PRETRG         = 230;
int AFTTRG         = 249;
double TICK_LEN    = 0.016;
 
double SAT_UP      = 200;
double SAT_LOW     = -445;
double BSL         = 50;
double RMS         = 40;
 
 
//calibration
int INT_LOW        = 126;
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
//ProtoDUNE
size_t CHANNEL     = 10900;
