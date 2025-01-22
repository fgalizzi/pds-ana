WF_FILE = "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/Daphne_DAQ/binaries/run_34330/channel_1.dat";
TEMPL_F = "Template.dat";
NOISE_F = "";
MUON_F  = "FFT_20241204_M4_AFE0_Ch3_Vgain1000_clean.root";
TRG_F   = "NoiseTd_20241204_M4_AFE0_Ch3_Vgain1000.dat";
// daphne caen csv 
DATA_FORMAT    = "esteban";
INVERT         = 1;
MEMORYDEPTH    = 1024;
N_WF           = 30000;
RES            = 14;
PREPULSE_TICKS = 124;
PRETRG         = 230;
AFTTRG         = 249;
TICK_LEN       = 0.016;
 
SAT_UP      = 1600;
SAT_LOW     = -445;
BSL         = 110;
RMS         = 40;
 
 
//calibration
INT_LOW        = 126;
INT_UP         = 155;
NBINS          = 300;
NMAXPEAKS      = 6;
MU0_LOW        = -56;
MU0_UP         = 41;
SPE_LOW        = 19;
SPE_UP         = 84;
S0_LOW         = 20;
S0_UP          = 180;
SC_LOW         = 5;
SC_UP          = 80;
FIT_LOW        = -200;
FIT_UP         = 1800;
HMIN           = -50;
HMAX           = 700;
SPE_AMPL       = 3.11725;
SPE_CHARGE     = 1883.09;
PEDESTAL       = 0;
//deconvolution
INT_PROMPT     = 1500;
ROLL           = 70;
AMP            = 70;
F_PROMPT       = 0.25;
AMP_LOW        = 550;
AMP_UP         = 2000;
DECO_SM        = 0.141051;
N2_            = 5e-06;
FIT_L          = 6.15;
FIT_U          = 18;
A_FAST         = 1643.64;
TAU_FAST       = 0.004;
A_SLOW         = 10.125;
TAU_SLOW       = 1.12;
SIGMA          = 0.03645;
T_0            = 0.042;
YERR           = 5;
//ProtoDUNE
CHANNEL        = 10900;
