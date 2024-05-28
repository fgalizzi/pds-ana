// This file contains constant declarations and
// the necessary #include (bottom)

const int MEMORYDEPTH = 5000;     //Number of samples per waveform
const int N_WF = 10000;             //Number of WF in file
const int RES = 14;                //Resolution of the digitizer
const int PREPULSE_TICKS = 1550;   //Template pre-pulse ticks
const double TICK_LEN = 0.002;    //In mu_s

const double SAT_UP = 110.;
const double SAT_LOW= -30.;
const double BSL = 18.;


//Calibration
const int INT_LOW = 1650;         //Lower limit of wf integration
const int INT_UP  = 3000;         //Upper fino a 3000
const int NBINS     = 2000;       // " " number of bins
const int NMAXPEAKS = 6;          // " " number of expected peaks
const int MU0_LOW = -2000.;     //
const int MU0_UP = 600.;       //This is the cut under which an event is classified as noise
const int SPE_LOW = 3960-1000;  // 0: 5000->7800
const int SPE_UP = 3960+1000;   // 1: 4200->6800
const int S0_LOW = 30;      //
const int S0_UP = 850;
const int SC_LOW = 20;      //
const int SC_UP = 580;
const double FIT_LOW = -3500;
const double FIT_UP = 25000;
const double HMIN = -4500;
const double HMAX = 45000;
//Deconvolution
const int INT_PROMPT = 7500;
const double F_PROMPT = 0.8;
/*
 const double C_FR = 1.;
 const double SPE_AMPL = 7.1;
 */

std::string WF_FILE ("hd_low_gain_no_transformer.dat");


#include "/Users/federico/dune-pd-ana/Header/G_Func.hpp"
#include "/Users/federico/dune-pd-ana/Header/G_Read.hpp"
#include "/Users/federico/dune-pd-ana/Header/G_WF.hpp"
#include "/Users/federico/dune-pd-ana/Header/G_Utility.hpp"


