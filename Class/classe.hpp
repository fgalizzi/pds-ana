//
//List at the bottom the name of the .cpp macro you use
//


#pragma once

class cla{
  public:
    std::string wf_file;
    std::string templ_f;
    std::string muon_f;
    std::string noise_f;
    std::string trg_f;
    std::vector<std::vector<double>> wfs;    

    std::string data_format;
    bool invert;

    int memorydepth;      //Number of samples per waveform
    int n_wf;             //Number of WF in file
    int res;              //Resolution of the digitizer
    int prepulse_ticks;   //Template pre-pulse ticks
    int pretrg;           //Lower limit of acceptance window in selftrigger studies
    int afttrg;           //Upper " " " "
    double tick_len;      //In mu_s
    double spe_ampl;      
    double spe_charge;
    double pedestal;

    double sat_up;        //Saturation level - used to select wfs
    double sat_low;       //Low saturation level
    double bsl;           //To select wfs with baseline in [-bsl;+bsl]  
    double rms;           //Baseline RMS - used in baseline subtraction

    //Calibration
    int int_low;          //Lower limit of wf integration
    int int_up;           //Upper
    int nbins;            //Number of bins
    int nmaxpeaks;        //Number of expected peaks in calibration run
    //Extra parameters for manual mode
    int mu0_low;          //Lower bound of 0 pe peak integral
    int mu0_up;           //Upper - below an event is classified as noise
    int spe_low;          //Lower bound of spe peak - select spe cadidate 
    int spe_up;           //Upper
    int s0_low;           //Lower bound of 0 pe peak sigma
    int s0_up;            //Upper
    int sc_low;           //Lower bound of sigma channel
    int sc_up;            //Upper
    double fit_low;       //Fit lower limit 
    double fit_up;        //Upper
    double hmin;          //Histo lower limit 
    double hmax;          //Upper
    TF1* fgaus;           //Fit function

    //Deconvolution
    int int_prompt;
    double f_prompt;
    double n2_;
    double deco_sm;
    double amp_low;
    double amp_up;
    double fit_l;
    double fit_u;
    double a_fast;
    double tau_fast;
    double a_slow;
    double tau_slow;
    double sigma;
    double t_0; 
  
    // ProtoDUNE
    size_t channel;  

    // DCR
    int win=20;
    int ite=2;
    double den=9.;

    //Subtract the baseline according to "prepulse_ticks"
    bool sub_bsl = true;
    int sub_bsl_mode = 1;
    //Saving data into files
    bool print = false;
    //Display single waveforms
    bool display = false;
    //Set calibration parameters manually
    bool manual = false;
    //Apply moving window on calib wfs
    bool mov_win = false;
    //Apply matched/wiener filters (only in some macros)
    bool apply_filter = false;
    //Select wfs also for post-trigger features in SPE()
    bool pde_selection = false;
    //Deco Muon
    bool FFUNC = true;
    bool FIX_CONST = false;

    //Macro
    void set();
    void read();
    void Persistence();
    void AverageWF();
    void LED_Analysis();
    void Filt_Analysis();
    void Full_Resolution();
    void ProtoDUNE_Calibration();
    void Pdhd_FFT();
    void Noise_PSD();
    void SPE();
    void Muon_PDHD();
    void Avg_Muon();
    void DCR();
    void Self_Trigger();
    void configDCR();
    void Jitter();
    void Saturation();
    void update();    
    void LoadFitParameters(TF1* f);
    void TooAnnoying();    
    //Constructor
    cla(){set();}

  private:
    std::string oldwf_file;
    int oldprepulse_ticks;
    size_t oldchannel;
};

#define hdf5torootclass_cxx
#include "ProtoduneHD/hdf5torootclass.h"
#include "ProtoduneHD/wffunctions2.h"

#include "../Header/G_Func.hpp"
#include "../Header/G_Read.hpp"
#include "../Header/G_WF.hpp"
#include "../Header/G_Utility.hpp"

#include "_c/set.cpp"
#include "_c/read.cpp"
#include "_c/Persistence.cpp"
#include "_c/AverageWF.cpp"
#include "_c/LED_Analysis.cpp"
#include "_c/Filt_Analysis.cpp"
#include "_c/Full_Resolution.cpp"
#include "_c/Noise_PSD.cpp"
#include "_c/Muon_PDHD.cpp"
#include "_c/Avg_Muon.cpp"
#include "_c/Self_Trigger.cpp"
#include "_c/SPE.cpp"
#include "_c/DCR.cpp"
#include "_c/configDCR.cpp"
#include "_c/Saturation.cpp"
#include "_c/update.cpp"
#include "_c/LoadFitParameters.cpp"
#include "_c/ProtoDUNE_Calibration.cpp"
#include "_c/Pdhd_FFT.cpp"
#include "_c/TooAnnoying.cpp"
#include "_c/Jitter.cpp"
