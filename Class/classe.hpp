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
    std::vector<std::vector<double>> wfs;    

    std::string data_format;
    bool invert;

    int memorydepth;     //Number of samples per waveform
    int n_wf;             //Number of WF in file
    int res;                //Resolution of the digitizer
    int prepulse_ticks;   //Template pre-pulse ticks
    double tick_len;    //In mu_s
    double spe_ampl;

    double sat_up;
    double sat_low;
    double bsl;

    //Calibration
    int int_low;         //Lower limit of wf integration
    int int_up;         //Upper fino a 3000
    int nbins;       // " " number of bins
    int nmaxpeaks;          // " " number of expected peaks
    int mu0_low;     //
    int mu0_up;       //Below an event is classified as noise
    int spe_low;  
    int spe_up;  
    int s0_low;      //
    int s0_up;
    int sc_low;      //
    int sc_up;
    double fit_low;
    double fit_up;
    double hmin;
    double hmax;
    TF1* fgaus;

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
  

    // DCR
    int win=10;
    int ite=2;
    double den=9.;

    //Saving data into files
    bool print = false;
    //Set calibration parameters manually
    bool manual = false;
    //Apply moving window on calib wfs
    bool mov_win = false;
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
    void Noise_PSD();
    void SPE();
    void upclass(cla A);
    void muonPDHD();
    void Avg_muon();
    void DCR();
    void SelfTrigger();
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
};

#include "../Header/G_Func.hpp"
#include "../Header/G_Read.hpp"
#include "../Header/G_WF.hpp"
#include "../Header/G_Utility.hpp"

#include "_c/set.cpp"
#include "_c/read.cpp"
#include "_c/Persistence.cpp"
#include "_c/AverageWF.cpp"
#include "_c/LED_Analysis.cpp"
#include "_c/Noise_PSD.cpp"
#include "_c/upclass.cpp"
#include "_c/muonPDHD.cpp"
#include "_c/Avg_muon.cpp"
#include "_c/SelfTrigger.cpp"
#include "_c/SPE.cpp"
#include "_c/DCR.cpp"
#include "_c/configDCR.cpp"
#include "_c/Saturation.cpp"
#include "_c/update.cpp"
#include "_c/LoadFitParameters.cpp"
#include "_c/TooAnnoying.cpp"
#include "_c/Jitter.cpp"
