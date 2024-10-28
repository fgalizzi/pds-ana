#ifndef CLASSE_HPP
#define CLASSE_HPP

#include <string_view>
#ifndef my_headers_hpp
  #define my_headers_hpp
  #include "../Header/G_Func.hpp"
  #include "../Header/G_Read.hpp"
  #include "../Header/G_WF.hpp"
  #include "../Header/G_Utility.hpp"
#endif // !my_headers_hpp

#include <vector>
#include <string>
#include <cstddef>

//ROOT 
#include <TEfficiency.h>
#include <TF1.h>
#include <TH1D.h>
#include <TSpectrum.h>



using namespace std;



#pragma once


class cla{
  public:
    std::string class_path;
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
    double spe_under;
    double spe_charge;
    double sigma_zero;
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
    TF1* fgaus = nullptr; //Fit function

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
    size_t channel;   // e.g. 11245 <-> ep 112 ch 45
    size_t ep;        // Daphne - Endpoint

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
    void Jitter();
    void ST_Analysis();
    void configDCR();
    void Saturation();
    void update();    
    void LoadFitParameters(TF1* f);
    void TooAnnoying();    
    
    //Loops: to repeat the analysis on many files : )
    void Loop_ST_Analysis();
    void Loop_RMS_Analysis();

    //Constructor
    cla(){set();}

  private:
    // old = indicator to decide whether re-read the file or not
    std::string oldwf_file;
    int oldprepulse_ticks;
    size_t oldchannel;
    
    std::vector<std::vector<double>> trg_wf;
    void set(); //Initialize the class according to const.hpp, plot syle, fit preferences
    void read();//Read the wf_file and store the waveforms in wfs
    vector<size_t> read_pdhd_ch_map(int mask=0);
    vector<string> read_chs(string ch_file_name);
    TF1* set_charge_fit_function(TH1D* hI, TH1D* hFind=nullptr, bool avoid_auto_peak=false);
    double spe_ampl_correction;
    double compute_spe_correction(TF1* f);
    
    // Self trigger stuff
    void self_histos(TH1D* h_all, TH1D* h_trg, std::vector<double>& int_wf);
    int sf_bins = 100;      //Number of bins in self-trigger histos
    double sf_hmin = -2.;   //Lower limit [pe] self-trigger histos
    double sf_hmax =  7.;   //Upper " " " "

    // When running many channels from the same run
    size_t calibration_run; //Needed create the root file with the results
    size_t channel_low;     //Lower channel to look at (included)
    size_t channel_up;      //Upper " "
    int mask = 0;           //LED mask
    bool avoid_auto_peak = true;
    int pspe_low, pspe_up;  //ProtoDUNE spe_low and spe_up initial guess

};


#endif // !CLASSE_HPP
