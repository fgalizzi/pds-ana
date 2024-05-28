//
//List at the bottom the name of the .cpp macro you use
//


#pragma once

class cla{
  public:
    std::string WF_FILE;
    std::vector<double> WFS;    

    std::string DATA_FORMAT;
    bool INVERT;

    int MEMORYDEPTH;     //Number of samples per waveform
    int N_WF;             //Number of WF in file
    int RES;                //Resolution of the digitizer
    int PREPULSE_TICKS;   //Template pre-pulse ticks
    double TICK_LEN;    //In mu_s

    double SAT_UP;
    double SAT_LOW;
    double BSL;

    //Calibration
    int INT_LOW;         //Lower limit of wf integration
    int INT_UP;         //Upper fino a 3000
    int NBINS;       // " " number of bins
    int NMAXPEAKS;          // " " number of expected peaks
    int MU0_LOW;     //
    int MU0_UP;       //Below an event is classified as noise
    int SPE_LOW;  
    int SPE_UP;  
    int S0_LOW;      //
    int S0_UP;
    int SC_LOW;      //
    int SC_UP;
    double FIT_LOW;
    double FIT_UP;
    double HMIN;
    double HMAX;
    
    //Deconvolution
    int INT_PROMPT;
    double F_PROMPT;


    // DCR
    int WIN=10;
    int ITE=2;
    double DEN=9.;

    //Saving data into files
    bool PRINT = false;
    //Set calibration parameters manually
    bool MANUAL = false;

    //Macro
    void set();
    void read();
    void Persistence();
    void AverageWF();
    void LED_Analysis();
    void Noise_PSD();
    void SPE();
    void DCR();
    void configDCR();
    void Saturation();
    void update();    
    void LoadFitParameters(TF1* fgaus);
    
    //Constructor
    cla(){set();}
};

#include "/Users/federico/dune-pd-ana/Header2/G_Func.hpp"
#include "/Users/federico/dune-pd-ana/Header2/G_Read.hpp"
#include "/Users/federico/dune-pd-ana/Header2/G_WF.hpp"
#include "/Users/federico/dune-pd-ana/Header2/G_Utility.hpp"

#include "_c/set.cpp"
#include "_c/read.cpp"
#include "_c/Persistence.cpp"
#include "_c/AverageWF.cpp"
#include "_c/LED_Analysis.cpp"
#include "_c/Noise_PSD.cpp"
#include "_c/SPE.cpp"
#include "_c/DCR.cpp"
#include "_c/configDCR.cpp"
#include "_c/Saturation.cpp"
#include "_c/update.cpp"
#include "_c/LoadFitParameters.cpp"
