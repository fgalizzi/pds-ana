#include "../classe.hpp"

// ****************************************************************
// Change this marco when you can easily loop over runs and you
// don't want to repeat the same commands every time!
// If the loop could be useful for the future, store it in Loops.cpp
// ****************************************************************

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////
string base_path = "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/Daphne_DAQ/";
///////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::TooAnnoying(){

  // Tuple input file, run, vgain, channel
  vector<tuple<string,int,int,int>> my_tuple;

  vector<int> runs;
  for(int i=9; i<30+9; i++){
    runs.push_back(33300+i);
  }

  int vgain = 100;
  vector<int> channels = {0,1,2,3,20,21,26,27};

  vector<string> ifiles;
  for(auto& run : runs){
    for(auto& ch : channels){
      string ifile = base_path+"run_"+to_string(run)+"/104/channel_"+to_string(ch)+".dat";
      my_tuple.push_back(make_tuple(ifile,run,vgain,ch));
    }
    vgain += 100;
  }

  // for(auto tuple : my_tuple) print_tuple(tuple);

  map<int,TFile> map_ch_FFTfile;
  for(auto& ch : channels){
    string filename= "FFT_VgainScan_Ch_"+to_string(ch)+".root";
    TFile tfile(filename.c_str(), "recreate");
    map_ch_FFTfile[ch] = tfile;
  }

  for(auto& tuple : my_tuple){
    wf_file = get<0>(tuple);
    read();
    TGraph* gNoise_spectral_density = build_avg_spectral_density(memorydepth,
      tick_len*memorydepth, tick_len, wfs, res);

    string gr_name = "VGain_"+to_string(get<2>(tuple));
    gNoise_spectral_density->SetName(gr_name.c_str());
    gNoise_spectral_density->SetTitle(gr_name.c_str());
    map_ch_FFTfile.second().Open();
    gNoise_spectral_density->write();
    map_ch_FFTfile.second().Close();
  }

/*  
  for(int i=0; i<ifiles.size(); i++){
    wf_file = ifiles[i];
    print = 1;
    trg_f = td_files[i];
    noise_f = "";
    muon_f = raw_files[i];
    Noise_PSD();
    print = 0;
    ite = 0;
    noise_f = trg_f;
    muon_f = clean_files[i];
    Noise_PSD();

  }
*/
}
