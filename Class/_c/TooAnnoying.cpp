#include "../classe.hpp"
#include <cstddef>
#include <numeric>
#include <string>
#include <vector>

// ****************************************************************
// Change this marco when you can easily loop over runs and you
// don't want to repeat the same commands every time!
// If the loop could be useful for the future, store it in Loops.cpp
// ****************************************************************

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////
string base_path   = "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/Daphne_DAQ/";
string fft_outfile = "FFT_VGainScans_Membrane.root";
string ana_outfile = "VGain_RMS_Membrane.csv";
///////////////////////////////////////////////////////////////////

double standard_deviation_vec_vec(vector<vector<double>>& wfs){
  double variance = 0.;
  double len = double(wfs[0].size());
  for(auto& wf : wfs){
    double mean = accumulate(wf.begin(),wf.end(),0.)/len;
    for(auto& e : wf) variance += (e-mean)*(e-mean);
  }
  return sqrt(variance / double(len*wfs.size()));
}

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

  TFile fft_ofile(fft_outfile.c_str(), "recreate");
  for(auto& ch : channels){
    string dir_name = "Ch_"+to_string(ch);
    fft_ofile.mkdir(dir_name.c_str());
  }

  for(auto& tuple : my_tuple){
    wf_file = get<0>(tuple);
    read();
    TGraph* gNoise_spectral_density = build_avg_spectral_density(memorydepth,
      tick_len*memorydepth, tick_len, wfs, res);

    string gr_name = "VGain_"+to_string(get<2>(tuple));
    gNoise_spectral_density->SetName(gr_name.c_str());
    gNoise_spectral_density->SetTitle(gr_name.c_str());
    fft_ofile.cd(to_string(get<3>(tuple)).c_str());
    gNoise_spectral_density->Write();
    if (print == true){
        vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 

        feature_value.push_back({"Run", get<1>(tuple)});
        feature_value.push_back({"VGain" , get<2>(tuple)});
        feature_value.push_back({"Ch", get<3>(tuple)});
        feature_value.push_back({"Rms", standard_deviation_vec_vec(wfs)});

        print_vec_pair_csv(ana_outfile, feature_value);
    }
  }
}
