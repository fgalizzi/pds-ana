#include "../../../Class/_c/class_include.hpp"
using namespace std;

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void NoiseVGain_ana(){
  // --------------------------------------------------------------
  // --- HARD CODE ------------------------------------------------
  string base_path   = "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/Daphne_DAQ/binaries/";
  string output_ana_folder  = "/eos/home-f/fegalizz/ColdBox_VD/December24/Daphne_DAQ/Noise_RMS_FFTs/";
  int LED = 0; // 0 or 1
  // Class settings
  auto a = cla();
  a.data_format = "esteban";
  a.n_wf = 30000;
  a.memorydepth = 1024;
  a.prepulse_ticks = 1023;
  a.tick_len = 0.016;
  a.res = 14;
  a.print = 1;

  // --- CODE ----------------------------------------------------
  string fft_outfile, ana_outfile;
  fft_outfile = output_ana_folder+Form("FFT_VGainScans_LED%i_Membrane.root", LED);
  ana_outfile = output_ana_folder+Form("VGain_RMS_LED%i_Membrane.csv", LED);
  
  vector<int> channels = {0,1,2,3,20,21,26,27};

  vector<int> runs;
  for(int i=9; i<300; i++){
    if((i+(LED+1))%2==0) runs.push_back(33300+i);
    if(i==66) i=10000;
  }

  int vgain = 100;

  // Tuple input file, run, vgain, channel
  vector<tuple<string,int,int,int>> my_tuple;
  vector<string> ifiles;
  for(auto& run : runs){
    for(auto& ch : channels){
      string ifile = base_path+"run_"+to_string(run)+"/channel_"+to_string(ch)+".dat";
      my_tuple.push_back(make_tuple(ifile,run,vgain,ch));
    }
    vgain += 100;
  }

  TFile fft_ofile(fft_outfile.c_str(), "recreate");
  for(auto& ch : channels){
    string dir_name = "Ch_"+to_string(ch);
    fft_ofile.mkdir(dir_name.c_str());
  }
  std::cout << "dd" << std::endl;
  // --- LOOP OVER FILES -----------------------------------------
  for(auto& tuple : my_tuple){
    a.wf_file = get<0>(tuple);
    a.read();
    TGraph* gNoise_spectral_density = build_avg_spectral_density(a.memorydepth,
      a.tick_len*a.memorydepth, a.tick_len, a.wfs, a.res);

    string gr_name = "VGain_"+to_string(get<2>(tuple));
    gNoise_spectral_density->SetName(gr_name.c_str());
    gNoise_spectral_density->SetTitle(gr_name.c_str());
    
    string dir_name = "Ch_"+to_string(get<3>(tuple));
    fft_ofile.cd(dir_name.c_str());
    gNoise_spectral_density->Write();
    if (a.print == true){
        vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 

        feature_value.push_back({"Run", get<1>(tuple)});
        feature_value.push_back({"VGain", get<2>(tuple)});
        feature_value.push_back({"Ch", get<3>(tuple)});
        feature_value.push_back({"Rms", standard_deviation_vec_vec(a.wfs)});

        print_vec_pair_csv(ana_outfile, feature_value);
    }
  }
}
