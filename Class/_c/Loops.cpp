#include "../classe.hpp"
#include <cstddef>
#include <utility>

// ****************************************************************
// Description
// ****************************************************************

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Loop_ST_Analysis(){
  //////////////////////////////////////////////////////////////////////
  // HARD CODE HERE
  //
  calibration_run = 9999; //Needed create the root file with the results
  ep = 112;
  mask = 1400;
  channel_low = ep*100;     //Lower channel to look at (included)
  channel_up  = (ep+1)*100;     //Upper " "
  //////////////////////////////////////////////////////////////////////

  vector<string> signal_files = read_chs("ch_sipm.txt");
  vector<string> self_files   = read_chs("ch_self.txt");
  vector<size_t> channels = read_pdhd_ch_map(mask);
  print = 1;
  
  for(size_t i=0; i<signal_files.size(); i++){
    wf_file = signal_files[i];
    trg_f = self_files[i];
    //Check wether the LED mask was good for this channel
    channel = size_t(extract_channel_from_filename(wf_file))+ep*100;
    int cnt = std::count(channels.begin(), channels.end(), channel);
    std::cout << cnt << std::endl;
    if (cnt == 0 || channel < channel_low || channel > channel_up) continue;
    std::cout << i << std::endl;
    
    //Calibrate and perform the self-trigger analisys
    LED_Analysis();
    LoadFitParameters(fgaus);
    ST_Analysis();
  }

}

// ****************************************************************
// Loop to compute the FFT and RMS of noise runs. The FFTs are
// stored as TGraph in a fft_outfile.root, the RMS in an
// ana_outfile.csv (see my_tuple)
// ****************************************************************
//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Loop_FFT_RMS_Analysis(){
///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////
  string base_path   = "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/Daphne_DAQ/";
  string fft_outfile = "FFT_VGainScans_LED1_Membrane.root";
  string ana_outfile = "VGain_RMS_LED1_Membrane.csv";
///////////////////////////////////////////////////////////////////
  

  // Tuple input file, run, vgain, channel
  vector<tuple<string,int,int,int>> my_tuple;

  vector<int> runs;
  for(int i=9; i<300+9; i++){
    if(i%2==0) runs.push_back(33300+i);
    if(i==66) i=10000;
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
    
    string dir_name = "Ch_"+to_string(get<3>(tuple));
    fft_ofile.cd(dir_name.c_str());
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

