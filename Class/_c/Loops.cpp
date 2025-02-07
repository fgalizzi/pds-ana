#include "../classe.hpp"

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

// ****************************************************************
// Loop to fit a calibration spectrum for different VBias.
// Useful to builv the Gain vs VBias plot and estrapolate the
// SiPM breakdown voltage and to build the cross-talk vs VBias plot.
// Suggestion: set "scan_sat_up" according to the highest VBias
// set in the scan.
// ****************************************************************
//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Loop_VBias_Scan(){
  // -------------------------------------------------------------
  // --- HARD CODE -----------------------------------------------
  // INPUT
  TString runs_folder        = "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/Daphne_DAQ/binaries/";
  TString output_ana_folder  = "/eos/home-g/gpiemont/ColdBox_VD/December24/Daphne_DAQ/FineBiasScan/";
  // Runs and corresponding Bias
  // Channels good for these runs
std::vector<double> biases = {751,758,765,772,779,786,793,800,807,814,821,828};
std::vector<int> runs = {34353,34352,34351,34350,34349,34348,34347,34346,34345,34344,34343,34342};
  int module = 3; // M1 (20,27), M2 (21,26), M3 (0,2), M4 (1,3)
  vector<int> channel_this_mask = {0, 2};
  // Initial sat_up for the scan
  double scan_sat_up = 1600;

  // CLASS SETTINGS
  display = 0;
  print   = 0;
  plot    = 0;

  // OUTPUT
  bool print_results = true;
  TString out_file = output_ana_folder+Form("VBias_Scan_Module_%i", module);
  TString out_root_file = out_file+".root";
  string out_csv_file(out_file+".csv");
  // --- END HARD CODE -------------------------------------------
 

  // --- CODE ----------------------------------------------------
  if(biases.size() != runs.size()){
    std::cout << "Biases and runs vectors must have the same size" << std::endl;
    return;
  }

  vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 
  std::vector<TString> files = {};
  for(auto& run : runs){
    files.push_back(runs_folder+"run_"+run);
  }

  TFile hf(out_root_file, "recreate");
  hf.mkdir("chargehistos");
  hf.cd("chargehistos");
  vector<TH1D*> h_charge_vec;

 
  std::cout << "files " << files.size() << std::endl;
  for(size_t i=0; i<files.size(); i++){
    for(auto&ch : channel_this_mask){
      wf_file = files[i]+"/channel_"+ch+".dat";
      std::cout << wf_file << std::endl;
      ifstream this_file(wf_file);
      if (!this_file.is_open()){
        std::cout << "File not found: " << wf_file << std::endl;
        this_file.close();
        continue;
      }
      std::cout << "\n\n\nReading file: " << wf_file << std::endl;  
      cout << wf_file << endl;
      sat_up = scan_sat_up;
      LED_Analysis();
      LoadFitParameters(fgaus);
      SPE();
      sat_up = spe_ampl*20;
      LED_Analysis();
      LoadFitParameters(fgaus);
      SPE();

      h_charge->SetTitle(Form("VBias_%i_ch_%i",int(biases[i]),ch));
      h_charge->SetName(Form("VBias_%i_ch_%i",int(biases[i]),ch));
      h_charge_vec.push_back(h_charge);

      int arap_ch;
      if (ch == 0 || ch == 1 || ch == 27 || ch == 26) arap_ch = 1;
      else arap_ch = 2; 
      
      feature_value.push_back({"Run", runs[i]});
      feature_value.push_back({"Membrane modules Channel", int(arap_ch)});
      feature_value.push_back({"Channel", double(ch)});
      feature_value.push_back({"Bias [dac]", biases[i]});
      // feature_value.push_back({"Bias [V]", bias_volts});
      // feature_value.push_back({"VGain", double(vgains[i])});
      feature_value.push_back({"Baseline", bsl});
      feature_value.push_back({"Prepulse ticks", double(prepulse_ticks)});
      feature_value.push_back({"Saturation up", sat_up});
      feature_value.push_back({"Int low", double(int_low)});
      feature_value.push_back({"Int up", double(int_up)});
      feature_value.push_back({"Gain", spe_charge});
      feature_value.push_back({"Err Gain", err_spe_charge});
      feature_value.push_back({"Spe ampl", spe_ampl});
      feature_value.push_back({"DR", pow(2,14)/spe_ampl});
      feature_value.push_back({"SNR", SNR});
      feature_value.push_back({"Err SNR", err_SNR});
      feature_value.push_back({"CX", cx});
      feature_value.push_back({"Err CX", err_cx});
      feature_value.push_back({"Avg #ph cx", avg_n_ph_cx});
      feature_value.push_back({"Err #ph cx", err_avg_n_ph_cx});
      feature_value.push_back({"Avg #ph", avg_n_photons});
      feature_value.push_back({"Avg #pe", avg_n_photoelectrons});
      
      if(print_results==true){
        std::cout << "\n\nPRINTING\n\n" << std::endl;
        print_vec_pair_csv(out_csv_file, feature_value);
      }
      
      // Reset the vector
      feature_value = {};
    }
  }
  
  std::cout << "\n\nOUT OF THE LOOP\n\n" << std::endl;
  if(print_results==true) hf.cd("chargehistos"); for(auto h : h_charge_vec) h->Write();
  hf.Close();
}


      
