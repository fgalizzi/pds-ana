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

///////////////////////////////////////////////////////////////////



std::vector<std::pair<std::string, std::vector<double>>> read_vec_pair_CSV(const std::string& filename) {
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  std::vector<std::pair<std::string, std::vector<double>>> data;
  std::string line;

  // Read the header line
  if (std::getline(infile, line)) {
    std::stringstream ss(line);
    std::string header;
    while (std::getline(ss, header, ',')) {
      if (!header.empty()) {
        data.emplace_back(header, std::vector<double>{});
      }
    }
  }

  // Read the data lines
  while (std::getline(infile, line)) {
    std::stringstream ss(line);
    std::string value;
    size_t col_index = 0;

    while (std::getline(ss, value, ',') && col_index < data.size()) {
      try {
        if (!value.empty()) {
          data[col_index].second.push_back(std::stod(value)); // Convert to double
        }
      } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid number at column " << col_index + 1 << std::endl;
        data[col_index].second.push_back(0.0); // Default value for invalid numbers
      }
      col_index++;
    }
  }

  infile.close();
  return data;
}



//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::TooAnnoying(){
  // -------------------------------------------------------------
  // --- HARD CODE -----------------------------------------------
  // INPUT
  TString runs_folder        = "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/Daphne_DAQ/binaries/";
  TString input_ana_folder   = "/eos/home-f/fegalizz/ColdBox_VD/December24/Daphne_DAQ/Noise_RMS_FFTs/";
  TString output_ana_folder  = "/eos/home-f/fegalizz/ColdBox_VD/December24/Daphne_DAQ/VGain_Scans/";
  // Runs and corresponding VGains
  std::vector<int> runs = {};
  std::vector<int> vgains = {};
  for(int i=33740; i<=33759; i++){
    runs.push_back(i);
    vgains.push_back((i-33740)*100+100);
  }
  // Channels good for these runs
  vector<int> channel_this_mask = {0, 2};
  double bias_dac   = 793;
  double bias_volts = 32;
  // File with the RMS of the channels
  string rms_result_file((input_ana_folder+"VGain_RMS_LED0_Membrane.csv").Data());

  // CLASS SETTINGS
  n_wf = 5000;
  display=0;
  print = 0;

  // OUTPUT
  bool print_results = true;
  TString out_root_file = output_ana_folder+"VGain_Scan_Module_Bias.root";
  string out_csv_file((output_ana_folder+"VGain_Scan_Module_Bias.csv").Data());
  // --- END HARD CODE -------------------------------------------
 

  // --- CODE ----------------------------------------------------
  vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 
  std::vector<TString> files = {};
  for(auto& run : runs){
    files.push_back(runs_folder+"run_"+run);
  }

  // string path(ana_folder.Data());

  TFile hf(out_root_file, "recreate");
  hf.mkdir("chargehistos");
  hf.cd("chargehistos");
  vector<TH1D*> h_charge_vec;

  vector<pair<string, vector<double>>> ch_rms = read_vec_pair_CSV(rms_result_file.c_str()); // Store the results of the analysis to be printed 
 
  std::cout << "files " << files.size() << std::endl;
  for(size_t i=0; i<files.size(); i++){
    for(auto&ch : channel_this_mask){
      for(size_t j=0; j<ch_rms[0].second.size(); j++){
        if(ch_rms[1].second[j] == vgains[i] && int(ch_rms[2].second[j]) == ch){
          bsl = 4.*ch_rms[3].second[j];
          sat_up = bsl*10;
        }
      }
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
      LED_Analysis();
      LoadFitParameters(fgaus);
      SPE();

      h_charge->SetTitle(Form("VGain_%i_ch_%i",vgains[i],ch));
      h_charge->SetName(Form("VGain_%i_ch_%i",vgains[i],ch));
      h_charge_vec.push_back(h_charge);
      feature_value.push_back({"Channel", double(ch)});
      feature_value.push_back({"Bias [dac]", bias_dac});
      feature_value.push_back({"Bias [V]", bias_volts});
      feature_value.push_back({"VGain", double(vgains[i])});
      feature_value.push_back({"Baseline", bsl});
      feature_value.push_back({"Prepulse ticks", double(prepulse_ticks)});
      feature_value.push_back({"Saturation up", sat_up});
      feature_value.push_back({"Int low", double(int_low)});
      feature_value.push_back({"Int up", double(int_up)});
      feature_value.push_back({"Gain", spe_charge});
      feature_value.push_back({"Spe ampl", spe_ampl});
      feature_value.push_back({"DR", pow(2,14)/spe_ampl});
      feature_value.push_back({"SNR", spe_charge/sigma_zero});
      feature_value.push_back({"RMS", ch_rms[3].second[i]});
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
      // if(print==true) print_vec_pair_csv(Form("%s/VGain_results_Ep_104_run_%d.csv",path.c_str(),runs[0]), feature_value);
      // cout << Form("%s/VGain_results_Ep_104_run_%d.csv",path.c_str(),runs[0]) << endl;
      
      // Reset the vector
      feature_value = {};
    }
  }
  
  std::cout << "\n\nOUT OF THE LOOP\n\n" << std::endl;
  if(print_results==true) hf.cd("chargehistos"); for(auto h : h_charge_vec) h->Write();
  hf.Close();
}
