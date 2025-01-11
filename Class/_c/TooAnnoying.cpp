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

void cla::TooAnnoying(){

//Scan tick
int a = int_up;
for(int i = a; i < a+30; i++){
  int_up = i;
  LED_Analysis();
}
  





/*  TString runs_folder = "/eos/home-f/fegalizz/ProtoDUNE_HD/VGain/files/run_";
  TString ana_folder  = "/eos/home-f/fegalizz/ProtoDUNE_HD/VGain/analysis";
  //std::vector<int> runs = {29947,30053,30101,30149,30263,30314,30377,30458,30539,30620,30703,
    //                        30783,30945,31187,31166};

  print=1;
  Noise_PSD();
  trg_f="";
}

/*
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
  TString runs_folder = "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/Daphne_DAQ/";
  TString ana_folder  = "/eos/home-f/fegalizz/ColdBox_VD/December24/Daphne_DAQ/";
  std::vector<int> runs = {};
  std::vector<int> vgains = {};
  for(int i=33740; i<=33759; i++){
    runs.push_back(i);
    vgains.push_back((i-33740)*100+100);
  }
  // std::vector<int> runs = {30620};
  std::vector<TString> files = {};
  vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 
  for(auto& run : runs){
    files.push_back(runs_folder+"run_"+run+"/104/");
  }
  vector<int> channel_this_mask = {0, 2};
  // vector<int> channel_this_mask = {01, 02, 03, 04, 06, 10, 11, 12, 13, 14, 15, 17, 22,
  //                 24, 25, 26, 30, 31, 32, 33, 34, 35, 36, 37, 41, 42, 43, 44, 45, 46, 47};

  string path(ana_folder.Data());

  TFile hf("attemptroot.root", "recreate");
  hf.mkdir("chargehistos");
  hf.cd("chargehistos");
  vector<TH1D*> h_charge_vec;

  string rms_result_file((ana_folder+"VGain_RMS_LED0_Membrane.csv").Data());
  vector<pair<string, vector<double>>> ch_rms = read_vec_pair_CSV(rms_result_file.c_str()); // Store the results of the analysis to be printed 
  print = 1;
  display=0;
  
  for(size_t i=0; i<files.size(); i++){
    for(auto&ch : channel_this_mask){
      n_wf = 11000;
      for(size_t j=0; j<ch_rms[0].second.size(); j++){
        if(ch_rms[1].second[j] == vgains[i] && int(ch_rms[2].second[j]) == ch){
          bsl = 4.*ch_rms[3].second[j];
          sat_up = bsl*10;
        }
      }
      wf_file = files[i]+"channel_"+ch+".dat";
      ifstream this_file(wf_file);
      if (!this_file.is_open()){
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
      feature_value.push_back({"VGain", double(vgains[i])});
      feature_value.push_back({"Spe ampl", spe_ampl});
      feature_value.push_back({"DR", pow(2,14)/spe_ampl});
      feature_value.push_back({"SNR", spe_charge/sigma_zero});
      print = 1;
      if(print==true) std::cout << "\n\nPRINTING\n\n" << std::endl;
      if(print==true) print_vec_pair_csv("attempts.csv", feature_value);
      // if(print==true) print_vec_pair_csv(Form("%s/VGain_results_Ep_104_run_%d.csv",path.c_str(),runs[0]), feature_value);
      feature_value = {};
      // cout << Form("%s/VGain_results_Ep_104_run_%d.csv",path.c_str(),runs[0]) << endl;
    }
  }
  
  //hf.cd("chargehistos"); for(auto h : h_charge_vec) h->Write();*/

 
  std::cout << "\n\nOUT OF THE LOOP\n\n" << std::endl;
  hf.cd("chargehistos"); for(auto h : h_charge_vec) h->Write();
  hf.Close();


}
*/
