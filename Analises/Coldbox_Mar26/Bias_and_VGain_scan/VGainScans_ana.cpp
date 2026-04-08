#include "TString.h"
#include "Utils.hpp"
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;
using json = nlohmann::json;

// map in the form { Module : breakdown voltage }
// Values taken from the spreasdheet adding 0.2 V to HPK
// modules and 0.1 V to FBK modules to account LAr-LN2
// temperature difference, when needed
map<int, pair<double,double>> breakdown_voltages = {
  {1, {41.4, 41.5}},
  {2, {41.5, 41.6}},
  {3, {26.9, 26.9}},
  {4, {26.8, 26.8}},
  {5, {42.7, 42.6}},
  {6, {42.0, 42.1}},
};

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void VGainScans_ana(cla& a, string jsonfile_module_config){
  gROOT->SetBatch(kTRUE);
  // --- ANA CONFIG -----------------------------------------------
  AnaConfig ana_config     = load_ana_config("config/ana_config.json");
  string runs_folder       = ana_config.runs_folder;
  string input_ana_folder  = ana_config.input_ana_folder;
  string output_ana_folder = ana_config.output_ana_folder;
  double allowed_bsl_rms   = ana_config.allowed_bsl_rms;
  bool   print_results     = ana_config.print_results;
  // Class settings 
  a.display     = ana_config.display;
  a.plot        = ana_config.plot;
  a.print       = ana_config.print;
  a.tick_len    = ana_config.tick_len;
  a.nbins       = ana_config.nbins;
  a.data_format = ana_config.data_format;

  
  // --- MODULE CONFIG --------------------------------------------
  ModuleConfig module_config = load_module_config(jsonfile_module_config);
  // Use configuration from JSON:
  a.memorydepth    = module_config.memorydepth;
  a.invert         = module_config.invert;
  a.n_wf           = module_config.n_wf;
  a.prepulse_ticks = module_config.prepulse_ticks;
  a.int_low        = module_config.int_low;
  a.int_up         = module_config.int_up;
  a.nmaxpeaks      = module_config.nmaxpeaks;

  int module                   = module_config.module;
  double bias_slope            = module_config.bias_slope;
  double bias_offset           = module_config.bias_offset;
  vector<int> module_channels  = module_config.module_channels;
  vector<int> biases           = module_config.daphne_biases;
  vector<int> vgains           = module_config.vgains;
  string folder_extension      = module_config.folder_extension;
  vector<string> led_folders   = module_config.led_folders;
  string rms_result_file       = input_ana_folder+module_config.rms_result_file;

  // --- ANALYSIS -------------------------------------------------
  // Loop over the biases and analise the corresponding run batches
  for(auto& bias : biases){
    string electronics = (module <= 2) ? "VD" : "HD";
    int ref_channel = -1;
    if (module > 4) {
      electronics = "PoF";
      ref_channel = 12;
    }
    string module_name = Form("M%i", module);
    string sub_folder = "/BiasVGainScan"+folder_extension;
  
    
    int AFE = module_channels[0]/8;

    // ostringstream oss;
    // oss << int(bias)
        // << "V"
        // << setw(2) << setfill('0') << static_cast<int>(std::round(bias * 100)) % 100;
    // string bias_str = oss.str();
    string zero = (bias < 1000) ? "0" : "";
    if (bias == 0) zero = "000";
    string bias_folder = runs_folder+electronics+sub_folder+"/afe"+to_string(AFE)+"bias_"+zero+to_string(bias)+"/";
    vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 

    string out_files_name = string(Form("M%i",module))+"/Bias_"+to_string(bias)+folder_extension;
    string out_root_file  = output_ana_folder+out_files_name+".root";
    string out_csv_file   = output_ana_folder+out_files_name+".csv";
    if (std::remove(out_csv_file.c_str()) == 0) {
      cout << "Old output file removed: " << out_csv_file << endl;
    }

    TFile hf(TString(out_root_file), "recreate");
    for(auto& channel : module_channels){
      hf.mkdir(Form("Ch_%i", channel));
    }

    vector<pair<string, vector<double>>> ch_rms = read_vec_pair_CSV(rms_result_file.c_str()); // Store the results of the analysis to be printed 
   
    // --- LOOP OVER THE VGAINS ------------------------------------
    for(auto& vgain : vgains){
      // --- LOOP OVER THE CHANNELS --------------------------------
      // We have two channels per module
      for(auto& channel : module_channels){
        double real_bias_value = bias*bias_slope + bias_offset;
        int ch_index = distance(module_channels.begin(),
                                find(module_channels.begin(), module_channels.end(), channel));
        double breakdown_voltage;
        if (ch_index % 2 == 0){
          breakdown_voltage = breakdown_voltages[module].first;
        } else {
          breakdown_voltage = breakdown_voltages[module].second;
        }
        
        hf.cd(Form("Ch_%i", channel));
        // Look for baseline RMS in the file with the RMS results
        for(size_t j=0; j<ch_rms[0].second.size(); j++){
          if( ch_rms[0].second[j] == vgain && ( int(ch_rms[2].second[j]) == channel || (int(ch_rms[2].second[j]) == ref_channel) ) ){
            a.bsl = allowed_bsl_rms*ch_rms[5].second[j];
            if (vgain == 1500 && (electronics == "VD")) a.bsl += 7.;// HARD CODE
                                                                    // Fix underestimated RMS for VD at VGain=1500
            break;
          }
        }
        cout << "Found RMS = " << a.bsl << endl;
       
        if (electronics == "HD" && bias == 874) {
          a.statlost_low = 0.002; a.statlost_up = 0.90;
        }
        else {
          a.statlost_low = 0.005; a.statlost_up = 0.99;
        }
        a.sat_up = a.bsl*20;
        a.sat_low = -a.sat_up/2.;

        // Load the file and check if it exists
        string o = (vgain < 1000) ? "0" : "";
        string led_folder = (led_folders.size() > 0) ? "/"+led_folders[ch_index] : "";
        a.wf_file = bias_folder+"vgain_"+o+to_string(vgain)+led_folder+"/channel_"+channel+".dat";
        ifstream this_file(a.wf_file);
        if (!this_file.is_open()){
          cout << "File not found: " << a.wf_file << endl;
          this_file.close();
          continue;
        }
        
        // The actual analysis
        cout << "Start LED analysis..." << endl;
        a.LED_Analysis();
        if (a.class_skip == 1){
          cout << "Skipping channel " << channel << " at VGain " << vgain << " due to lack of calibration wfs" << endl;
          a.class_skip = 0;
          continue;
        }
        cout << "...end" << endl;
        a.LoadFitParameters(a.fgaus);
        a.h_charge->SetTitle(Form("Bias_%i_VGain_%i", bias, vgain));
        a.h_charge->SetName(Form("Bias_%i_VGain_%i",  bias, vgain));
        if (a.SNR < 2 || a.SNR > 20) continue;

        cout << "Start SPE..." << endl;
        a.SPE();
        if (a.class_skip == 1){
          cout << "Skipping channel " << channel << " at VGain " << vgain << " due to lack of SPE candidates" << endl;
          a.class_skip = 0;
          continue;
        }
        cout << "...end" << endl;
      
        // --- OUTPUT ------------------------------------------------
        feature_value.push_back({"Module", double(module)});
        feature_value.push_back({"DAPHNE Channel", double(channel)});
        feature_value.push_back({"Bias [V]", bias});
        feature_value.push_back({"Real Bias [V]", real_bias_value});
        feature_value.push_back({"OV [V]", real_bias_value - breakdown_voltage});
        feature_value.push_back({"VGain", vgain});
        feature_value.push_back({"Baseline", a.bsl});
        feature_value.push_back({"Prepulse ticks", double(a.prepulse_ticks)});
        feature_value.push_back({"Saturation up", a.sat_up});
        feature_value.push_back({"Int low", double(a.int_low)});
        feature_value.push_back({"Int up", double(a.int_up)});
        feature_value.push_back({"Gain", a.spe_charge});
        feature_value.push_back({"Err Gain", a.err_spe_charge});
        feature_value.push_back({"Spe ampl", a.spe_ampl});
        feature_value.push_back({"DR", pow(2,14)/a.spe_ampl});
        feature_value.push_back({"SNR", a.SNR});
        feature_value.push_back({"Err SNR", a.err_SNR});
        feature_value.push_back({"RMS", a.bsl/allowed_bsl_rms});
        feature_value.push_back({"CX", a.cx});
        feature_value.push_back({"Err CX", a.err_cx});
        feature_value.push_back({"Avg #ph cx", a.avg_n_ph_cx});
        feature_value.push_back({"Err #ph cx", a.err_avg_n_ph_cx});
        feature_value.push_back({"Avg #ph", a.avg_n_photons});
        feature_value.push_back({"Avg #pe", a.avg_n_photoelectrons});
        feature_value.push_back({"StatLoss [%]", double(1-(a.h_charge->GetEntries()/a.n_wf))*100.});
        feature_value.push_back({"CalibrationWFS", double(a.h_charge->GetEntries())});
        feature_value.push_back({"TotalWFS", double(a.n_wf)});
        
        if(print_results==true){
          cout << "\n\nPRINTING\n\n" << endl;
          print_vec_pair_csv(out_csv_file, feature_value);
          a.h_charge->Write(); a.h_charge->Delete();
        }
        
        // Reset the vector
        feature_value = {};
      }
    }
    
    cout << "\n\nOUT OF THE LOOP\n\n" << endl;
    hf.Close();
  }
}
