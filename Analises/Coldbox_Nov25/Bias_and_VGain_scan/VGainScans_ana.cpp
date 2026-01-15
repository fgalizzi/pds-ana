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
  // --- ANA CONFIG -----------------------------------------------
  AnaConfig ana_config     = load_ana_config("config/ana_config.json");
  string runs_folder       = ana_config.runs_folder;
  string input_ana_folder  = ana_config.input_ana_folder;
  string rms_result_file   = input_ana_folder+ana_config.rms_result_file;
  string output_ana_folder = ana_config.output_ana_folder;
  double allowed_bsl_rms   = ana_config.allowed_bsl_rms;
  bool   print_results     = ana_config.print_results;
  string bias_fit_csv      = ana_config.bias_fit_csv;
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
  vector<int> module_channels  = module_config.module_channels;
  vector<double> biases        = module_config.daphne_biases;
  vector<int> vgains           = module_config.vgains;
  string modules_in_foldername = module_config.modules_in_foldername;
  string led_afe_extension     = module_config.led_afe_extension;
  string custom_vgain_folder   = module_config.custom_vgain_folder;

  auto bias_fit_table = read_bias_fit_csv(bias_fit_csv);
 
  // --- ANALYSIS -------------------------------------------------
  // Loop over the biases and analise the corresponding run batches
  for(auto& bias : biases){
    string module_name = Form("M%i", module);
    string sub_folder = "VGAIN_SCAN/";
    if (modules_in_foldername != ""){
      module_name = modules_in_foldername;
      sub_folder = "debugs/";
    }

    ostringstream oss;
    oss << int(bias)
        << "V"
        << setw(2) << setfill('0') << static_cast<int>(std::round(bias * 100)) % 100;
    string bias_str = oss.str();
    string bias_folder = runs_folder+sub_folder+module_name+"/vgain_scan_"+module_name+"_DVbias_"+bias_str+led_afe_extension+"/";
    if (custom_vgain_folder != ""){
      bias_folder = custom_vgain_folder;
    }
    vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 

    string out_files_name = string(Form("M%i",module))+"/Bias_"+bias_str+led_afe_extension;
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
        int AFE = channel/8;
        double real_bias_value;
        int ch_index = distance(module_channels.begin(),
                                find(module_channels.begin(), module_channels.end(), channel));
        double breakdown_voltage;
        if (ch_index % 2 == 0){
          breakdown_voltage = breakdown_voltages[module].first;
        } else {
          breakdown_voltage = breakdown_voltages[module].second;
        }

        if (module < 5){
          real_bias_value = real_bias(AFE, "DaphnetoMult", bias,  bias_fit_table);
        } else {
          real_bias_value = breakdown_voltage + 5.;
        }
        hf.cd(Form("Ch_%i", channel));
        // Look for baseline RMS in the file with the RMS results
        // We take channel 17 as reference
        for(size_t j=0; j<ch_rms[0].second.size(); j++){
          if( ch_rms[0].second[j] == vgain && int(ch_rms[2].second[j]) == 17){
            a.bsl = allowed_bsl_rms*ch_rms[5].second[j];
            if (module == 4) a.bsl *= 1.5;
            break;
          }
        }
        cout << a.bsl << endl;
        
        a.sat_up = a.bsl*20;
        a.sat_low = -a.sat_up/2.;

        // Load the file and check if it exists
        a.wf_file = bias_folder+"vgain_"+vgain+"/channel_"+channel+".dat";
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
        a.h_charge->SetTitle(Form("Bias_%s_VGain_%i", bias_str.c_str(), vgain));
        a.h_charge->SetName(Form("Bias_%s_VGain_%i",  bias_str.c_str(), vgain));
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
