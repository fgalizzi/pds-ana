#include "Utils.hpp"

using namespace std;
using json = nlohmann::json;

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void FineBiasScan_ana(cla& a, string jsonfile_module_config){
  // --- ANA CONFIG -----------------------------------------------
  AnaConfig ana_config = load_ana_config("config/ana_config.json");
  TString runs_folder = ana_config.runs_folder;
  string input_ana_folder = ana_config.input_ana_folder;
  string rms_result_file = input_ana_folder+ana_config.rms_result_file;
  TString output_ana_folder = ana_config.output_ana_folder;
  double allowed_bsl_rms = ana_config.allowed_bsl_rms;
  bool   print_results = ana_config.print_results;
  // Class settings 
  a.display = ana_config.display;
  a.print= ana_config.print;
  a.plot = ana_config.plot;
  a.memorydepth = ana_config.memorydepth;
  a.tick_len = ana_config.tick_len;
  a.nbins = ana_config.nbins;
  a.nmaxpeaks = ana_config.nmaxpeaks;
  a.data_format = ana_config.data_format;
  
  // --- MODULE CONFIG --------------------------------------------
  ModuleConfig module_config = load_module_config(jsonfile_module_config);
  // Use configuration from JSON:
  a.invert = module_config.invert;
  a.n_wf = module_config.n_wf;
  a.prepulse_ticks = module_config.prepulse_ticks;
  a.int_low = module_config.int_low;
  a.int_up = module_config.int_up;
  double scan_sat_up = module_config.scan_sat_up; // Initial sat_up for the scan

  int module = module_config.module;
  vector<int> module_channels = module_config.module_channels;
  vector<double> biases = module_config.biases; 
  vector<int> runs = module_config.runs;

  if (biases.size() != runs.size()){
    std::cout << "Biases and runs vectors must have the same size" << std::endl;
    return;
  }

  // --- HARD CODE -----------------------------------------------
  // OUTPUT
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

  vector<pair<string, vector<double>>> ch_rms = read_vec_pair_CSV(rms_result_file.c_str()); // Store the results of the analysis to be printed 

  // --- LOOP OVER FILES == RUNS -----------------------------------
  for(size_t i=0; i<files.size(); i++){
    // --- LOOP OVER CHANNELS --------------------------------------
    for(auto&ch : module_channels){
      a.wf_file = files[i]+"/channel_"+ch+".dat";
      ifstream this_file(a.wf_file);
      if (!this_file.is_open()){
        std::cout << "File not found: " << a.wf_file << std::endl;
        this_file.close();
        continue;
      }
      
      // Look for baseline RMS in the file with the RMS results
      for(size_t j=0; j<ch_rms[0].second.size(); j++){
        if(ch_rms[1].second[j] == 1000 && int(ch_rms[2].second[j]) == ch){
          a.bsl = allowed_bsl_rms*ch_rms[3].second[j];
          a.sat_up = a.bsl*10;
        }
      }

      // The actual analysis
      a.sat_up = scan_sat_up;
      a.sat_low = -0.1*a.sat_up;
      a.LED_Analysis();
      a.LoadFitParameters(a.fgaus);
      a.SPE();
      a.sat_up = a.spe_ampl*20;
      a.sat_low = -0.1*a.sat_up;
      a.LED_Analysis();
      a.LoadFitParameters(a.fgaus);
      a.SPE();

      // --- OUTPUT ------------------------------------------------
      a.h_charge->SetTitle(Form("VBias_%i_ch_%i",int(biases[i]),ch));
      a.h_charge->SetName(Form("VBias_%i_ch_%i",int(biases[i]),ch));
      h_charge_vec.push_back(a.h_charge);

      double bias_volts, err_bias_volts;
      DAC_to_Volt(module, biases[i], bias_volts, err_bias_volts);
      
      int arap_ch;
      if (ch == 0 || ch == 1 || ch == 27 || ch == 26) arap_ch = 1;
      else arap_ch = 2; 
      
      feature_value.push_back({"Run", runs[i]});
      feature_value.push_back({"ARAPUCA Channel", int(arap_ch)});
      feature_value.push_back({"DAPHNE Channel", double(ch)});
      feature_value.push_back({"Bias [dac]", biases[i]});
      feature_value.push_back({"Bias [V]", bias_volts});
      feature_value.push_back({"Err Bias [V]", err_bias_volts});
      feature_value.push_back({"VGain", double(1000)});
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
      feature_value.push_back({"CX", a.cx});
      feature_value.push_back({"Err CX", a.err_cx});
      feature_value.push_back({"Avg #ph cx", a.avg_n_ph_cx});
      feature_value.push_back({"Err #ph cx", a.err_avg_n_ph_cx});
      feature_value.push_back({"Avg #ph", a.avg_n_photons});
      feature_value.push_back({"Avg #pe", a.avg_n_photoelectrons});
      
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
