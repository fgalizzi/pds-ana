#include "Utils.hpp"

using namespace std;
using json = nlohmann::json;

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void VGainScans_ana(cla& a, string jsonfile_module_config){
  // --- ANA CONFIG -----------------------------------------------
  AnaConfig ana_config = load_ana_config("config/ana_config.json");
  string runs_folder = ana_config.runs_folder;
  string input_ana_folder = ana_config.input_ana_folder;
  string rms_result_file = input_ana_folder+ana_config.rms_result_file;
  string output_ana_folder = ana_config.output_ana_folder;
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

  int module = module_config.module;
  vector<int> module_channels = module_config.module_channels;
  vector<double> v_brs = module_config.v_brs;
  vector<double> err_v_brs = module_config.err_v_brs;
  vector<int> vgains = module_config.vgains;
  vector<double> bias_dacs = module_config.bias_dacs;
  vector<vector<int>> run_batches = module_config.run_batches;
 
  // --- ANALYSIS -------------------------------------------------
  // Loop over the biases and analise the corresponding run batches
  for(size_t idx_bias=0; idx_bias<bias_dacs.size(); idx_bias++){
    double bias_dac = bias_dacs[idx_bias]; 
    std::vector<int> runs = run_batches[idx_bias];
    vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 

    string out_files_name = Form("Module_%i_Bias_%i", module, int(bias_dac));
    string out_root_file = output_ana_folder+out_files_name+".root";
    string out_csv_file  = output_ana_folder+out_files_name+".csv";

    std::vector<TString> files = {};
    for(auto& run : runs){
      files.push_back(runs_folder+"run_"+run);
    }

    TFile hf(TString(out_root_file), "recreate");
    hf.mkdir("chargehistos");
    hf.cd("chargehistos");
    vector<TH1D*> h_charge_vec;

    vector<pair<string, vector<double>>> ch_rms = read_vec_pair_CSV(rms_result_file.c_str()); // Store the results of the analysis to be printed 
   
    // --- LOOP OVER THE FILES == RUNS ----------------------------
    for(size_t idx_file=0; idx_file<files.size(); idx_file++){
      // --- LOOP OVER THE CHANNELS --------------------------------
      // We have two channels per module
      for(size_t idx_channel=0; idx_channel<module_channels.size(); idx_channel++){
        double bias_volt, overvoltage, err_bias_volt, err_overvoltage;
        give_me_Bias_OV_and_errors(module, bias_dac, v_brs[idx_channel], err_v_brs[idx_channel],
                                   bias_volt, overvoltage, err_bias_volt, err_overvoltage);
       
        // Look for baseline RMS in the file with the RMS results
        for(size_t j=0; j<ch_rms[0].second.size(); j++){
          if(ch_rms[1].second[j] == vgains[idx_file] && int(ch_rms[2].second[j]) == module_channels[idx_channel]){
            a.bsl = allowed_bsl_rms*ch_rms[3].second[j];
          }
        }
       
        // Hard Code: we don't have noise runs with vgain==0
        if (a.bsl>1000){
          a.bsl = 200;
        }
        a.sat_up = a.bsl*10;

        // Load the file and check if it exists
        a.wf_file = files[idx_file]+"/channel_"+module_channels[idx_channel]+".dat";
        ifstream this_file(a.wf_file);
        if (!this_file.is_open()){
          std::cout << "File not found: " << a.wf_file << std::endl;
          this_file.close();
          continue;
        }
        
        // The actual analysis
        a.LED_Analysis();
        a.LoadFitParameters(a.fgaus);
        a.SPE();
      
        // --- OUTPUT ------------------------------------------------
        int daphne_channel = module_channels[idx_channel];
        int arap_ch;
        if (daphne_channel == 0 || daphne_channel == 1
            || daphne_channel == 27 || daphne_channel == 26) arap_ch = 1;
        else arap_ch = 2;

        a.h_charge->SetTitle(Form("Ch_%i_Bias_%.2f_VGain_%i", module_channels[idx_channel], bias_volt, vgains[idx_file]));
        a.h_charge->SetName(Form("Ch_%i_Bias_%.2f_VGain_%i", module_channels[idx_channel], bias_volt, vgains[idx_file]));
        h_charge_vec.push_back(a.h_charge);
        feature_value.push_back({"Run", double(runs[idx_file])});
        feature_value.push_back({"ARAPUCA Channel", double(arap_ch)});
        feature_value.push_back({"DAPHNE Channel", double(daphne_channel)});
        feature_value.push_back({"Bias [dac]", bias_dac});
        feature_value.push_back({"Bias [V]", bias_volt});
        feature_value.push_back({"Err Bias [V]", err_bias_volt});
        feature_value.push_back({"OV [V]", overvoltage});
        feature_value.push_back({"Err OV [V]", err_overvoltage});
        feature_value.push_back({"VGain", double(vgains[idx_file])});
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
}
