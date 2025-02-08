#include "../../../Class/_c/class_include.hpp"
#include <nlohmann/json.hpp>
using namespace std;
using json = nlohmann::json;

struct AnaConfig {
  string runs_folder;
  string input_ana_folder;
  string rms_result_file;
  string output_ana_folder;
  double allowed_bsl_rms;
  int display;
  int print;
  int plot;
  int memorydepth;
  double tick_len;
  int nbins;
  int nmaxpeaks;
  string data_format;
};

AnaConfig load_ana_config(const std::string& filename) {
  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error("Could not open config file: " + filename);
  }

  json j;
  file >> j;
  AnaConfig config;
  config.runs_folder = j.at("runs_folder").get<std::string>();
  config.input_ana_folder = j.at("input_ana_folder").get<std::string>();
  config.rms_result_file = j.at("rms_result_file").get<std::string>();
  config.output_ana_folder = j.at("output_ana_folder").get<std::string>();
  config.allowed_bsl_rms = j.at("allowed_bsl_rms").get<double>();
  config.display = j.at("display").get<int>();
  config.print = j.at("print").get<int>();
  config.plot = j.at("plot").get<int>();
  config.memorydepth = j.at("memorydepth").get<int>();
  config.tick_len = j.at("tick_len").get<double>();
  config.nbins = j.at("nbins").get<int>();
  config.nmaxpeaks = j.at("nmaxpeaks").get<int>();
  config.data_format = j.at("data_format").get<std::string>();
  return config;
}


struct ModuleConfig {
  int module;
  std::vector<int> module_channels;
  std::vector<double> v_brs;
  std::vector<double> err_v_brs;
  std::vector<int> vgains;
  std::vector<double> bias_dacs;
  std::vector<std::vector<int>> run_batches;
};

ModuleConfig load_module_config(const std::string& filename) {
  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error("Could not open config file: " + filename);
  }
  
  json j;
  file >> j;
  
  ModuleConfig config;
  config.module = j.at("module").get<int>();
  config.module_channels = j.at("module_channels").get<std::vector<int>>();
  config.v_brs = j.at("v_brs").get<std::vector<double>>();
  config.err_v_brs = j.at("err_v_brs").get<std::vector<double>>();
  config.vgains = j.at("vgains").get<std::vector<int>>();
  config.bias_dacs = j.at("bias_dacs").get<std::vector<double>>();
  config.run_batches = j.at("run_batches").get<std::vector<std::vector<int>>>();
  
  return config;
}



#include "config/module_1.hpp"

void give_me_Bias_OV_and_errors(int module, double bias_dac, double v_br,
                                double err_v_br, double& bias_volt,
                                double &overvoltage, double &err_bias_volt, double &err_overvoltage){

  double err_volt_m1_m2 = 0.03; // to have chi2 ~1
  double err_volt_m3_m4 = 0.07; // to have chi2 ~1
  vector<double> dacs, volts;  
  double err_volt;
  if(module==1 || module==2 || module==11){
    dacs = {1148, 1161, 1174, 1187, 1200};
    volts = {45.06, 45.54, 46.09, 46.55, 47.03};
    err_volt = err_volt_m1_m2;
  } else {
    dacs = {754, 767, 780, 793, 806};
    volts = {30.54, 31.11, 31.64, 32.01, 32.52};
    err_volt = err_volt_m3_m4;
  }
  vector<double> err_volts(volts.size(), err_volt);
  vector<double> err_dacs(dacs.size(), 0.);

  TF1* f1 = new TF1 ( "f1", "[0]+[1]*x"); // y = q + m*x

  f1->SetParameter(0, 0.);
  f1->SetParameter(1, (volts[volts.size()-1]-volts[0])/(dacs[dacs.size()-1]-dacs[0]));

  TGraphErrors* g_Volt_DAC = new TGraphErrors(dacs.size(), &dacs[0], &volts[0],
                                              &err_dacs[0], &err_volts[0]);
  
  TFitResultPtr r1 = g_Volt_DAC->Fit(f1 ,"S");
  double this_bias_dac[1] = {bias_dac};
  double err_this_bias_volt[1];
  r1->GetConfidenceIntervals(1, 1, 1, this_bias_dac, err_this_bias_volt, 0.683, false);

  error_propagation(bias_volt, err_bias_volt, v_br, err_v_br, "sub");
  bias_volt = f1->Eval(this_bias_dac[0]);
  err_bias_volt = err_this_bias_volt[0];
  overvoltage = bias_volt - v_br;
  err_overvoltage = error_propagation(bias_volt, err_bias_volt, v_br, err_v_br, "sub");

  return;
}



//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void VGainScans_ana(cla& a, string jsonfile_module_config){
  AnaConfig ana_config = load_ana_config("config/ana_config.json");
  string runs_folder = ana_config.runs_folder;
  string input_ana_folder = ana_config.input_ana_folder;
  string rms_result_file = ana_config.rms_result_file;
  string output_ana_folder = ana_config.output_ana_folder;
  double allowed_bsl_rms = ana_config.allowed_bsl_rms;
  // Class settings 
  a.display = ana_config.display;
  a.print= ana_config.print;
  a.plot = ana_config.plot;
  a.memorydepth = ana_config.memorydepth;
  a.tick_len = ana_config.tick_len;
  a.nbins = ana_config.nbins;
  a.nmaxpeaks = ana_config.nmaxpeaks;
  a.data_format = ana_config.data_format;

  
  // Load configuration from JSON file
  ModuleConfig module_config = load_module_config(jsonfile_module_config);
  // Use configuration from JSON:
  int module = module_config.module;
  vector<int> module_channels = module_config.module_channels;
  vector<double> v_brs = module_config.v_brs;
  vector<double> err_v_brs = module_config.err_v_brs;
  vector<int> vgains = module_config.vgains;
  vector<double> bias_dacs = module_config.bias_dacs;
  vector<vector<int>> run_batches = module_config.run_batches;
  
  // Loop over the biases and analise the corresponding run batches
  for(size_t idx_bias=0; idx_bias<bias_dacs.size(); idx_bias++){
    double bias_dac = bias_dacs[idx_bias]; 
    std::vector<int> runs = run_batches[idx_bias];
    vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 

    string out_files_name = Form("Module_%i_Bias_%i_VGain_3RMS", module, int(bias_dac));
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
   
    std::cout << "files " << files.size() << std::endl;
    for(size_t idx_file=0; idx_file<files.size(); idx_file++){
      for(size_t idx_channel=0; idx_channel<module_channels.size(); idx_channel++){
        double bias_volt, overvoltage, err_bias_volt, err_overvoltage;
        give_me_Bias_OV_and_errors(module, bias_dac, v_brs[idx_channel], err_v_brs[idx_channel],
                                   bias_volt, overvoltage, err_bias_volt, err_overvoltage);
        
        for(size_t j=0; j<ch_rms[0].second.size(); j++){
          if(ch_rms[1].second[j] == vgains[idx_file] && int(ch_rms[2].second[j]) == module_channels[idx_channel]){
            a.bsl = allowed_bsl_rms*ch_rms[3].second[j];
            a.sat_up = a.bsl*10;
          }
        }
        
        a.wf_file = files[idx_file]+"/channel_"+module_channels[idx_channel]+".dat";
        std::cout << a.wf_file << std::endl;
        ifstream this_file(a.wf_file);
        if (!this_file.is_open()){
          std::cout << "File not found: " << a.wf_file << std::endl;
          this_file.close();
          continue;
        }
        std::cout << "\n\n\nReading file: " << a.wf_file << std::endl;  
        cout << a.wf_file << endl;
        a.LED_Analysis();
        a.LoadFitParameters(a.fgaus);
        a.SPE();

        a.h_charge->SetTitle(Form("Ch_%i_Bias_%.2f_VGain_%i", module_channels[idx_channel], bias_volt, vgains[idx_file]));
        a.h_charge->SetName(Form("Ch_%i_Bias_%.2f_VGain_%i", module_channels[idx_channel], bias_volt, vgains[idx_file]));
        h_charge_vec.push_back(a.h_charge);
        feature_value.push_back({"Run", double(runs[idx_file])});
        feature_value.push_back({"Channel", double(module_channels[idx_channel])});
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
        feature_value.push_back({"RMS", ch_rms[3].second[idx_file]});
        feature_value.push_back({"CX", a.cx});
        feature_value.push_back({"Err CX", a.err_cx});
        feature_value.push_back({"Avg #ph cx", a.avg_n_ph_cx});
        feature_value.push_back({"Err #ph cx", a.err_avg_n_ph_cx});
        feature_value.push_back({"Avg #ph", a.avg_n_photons});
        feature_value.push_back({"Avg #pe", a.avg_n_photoelectrons});
        
        if(a.print==true){
          std::cout << "\n\nPRINTING\n\n" << std::endl;
          print_vec_pair_csv(out_csv_file, feature_value);
        }
        
        // Reset the vector
        feature_value = {};
      }
    }
    
    std::cout << "\n\nOUT OF THE LOOP\n\n" << std::endl;
    if(a.print==true) hf.cd("chargehistos"); for(auto h : h_charge_vec) h->Write();
    hf.Close();
  }
}
