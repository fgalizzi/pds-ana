#include "Utils.hpp"

using namespace std;
using json = nlohmann::json;

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void Calibration_ana(cla& a, string jsonfile_module_config){
  // --- ANA CONFIG -----------------------------------------------
  AnaConfig ana_config = load_ana_config("config/ana_config.json");
  string runs_folder = ana_config.runs_folder;
  string input_ana_folder = ana_config.input_ana_folder;
  string rms_result_file = input_ana_folder+ana_config.rms_result_file;
  string channel_map_file = ana_config.channel_map_file;
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

  vector<int> channels = module_config.channels;
  vector<int> runs = module_config.runs;
 
  // --- ANALYSIS -------------------------------------------------
  vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 

  string out_files_name = Form("TimeResolutionCalib_Run_%i", runs[0]);
  string out_root_file = output_ana_folder+out_files_name+".root";
  string out_csv_file  = output_ana_folder+out_files_name+".csv";

  std::vector<TString> files = {};

  TFile hf(TString(out_root_file), "recreate");
  hf.mkdir("chargehistos");
  hf.cd("chargehistos");
  vector<TH1D*> h_charge_vec;

  vector<pair<string, vector<double>>> ch_rms = read_vec_pair_CSV(rms_result_file.c_str()); // Store the results of the analysis to be printed
  map<int,double> offlinech_to_rms = {};
  for(size_t j=0; j<ch_rms[0].second.size(); j++){
    offlinech_to_rms[int(ch_rms[0].second[j])] = ch_rms[3].second[j];
  }
  vector<pair<string, vector<double>>> channel_map = read_vec_pair_CSV(channel_map_file.c_str()); // Store the results of the analysis to be printed
  map<int, int> daphnech_to_offlinech = {};
  for(size_t j=0; j<channel_map[0].second.size(); j++){
    daphnech_to_offlinech[int(channel_map[0].second[j]*100+channel_map[1].second[j])] = int(channel_map[2].second[j]);
  }

  // --- LOOP OVER THE FILES == RUNS ----------------------------
  for(auto run : runs){
    // --- LOOP OVER THE CHANNELS --------------------------------
    for(auto channel : channels){
      // Load the file and check if it exists
      a.n_wf = module_config.n_wf;
      a.channel = channel;
      a.wf_file = runs_folder+"processed_merged_run_"+run+"_structured.hdf5";
      ifstream this_file(a.wf_file);
      if (!this_file.is_open()){
        std::cout << "File not found: " << a.wf_file << std::endl;
        this_file.close();
        continue;
      }
      std::cout << "Run " << run << " ch " << channel
                << "\n" << a.wf_file << std::endl;

      a.bsl = allowed_bsl_rms*offlinech_to_rms[daphnech_to_offlinech[channel]];
      a.sat_up = a.bsl*10;

      // The actual analysis
      a.LED_Analysis();
      a.LoadFitParameters(a.fgaus);
      a.SPE();

      // --- OUTPUT ------------------------------------------------
      a.h_charge->SetTitle(Form("DaphneCh_%i", channel));
      a.h_charge->SetName(Form("DaphneCh_%i", channel));
      h_charge_vec.push_back(a.h_charge);
      feature_value.push_back({"Run", double(run)});
      feature_value.push_back({"DaphneCh", double(channel)});
      feature_value.push_back({"Baseline", a.bsl});
      feature_value.push_back({"Prepulse ticks", double(a.prepulse_ticks)});
      feature_value.push_back({"Saturation up", a.sat_up});
      feature_value.push_back({"Int low", double(a.int_low)});
      feature_value.push_back({"Int up", double(a.int_up)});
      feature_value.push_back({"Gain", a.spe_charge});
      feature_value.push_back({"Err Gain", a.err_spe_charge});
      feature_value.push_back({"SpeAmpl", a.spe_ampl});
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
