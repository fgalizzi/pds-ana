#include "../../../Class/_c/class_include.hpp"
#include "ana_parameters.hpp"

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void BiasGainScan_ana(){
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
  auto a = cla();
  a.n_wf = 5000;
  a.display=0;
  a.print = 0;

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
          a.bsl = 4.*ch_rms[3].second[j];
          a.sat_up = a.bsl*10;
        }
      }
      a.wf_file = files[i]+"/channel_"+ch+".dat";
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

      a.h_charge->SetTitle(Form("VGain_%i_ch_%i",vgains[i],ch));
      a.h_charge->SetName(Form("VGain_%i_ch_%i",vgains[i],ch));
      h_charge_vec.push_back(a.h_charge);
      feature_value.push_back({"Channel", double(ch)});
      feature_value.push_back({"Bias [dac]", bias_dac});
      feature_value.push_back({"Bias [V]", bias_volts});
      feature_value.push_back({"VGain", double(vgains[i])});
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
      feature_value.push_back({"RMS", ch_rms[3].second[i]});
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
