#include "../classe.hpp"

// ****************************************************************
// Change this marco when you can easily loop over runs and you
// don't want to repeat the same commands every time!
// If the loop could be useful for the future, store it in Loops.cpp
// ****************************************************************

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////


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
void cla::TooAnnoying(){
  // -------------------------------------------------------------
  // --- HARD CODE -----------------------------------------------
  // INPUT
  string runs_folder        = "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/Daphne_DAQ/binaries/";
  string input_ana_folder   = "/eos/home-f/fegalizz/ColdBox_VD/December24/Daphne_DAQ/Noise_RMS_FFTs/";
  
  // Module, channels, bias_dac, runs and corresponding VGains
  // Take from ~/pds-ana/Analises/Coldbox_Dec24/Bias_and_VGain_scan/Run_Bias_VGain_correspondence.hpp
int module = 2;
std::vector<int> module_channels = {21, 26};
std::vector<double> v_brs     = {42.71, 42.59};
std::vector<double> err_v_brs = {0.05, 0.06};
std::vector<int> vgains = {0,    100,  200,  300,  400,  500,  600,  700,  800,  900,
                           1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,
                           2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900,
                           3000};
std::vector<double> bias_dacs = {1200, 1187, 1174, 1161, 1148};
std::vector<std::vector<int>> run_batches = {{34230, 34231, 34232, 34233, 34234, 34235, 34236, 34237, 34238, 34239,
                                      33380, 33381, 33382, 33383, 33384, 33385, 33386, 33387, 33388, 33389,
                                      33390, 33392, 33393, 33394, 33395, 33396, 33397, 33398, 33399, 33400,
                                      34043},
                                     {34240, 34241, 34242, 34243, 34244, 34245, 34246, 34247, 34248, 34249,
                                      33401, 33402, 33403, 33404, 33405, 33406, 33407, 33408, 33409, 33410,
                                      33411, 33413, 33414, 33415, 33416, 33417, 33418, 33419, 33420, 33421, 
                                      34044},
                                     {34250, 34251, 34252, 34253, 34254, 34255, 34256, 34257, 34258, 34259,
                                      33660, 33661, 33662, 33663, 33464, 33665, 33666, 33667, 33668, 33669, 
                                      33670, 33671, 33672, 33673, 33674, 33675, 33676, 33677, 33678, 33679,
                                      34045},
                                     {34260, 34261, 34262, 34263, 34264, 34265, 34266, 34267, 34268, 34269,
                                      33680, 33681, 33682, 33683, 33684, 33685, 33686, 33687, 33688, 33689, 
                                      33690, 33691, 33692, 33693, 33694, 33695, 33696, 33697, 33698, 33699,
                                      34046},
                                     {34270, 34271, 34272, 34273, 34274, 34275, 34276, 34277, 34278, 34279,
                                      33700, 33701, 33702, 33703, 33704, 33705, 33706, 33707, 33708, 33709, 
                                      33710, 33711, 33712, 33713, 33714, 33715, 33716, 33717, 33718, 33719,
                                      34047}};
  
  // File with the RMS of the channels. bsl = allowed_bsl_rms*channel_rms
  string rms_result_file = input_ana_folder+"VGain_RMS_LED0_Membrane.csv";
  double allowed_bsl_rms = 3.;

  // CLASS SETTINGS
  display = 0;
  print   = 0;
  plot    = 0;

  // OUTPUT
  bool print_results = true;
  string output_ana_folder  = "/eos/home-f/fegalizz/ColdBox_VD/December24/Daphne_DAQ/VGain_Scans/";

  // --- END HARD CODE -------------------------------------------
  // -------------------------------------------------------------
 

  // --- CODE ----------------------------------------------------
  
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
            bsl = allowed_bsl_rms*ch_rms[3].second[j];
            sat_up = bsl*10;
          }
        }
        
        wf_file = files[idx_file]+"/channel_"+module_channels[idx_channel]+".dat";
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

        h_charge->SetTitle(Form("Ch_%i_Bias_%.2f_VGain_%i", module_channels[idx_channel], bias_volt, vgains[idx_file]));
        h_charge->SetName(Form("Ch_%i_Bias_%.2f_VGain_%i", module_channels[idx_channel], bias_volt, vgains[idx_file]));
        h_charge_vec.push_back(h_charge);
        feature_value.push_back({"Run", double(runs[idx_file])});
        feature_value.push_back({"Channel", double(module_channels[idx_channel])});
        feature_value.push_back({"Bias [dac]", bias_dac});
        feature_value.push_back({"Bias [V]", bias_volt});
        feature_value.push_back({"Err Bias [V]", err_bias_volt});
        feature_value.push_back({"OV [V]", overvoltage});
        feature_value.push_back({"Err OV [V]", err_overvoltage});
        feature_value.push_back({"VGain", double(vgains[idx_file])});
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
        feature_value.push_back({"RMS", ch_rms[3].second[idx_file]});
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
}
