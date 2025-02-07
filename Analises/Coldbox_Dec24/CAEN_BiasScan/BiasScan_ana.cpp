#include "../../../Class/_c/class_include.hpp"
using namespace std;

//------- Macro ---------------------------------------------------
void Loop_VBias_Scan_Caen(){
  // -------------------------------------------------------------
  // --- HARD CODE -----------------------------------------------

  std::vector<double> biases = {};
  std::vector<TString> bias_volt = {};
  int module = 3;
  int day = 3; //Day of data taking 3 or 4

  // INPUT
  TString runs_folder = Form("/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/CAEN/M%i/VBias_Scans/2024120%i/", module, day);
  TString output_ana_folder = "/eos/home-g/gpiemont/ColdBox_VD/December24/CAEN/";

  TString name_file;
  if(module == 1){
	 if(day == 3) name_file = Form("M%i_LED6p5_", module);
	 else name_file = Form("M%i_LED6p30_", module);
  }
  if(module == 2) {
	  if(day == 3) name_file = Form("M%i_LED6p35_", module);
	  else name_file = Form("M%i_LED6p2_", module);
  if(module == 3) name_file = Form("M%i_LED7p25_", module);
  if(module == 4) name_file = Form("M%i_LED6p75_", module);
 
  
  if (module == 1 || module == 2){
          biases = {45, 45.5, 46, 46.5, 47};
          bias_volt = {"45p0V", "45p5V", "46p0V", "46p5V", "47p0V"};
  }
  if (module == 3 || module == 4){
  	if(day == 3){
          biases = {31, 31.5, 32, 32.5, 33, 33.5};
          bias_volt = {"31p0V", "31p5V", "32p0V", "32p5V", "33p0V", "33p5V"};
  	}
  	else{
          biases = {29, 29.5, 30, 30.5, 31, 31.5, 32, 32.5, 33, 33.5};
          bias_volt = {"29p0", "29p5", "30p0", "30p5","31p0V", "31p5V", "32p0V", "32p5V", "33p0V", "33p5V"};
  	}
  }

  vector<int> channel_this_mask = {0, 1}; //CAEN channels
  // Initial sat_up for the scan
  double scan_sat_up = 1600;

  // CLASS SETTINGS
  display = 0;
  print   = 0;
  plot    = 0;

  // OUTPUT
  bool print_results = true;
  TString out_file = output_ana_folder+Form("VBias_Scan_%i_Module_%i", day, module);
  TString out_root_file = out_file+".root";
  string out_csv_file(out_file+".csv");
  // --- END HARD CODE -------------------------------------------


  // --- CODE ----------------------------------------------------

  vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed

  std::vector<TString> files = {};
  for(auto& bias_volt : bias_volt){
    files.push_back(name_file+bias_volt);
  }

  if(biases.size() != files.size()){
    return;
  }

  TFile hf(out_root_file, "recreate");
  hf.mkdir("chargehistos");
  hf.cd("chargehistos");
  vector<TH1D*> h_charge_vec;


  std::cout << "files " << files.size() << std::endl;
  for(size_t i=0; i<files.size(); i++){
    for(auto&ch : channel_this_mask){
      wf_file = runs_folder+files[i]+"_Ch"+ch+".dat";
      std::cout << wf_file << std::endl;
      ifstream this_file(wf_file);
      if (!this_file.is_open()){
        std::cout << "File not found: " << wf_file << std::endl;
        this_file.close();
        continue;
      }
      std::cout << "\n\n\nReading file: " << wf_file << std::endl;
      cout << wf_file << endl;
      sat_up = scan_sat_up;
      LED_Analysis();
      LoadFitParameters(fgaus);
      SPE();
      sat_up = spe_ampl*20;
      LED_Analysis();
      LoadFitParameters(fgaus);
      SPE();

      h_charge->SetTitle(Form("VBias_%i_ch_%i",int(biases[i]),ch));
      h_charge->SetName(Form("VBias_%i_ch_%i",int(biases[i]),ch));
      h_charge_vec.push_back(h_charge);

      /*int arap_ch;
      if (ch == 0 || ch == 1 || ch == 27 || ch == 26) arap_ch = 1;
      else arap_ch = 2;
      */

      //feature_value.push_back({"Membrane modules Channel", int(arap_ch)});
      feature_value.push_back({"CAEN Channel", double(ch)});
      //feature_value.push_back({"Bias [dac]", biases[i]});
      feature_value.push_back({"Bias [V]", biases[i]});
      // feature_value.push_back({"VGain", double(vgains[i])});
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
 
