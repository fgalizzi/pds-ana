#include "../classe.hpp"

void cla::TooAnnoying(){
  TString runs_folder = "/eos/home-f/fegalizz/ProtoDUNE_HD/VGain/files/run_";
  TString ana_folder  = "/eos/home-f/fegalizz/ProtoDUNE_HD/VGain/analysis";
  //std::vector<int> runs = {29947,30053,30101,30149,30263,30314,30377,30458,30539,30620,30703,
    //                        30783,30945,31187,31166};
  std::vector<int> runs = {30620};
  std::vector<TString> files = {};
  vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 
  for(auto& run : runs){
    files.push_back(runs_folder+run+"/104");
  }
  vector<int> channel_this_mask = {01, 02, 03, 04, 06, 10, 11, 12, 13, 14, 15, 17, 22,
                  24, 25, 26, 30, 31, 32, 33, 34, 35, 36, 37, 41, 42, 43, 44, 45, 46, 47};

  string path(ana_folder.Data());

  TFile hf(Form("%s/VGain_results_Ep_104_run_%d.root",path.c_str(),runs[0]), "recreate");
  hf.mkdir("chargehistos");
  hf.cd("chargehistos");
  //vector<TH1D*> h_charge_vec;
  
  print = 1;
  display=0;
  
  for(size_t i=0; i<files.size(); i++){
    for(auto&ch : channel_this_mask){
      n_wf = 11000;
      wf_file = files[i]+"/channel_"+ch+".dat";
      ifstream this_file(wf_file);
      if (!this_file.is_open()){
        this_file.close();
        continue;
      }
      std::cout << "\t\tReading file: " << wf_file << std::endl;  
      // cout << wf_file << endl;
      bsl = spe_ampl*1.5;
      LED_Analysis();
      LoadFitParameters(fgaus);
      SPE();
      //h_charge_vec.push_back(h_charge);
      feature_value.push_back({"Channel", double(ch)});
      feature_value.push_back({"Spe ampl", spe_ampl});
      feature_value.push_back({"DR", pow(2,14)/(spe_ampl-spe_under)});
      feature_value.push_back({"SNR", spe_charge/sigma_zero});
      if(print==true)h_charge->Write();
      if(print==true) print_vec_pair_csv(Form("%s/VGain_results_Ep_104_run_%d.csv",path.c_str(),runs[0]), feature_value);
      feature_value = {};
      // cout << Form("%s/VGain_results_Ep_104_run_%d.csv",path.c_str(),runs[0]) << endl;
    }
  }
  
  //hf.cd("chargehistos"); for(auto h : h_charge_vec) h->Write();
  hf.Close();
}
