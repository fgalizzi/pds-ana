#include "../classe.hpp"

void cla::TooAnnoying(){
  TString runs_folder = "/eos/home-f/fegalizz/ProtoDUNE_HD/VGain/files/";
  TString ana_folder  = "/eos/home-f/fegalizz/ProtoDUNE_HD/VGain/analysis";
  std::vector<int> runs = {29947,30004};
  std::vector<TString> files = {};
  vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 
  for(auto& run : runs){
    files.push_back(runs_folder+run+"/104");
  }
  vector<int> channel_this_mask = {01, 02, 03, 04, 06, 10, 11, 12, 13, 14, 15, 17, 22,
                  24, 25, 26, 30, 31, 32, 33, 34, 35, 36, 37, 41, 42, 43, 44, 45, 46, 47};


  print = 1;
  display=1;
  
  for(size_t i=0; i<files.size(); i++){
    for(auto&ch : channel_this_mask){
      n_wf = 2000;
      wf_file = files[i]+"/channel_"+ch+".dat";
      // cout << wf_file << endl;
      bsl = spe_ampl*1.2;
      LED_Analysis();
      LoadFitParameters(fgaus);
      SPE();
      feature_value.push_back({"Channel", double(ch)});
      feature_value.push_back({"Spe ampl", spe_ampl});
      feature_value.push_back({"DR", pow(2,14)/(spe_ampl+spe_under)});
      feature_value.push_back({"SNR", spe_charge/sigma_zero});
      string path(ana_folder.Data());
      // cout << Form("%s/VGain_results_Ep_104.csv",path.c_str()) << endl;
      if(print==true) print_vec_pair_csv(Form("%s/VGain_results_Ep_104.csv",path.c_str()), feature_value);
    }
  }
}
