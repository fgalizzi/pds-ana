#include "../classe.hpp"

void cla::TooAnnoying(){

//Scan tick
int a = int_up;
for(int i = a; i < a+30; i++){
  int_up = i;
  LED_Analysis();
}
  





/*  TString runs_folder = "/eos/home-f/fegalizz/ProtoDUNE_HD/VGain/files/run_";
  TString ana_folder  = "/eos/home-f/fegalizz/ProtoDUNE_HD/VGain/analysis";
  //std::vector<int> runs = {29947,30053,30101,30149,30263,30314,30377,30458,30539,30620,30703,
    //                        30783,30945,31187,31166};

  std::vector<int> runs = {30620};

  //std::sort(runs.begin(), runs.end());
  //std::vector<int> runs = {30620};
  std::vector<TString> files = {};
  vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 
  for(auto& run : runs){
    files.push_back(runs_folder+run+"/104");
  }
//  vector<int> channel_this_mask = {01, 02, 03, 04, 06, 10, 11, 12, 13, 14, 15, 17, 22,
//                  24, 25, 26, 30, 31, 32, 33, 34, 35, 36, 37, 41, 42, 43, 44, 45, 46, 47};

  vector<int> channel_this_mask = {30};
  string path(ana_folder.Data());

  //vector<TH1D*> h_charge_vec;
  //bsl = 500;
  //sat_up = 4000;
  //sat_low = -4000;
  print = 1;
  display=0;
  
  for(size_t i=0; i<files.size(); i++){
    for(auto&ch : channel_this_mask){
      TFile hf(Form("%s/VGain_results_Ep_104_ch_30.root",path.c_str()), "update");
      hf.mkdir("chargehistos");

      //wf_file = files[i]+"/channel_"+ch+".dat";
      ifstream this_file(wf_file);
      if (!this_file.is_open()){
        this_file.close();
        continue;
      }
      std::cout << "\t\tReading file: " << wf_file << std::endl;  
      // cout << wf_file << endl;
      //bsl = spe_ampl*2.5;
      //sat_up = spe_ampl*10;
      //sat_low= -spe_ampl*10;
      LED_Analysis();
      if (class_skip == 1){
        class_skip=0;
        continue;
      }
      LoadFitParameters(fgaus);
      if (class_skip == 1){
        class_skip=0;
        continue;
      }
      SPE();
      if (class_skip == 1){
        class_skip=0;
        continue;
      }
      //h_charge_vec.push_back(h_charge);
      feature_value.push_back({"Run", double(runs[i])});
      feature_value.push_back({"Channel", double(ch)});
      feature_value.push_back({"Spe ampl", spe_ampl});
      feature_value.push_back({"DR", pow(2,14)/(spe_ampl-spe_under)});
      feature_value.push_back({"SNR", spe_charge/sigma_zero});
      feature_value.push_back({"Good", 1});

      h_charge->SetName(Form("Run_%d_ch_%d",runs[i],ch));
      if(print==true){
        h_charge->Write();
        print_vec_pair_csv(Form("%s/VGain_results_Ep_104_ch_30.csv",path.c_str()), feature_value);
      }
      feature_value = {};
      // cout << Form("%s/VGain_results_Ep_104_run_%d.csv",path.c_str(),runs[0]) << endl;
      hf.Close();
    }
  }
  
  //hf.cd("chargehistos"); for(auto h : h_charge_vec) h->Write();*/
}
