#include "../../../Class/_c/class_include.hpp"
using namespace std;

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void NoiseVGain_ana(){
  // --------------------------------------------------------------
  // --- HARD CODE ------------------------------------------------
  int module = 3;
  int day = 4; //3 or 4
  string base_path = Form("/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/December2024run/CAEN/M%i/Noise/2024120%i/", module, day);
  string output_ana_folder  = "/eos/home-g/gpiemont/ColdBox_VD/December24/CAEN/Noise/";
  // Class settings
  auto a = cla();
  a.data_format = "caen";
  a.n_wf = -1; // Full statistics
  a.memorydepth = 5000;
  a.prepulse_ticks = 2619;
  a.tick_len = 0.016;
  a.res = 14;
  a.print = 1;

  // --- CODE ----------------------------------------------------
  string fft_outfile, ana_outfile;
  fft_outfile = output_ana_folder+Form("FFT_VGainScans_Day%i_Membrane%i.root", day, module);
  ana_outfile = output_ana_folder+Form("VGain_RMS_Day%i_Membrane%i.csv", day, module);
  
  vector<int> channels = {0,1};

 /* vector<int> runs;
  for(int i=9; i<300; i++){
    if((i+(LED+1))%2==0) runs.push_back(33300+i);
    if(i==66) i=10000;
  }
*/
  int vgain = 1800;
  // Tuple input file, run, vgain, channel
  vector<tuple<string,int>> my_tuple;
  vector<string> ifiles;
  string led;
  if(day == 4){
  	if (module == 3) led = "_LED7p25_";
  	else led = "_LED6p75_";
  }
  else led = "_";
  //for(auto& run : runs){
    for(auto& ch : channels){    
  string ifile = base_path+Form("M%i", module)+led+Form("Noise_Ch%i.dat", ch);
      my_tuple.push_back(make_tuple(ifile,ch));
    }
    //vgain += 100;
  //}

  TFile fft_ofile(fft_outfile.c_str(), "recreate");
  for(auto& ch : channels){
    string dir_name = "Ch_"+to_string(ch);
    fft_ofile.mkdir(dir_name.c_str());
  }
  std::cout << "dd" << std::endl;
  // --- LOOP OVER FILES -----------------------------------------
  for(auto& tuple : my_tuple){
    a.wf_file = get<0>(tuple);
    a.read();
    TGraph* gNoise_spectral_density = build_avg_spectral_density(a.memorydepth,
      a.tick_len*a.memorydepth, a.tick_len, a.wfs, a.res);

    string dir_name = "Ch_"+to_string(get<1>(tuple));
    fft_ofile.cd(dir_name.c_str());
    gNoise_spectral_density->Write();
    if (a.print == true){
        vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 

        feature_value.push_back({"Ch", get<1>(tuple)});
        feature_value.push_back({"Rms", standard_deviation_vec_vec(a.wfs)});

        print_vec_pair_csv(ana_outfile, feature_value);
    }
  }
}
