#include "../classe.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

void cla::read(){
  if(int(wfs.size())!=n_wf || oldwf_file!=wf_file || oldprepulse_ticks != prepulse_ticks ||
      oldchannel!=channel || int(wfs[0].size())!=memorydepth){
    //Store old values to decide whether to re-read the file 
    oldwf_file = wf_file;
    oldprepulse_ticks = prepulse_ticks;
    oldchannel = channel;

    ite = 0; // For the macros whose behaviour must be different at the first iteration

    wfs.erase(wfs.begin(),wfs.end());

    bool consistency_checks = true;

    if (prepulse_ticks > memorydepth){
      std::cout << "\n\nprepulse_ticks > memorydepth !!\n\n" << std::endl;
      std::cout << "Setting prepulse_ticks = memorydepth" << std::endl;
      prepulse_ticks = memorydepth;
      consistency_checks = false;
    }
    if (int_low > int_up){
      std::cout << "\n\nint_low > int_up !!\n\n" << std::endl;
      std::cout << "Setting int_low = int_up" << std::endl;
      int_low = int_up-1;
      consistency_checks = false;
    }
    if (prepulse_ticks+int_up > memorydepth){
      std::cout << "\n\nprepulse_ticks+int_up > memorydepth !!\n\n" << std::endl;
      std::cout << "Setting int_up = 0" << std::endl;
      int_up = 0;
      consistency_checks = false;
    }

    if (consistency_checks == false){
      std::cout << "\n\n\n---Consistecy checks failed---\n\n\n" << std::endl;
    }
    
    if(data_format == "caen")    CAEN_WF_Binary(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "daphne")  CompleteWF_Binary(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "esteban") CompleteWF_Binary_Swap(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "csv")     CSV_WF_Binary(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "csvd")    CSV_double_WF_Binary(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "pdhd")    PDHD_ch_wfs(wf_file, wfs, channel, n_wf); 
    if(data_format == "hdf5")    StructuredWaveformSetReader(wf_file, wfs, channel, n_wf);
    if(data_format == "eth")     StructuredEthWaveformSetReader(wf_file, wfs, channel, n_wf);
  
   
    //Subtract the baseline and invert the wfs according to "invert"
    if(sub_bsl == true){
      if(sub_bsl_mode == 1) SubBaseline(wfs, prepulse_ticks, invert);
      if(sub_bsl_mode == 2) SubBaseline2(wfs, rms, invert);
    }
  }

}

void cla::waveforms_from_multifile(const vector<string>& files){
  vector<vector<vector<double>>> wfs_temp;
  for (auto& file : files){
    wfs.clear();
    wf_file = file;
    read();
    wfs_temp.push_back(wfs);
  }
  vectorVector_to_vector(wfs, wfs_temp);
  n_wf = wfs.size();
  wfs_temp.clear();
}

//Read the ProtoDUNE-HD channel-map
vector<size_t> cla::read_pdhd_ch_map(int mask){
  string file;
  if(mask==0) file = class_path+"/Class/ProtoduneHD/channelmap.txt";
  if(mask!=0) file = class_path+"/Class/ProtoduneHD/channel_led.txt";
  
  ifstream file_map(file);
  string line;
  stringstream ssmap;
  short dpch, ch, mk;
  vector<size_t> channels;

  if (file_map.is_open()){
    while (getline(file_map, line)) {
      ssmap.clear();
      ssmap.str(line);

      if(mask==0){
        while (ssmap >> dpch >> ch) {
          channels.push_back(size_t(dpch));
        }
      }

      if(mask!=0){
        while (ssmap >> dpch >> mk) {
          if(mask==mk) channels.push_back(size_t(dpch));
        }
      }

    }
  }
  else std::cout << "The PDHD channel map is not here: " << file << std::endl;

  return channels;
}

// Read file.txt with single column (channel name)
//**********************************************************
vector<string> cla::read_chs(string ch_file_name){
//**********************************************************
  string file_folder = class_path+"/Class/ProtoduneHD/";
  ifstream ch_file(file_folder+ch_file_name);
  string line;
  vector<string> chs;

  if (ch_file.is_open()){
    while (getline(ch_file, line)) {
      chs.push_back(line);
    }
  }

  return chs;
}
