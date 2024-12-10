#include "../classe.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

void cla::read(){
  if(wfs.size()!=n_wf || oldwf_file!=wf_file || oldprepulse_ticks != prepulse_ticks ||
      oldchannel!=channel || wfs[0].size()!=memorydepth){
    //Store old values to decide whether to re-read the file 
    oldwf_file = wf_file;
    oldprepulse_ticks = prepulse_ticks;
    oldchannel = channel;

    ite = 0; // For the macros whose behaviour must be different at the first iteration

    wfs.erase(wfs.begin(),wfs.end());
    
    if(data_format == "caen")    CAEN_WF_Binary(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "daphne")  CompleteWF_Binary(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "esteban") CompleteWF_Binary_Swap(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "csv")     CSV_WF_Binary(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "csvd")    CSV_double_WF_Binary(wf_file, wfs, n_wf, memorydepth);
    
    // ProtoDUNE-HD: update n_wf because the reading function stops automatically
    if(data_format == "pdhd"){
      PDHD_ch_wfs(wf_file, wfs, channel, n_wf); 
      n_wf = wfs.size();
    }
   
    //Subtract the baseline and invert the wfs according to "invert"
    if(sub_bsl == true){
      if(sub_bsl_mode == 1) SubBaseline(wfs, prepulse_ticks, invert);
      if(sub_bsl_mode == 2) SubBaseline2(wfs, rms, invert);
    }
  }

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
