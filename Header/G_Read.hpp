//
//  G_Read.hpp
//
//  Created by Federico Galizzi on 28/09/23
//

#ifndef G_Read_hpp
#define G_Read_hpp

#include <stdio.h>

// Read Binary file with #WF=WFs len-tick long
//**********************************************************
void CompleteWF_Binary(std::string fileName, vector<vector<double>>& y, int WFs, int len){
//**********************************************************
  y.resize(WFs, vector<double>(len));
  double t;
  std::ifstream file;

  file.open( fileName, std::ios::binary );
  if (file.is_open()) std::cout << "Running on " << fileName << std::endl;
  else {
    std::cout << "Error opening file" << std::endl;
    return;
  }
  
  for (int n_wf=0; n_wf < WFs; n_wf++) {
    for (int i=0; i < len; i++) {
      file.read( reinterpret_cast<char*>(&t), sizeof(t) );
      y[n_wf][i] = t;
    }
  }
  std::cout << "The file has been correctly read \t \t" << std::endl;
}

// Open a decoded ProtoDUNE .root file 
//**********************************************************
void PDHD_ch_wfs(std::string fileName, vector<vector<double>>& y, int this_ch, int WFs){
//**********************************************************
  int len = 1024;
  int wf_counter = 0;
  y.resize(WFs, vector<double>(len));

  std::ifstream file_map("/Users/federico/pds-ana/Class/ProtoduneHD/channelmap.txt"); // reading channel map
  Short_t sl, lk, dpch, ch;
  std::stringstream ssmap;

  std::map<Short_t, Short_t> detmap;
  std::map<Short_t, Short_t> invdetmap;

  std::string line;

  if (file_map.is_open()){
      while (getline(file_map, line)){
          ssmap.clear();
          ssmap.str(line);

          while (ssmap >> dpch >> ch){
              detmap[dpch] = ch;
              invdetmap[ch] = dpch;
          }
      }
      file_map.close();
  }
  else{
    std::cerr << "Unable to open file!!!" << std::endl;
    file_map.close();
  }

  wffunctions bs;

  map<TString, TString> filename = {{"file1", fileName}};

  for (auto f : filename){

     TChain *t[] = {NULL};

     t[0] = new TChain();
     t[0]->Add(Form("%s?#raw_waveforms", f.second.Data()));
     // t[0]->SetImplicitMT(true);
     Long64_t nentries = t[0]->GetEntries();
     hdf5torootclass event(t[0]);

     cout << "\nFile open -> " << f.second << "\tentries: " << nentries << endl;

       for (size_t ievt=0; ievt<nentries && wf_counter<WFs; ievt++){ // loop over entries in root file

          event.GetEntry(ievt);

          if (event.channel == this_ch){
            if (event.is_fullstream) continue;

            bs.setADCvector(event.adcs); // setting the adc vector to use function

            for (int i = 0; i < len; i++)  y[wf_counter][i] = event.adcs->at(i);
            wf_counter++;
          }
       }
  }
  std::cout << "The file has been correctly read \t \t" << std::endl;
  std::cout << wf_counter << std::endl;
}

// Read CAEN-Wavedump Binary file with #WF=WFs len-tick long
//**********************************************************
void CAEN_WF_Binary(std::string fileName, vector<vector<double>>& y, int WFs, int len){
//**********************************************************
  y.resize(WFs, std::vector<double>(len));
  int skip;
  uint16_t t;
  std::ifstream file;

  file.open( fileName, std::ios::binary );
  if (file.is_open()) std::cout << "Running on " << fileName << std::endl;
  else {
    std::cout << "Error opening file" << std::endl;
    return;
  }
  
  for (int n_wf=0; n_wf < WFs; n_wf++) {
    for (int i=0; i < 6; i++) file.read( reinterpret_cast<char*>(&skip), sizeof(skip) );
    for (int i=0; i < len; i++) {
      file.read( reinterpret_cast<char*>(&t), sizeof(t) );
      y[n_wf][i] = double(t);
    }
  }
  std::cout << "The file has been correctly read \t \t" << std::endl;
}

//*********************************************
void CompleteWF_Binary_Swap(std::string fileName, vector<vector<double>>& y, int WFs, int len){
//*********************************************
  y.resize(WFs, vector<double>(len)); 
  int skip;
  uint16_t t;
  uint8_t t_8_msb, t_8_lsb;
  std::ifstream file;

file.open( fileName, std::ios::binary );
  if (!file.is_open()) {
    std::cout << "Error opening file " << std::endl;
    return;
  }
  std::cout << "Running on " << fileName << std::endl;
  
  for (int n_wf=0; n_wf<WFs; n_wf++) {
    for (int i=0; i<len; i++) {
      file.read( reinterpret_cast<char*>( &t ), sizeof( t ) );
      //t = OSSwapBigToHostInt16(t);
      y[n_wf][i]= double(t);
    }
  }
  
  std::cout << "The file has been correctly read \t \t" << std::endl;
  
}

// Read CSV file with #WF=WFs len-tick long
//**********************************************************
void CSV_WF_Binary(std::string fileName, vector<vector<double>>& y, int WFs, int len){
//**********************************************************
  y.resize(WFs, vector<double>(len));
  int t;
  std::ifstream file;

  file.open( fileName, std::ios::in);
  if (file.is_open()) std::cout << "Running on " << fileName << std::endl;
  else {
    std::cout << "Error opening file" << std::endl;
    return;
  }
  
  for (int n_wf=0; n_wf < WFs; n_wf++) {
    for (int i=0; i < len; i++) {
     file >> t; 
      y[n_wf][i] = double(t);
    }
  }
  std::cout << "The file has been correctly read \t \t" << std::endl;
}

// Read a .txt and put in a vector
//**********************************************************
void CompleteWF(std::string fileName, vector<double>& y){
//**********************************************************
  double t;
  
  std::fstream file (fileName.c_str());
  if (file.is_open()) std::cout << "Running on " << fileName << std::endl;
  else {
    std::cout << "Error opening file" << std::endl;
    return;
  }
  while (true) {
    file >> t;
    if (file.eof() == true ){
      cout << "The file has been correctly read \t \t" <<std::endl;
      return;}
    y.push_back(t);
  }
  return;
}

#endif /* G_Read_hpp */
