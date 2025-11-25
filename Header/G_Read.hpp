///
//  G_Read.hpp
//
//  Created by Federico Galizzi on 28/09/23
//

#ifndef G_Read_hpp
#define G_Read_hpp


#include "RtypesCore.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <H5Cpp.h>
#ifndef hdf5torootclass_cxx
#include "../Class/ProtoduneHD/wffunctions2.h"
#include "../Class/ProtoduneHD/hdf5torootclass.h"
#endif // !hdf5torootclass_cxx
using namespace std;

// ---------------------------------------------------------
// --- METHODS TO READ WAVEFORMS FROM DIFFERENT SOURCES ----
// ---------------------------------------------------------

// Read Binary file with #WF=WFs len-tick long
//**********************************************************
void CompleteWF_Binary(std::string fileName, vector<vector<double>>& y, int WFs, int len){
//**********************************************************
  if (WFs < 0) WFs = 100000;
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
      file.read( reinterpret_cast<char*>(&t), sizeof(t));

      if (file.eof()) {
        std::cout << "End of file reached" << std::endl;
        WFs = n_wf+1;
        break;
      }

      y[n_wf][i] = t;
    }
  }
  y.resize(WFs);
  std::cout << "The file has been correctly read \t \t" << std::endl;
  return;
}

// Single wf version
//**********************************************************
void CompleteWF_Binary(std::string fileName, vector<double>& y, int len){
//**********************************************************
  y.resize(len);
  double t;
  std::ifstream file;

  file.open( fileName, std::ios::binary );
  if (file.is_open()) std::cout << "Running on " << fileName << std::endl;
  else {
    std::cout << "Error opening file" << std::endl;
    return;
  }
  
  for (int i=0; i < len; i++) {
    file.read( reinterpret_cast<char*>(&t), sizeof(t) );
     y[i] = t;
  }
  
  std::cout << "The file has been correctly read \t \t" << std::endl;
}
// Open a decoded ProtoDUNE .root file 
//**********************************************************
void PDHD_ch_wfs(std::string fileName, vector<vector<double>>& y, int this_ch, int& WFs){
//**********************************************************
  int len = 1024;
  int wf_counter = 0;
  std::vector<double> y2;
  y2.resize(len, 0);

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

       for (Long64_t ievt=0; ievt<nentries && wf_counter<WFs; ievt++){ // loop over entries in root file

          event.GetEntry(ievt);

          if (event.channel == this_ch){
           // if (event.is_fullstream) continue;

            bs.setADCvector(event.adcs); // setting the adc vector to use function

            for (int i = 0; i < len; i++)  y2[i] = event.adcs->at(i);
            y.push_back(y2);
            wf_counter++;
          }
       }
  }
  std::cout << "The file has been correctly read \t \t" << std::endl;
  std::cout << wf_counter << std::endl;
  WFs = wf_counter;
}

// Read CAEN-Wavedump Binary file with #WF=WFs len-tick long
//**********************************************************
void CAEN_WF_Binary(std::string fileName, vector<vector<double>>& y, int& WFs, int len){
//**********************************************************
  if (WFs < 0) WFs = 100000;
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
      
      if (file.eof()) {
        std::cout << "End of file reached" << std::endl;
        WFs = n_wf+1;
        break;
      }

      y[n_wf][i] = double(t);
    }
  }
  y.resize(WFs);
  std::cout << "The file has been correctly read \t \t" << std::endl;
  return;
}

//*********************************************
void CompleteWF_Binary_Swap(std::string fileName, vector<vector<double>>& y, int& WFs, const int& len){
//*********************************************
  if (WFs < 0) WFs = 100000;
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

      if (file.eof()) {
        std::cout << "End of file reached" << std::endl;
        WFs = n_wf+1;
        break;
      }

      //t = OSSwapBigToHostInt16(t);
      y[n_wf][i]= double(t);
    }
  }
  
  y.resize(WFs);
  std::cout << "The file has been correctly read \t \t" << std::endl;
  return; 
}

// Read CSV file with #WF=WFs len-tick long
//**********************************************************
void CSV_WF_Binary(std::string fileName, vector<vector<double>>& y, int& WFs, int len){
//**********************************************************
  if (WFs < 0) WFs = 100000;
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
      if (file.eof()) {
        std::cout << "End of file reached" << std::endl;
        WFs = n_wf+1;
        break;
      }
      y[n_wf][i] = double(t);
    }
  }
  y.resize(WFs);
  std::cout << "The file has been correctly read \t \t" << std::endl;
  return;
}

// Read CSV file with #WF=WFs len-tick long
//**********************************************************
void CSV_double_WF_Binary(std::string filename, vector<vector<double>>& y, int& n_wf, int len){
//**********************************************************
  std::ifstream file(filename);

  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  std::string line;
  int row_count = 0;

  while (std::getline(file, line) && row_count < n_wf) {
    std::vector<double> row;
    std::stringstream ss(line);
    std::string value;
    int col_count = 0;

    while (std::getline(ss, value, ',') && col_count < len) {
      try {
        row.push_back(std::stod(value)); // Convert string to double
      } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid number in file at row " << row_count + 1
                  << ", column " << col_count + 1 << std::endl;
        row.push_back(0.0); // Default value for invalid numbers
      }
      col_count++;
    }

    if (int(row.size()) != len) {
      std::cerr << "Warning: Row " << row_count + 1 << " does not have "
                << len << " elements. Filling with zeros." << std::endl;
      while (int(row.size()) < len) {
        row.push_back(0.0);
      }
    }

    if(row_count>0) y.push_back(row);
    row_count++;
  }

  if (row_count < n_wf) {
    n_wf = row_count-1;
    std::cerr << "Warning: File has fewer rows (" << row_count
              << ") than expected (" << n_wf << ")." << std::endl;
  }

  file.close();
  return;
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

void StructuredWaveformSetReader(const std::string fileName,
                                 std::vector<std::vector<double>>& wfs,
                                 const int& daphne_channel,
                                 int& n_wfs) {
  // Open the file
  H5::H5File file(fileName, H5F_ACC_RDONLY);
  int endpoint    = int(daphne_channel/100);
  uint8_t channel = uint8_t(daphne_channel%100);
  std::cout << "ep " << endpoint << " ch " << int(channel) << std::endl;

  // Open datasets
  auto ds_adcs       = file.openDataSet("adcs");
  auto ds_endpoints  = file.openDataSet("endpoints");
  auto ds_channels   = file.openDataSet("channels");

  // Get dimensions of adcs: [n_waveforms, n_samples]
  H5::DataSpace dsp_adcs = ds_adcs.getSpace();
  hsize_t dims[2];
  dsp_adcs.getSimpleExtentDims(dims);
  size_t n_waveforms = dims[0];
  size_t n_samples   = dims[1];

  // Read endpoints and channels
  std::vector<int32_t> endpoints(n_waveforms);
  ds_endpoints.read(endpoints.data(), H5::PredType::NATIVE_INT32);

  std::vector<uint8_t> channels(n_waveforms);
  ds_channels.read(channels.data(), H5::PredType::NATIVE_UINT8);

  // Filter waveform indices
  std::vector<size_t> idx_wfs;
  for (size_t i = 0; i < n_waveforms; ++i) {
    if (endpoints[i] == endpoint && channels[i] == channel)
      idx_wfs.push_back(i);
  }

  int n_chwfs = static_cast<int>(idx_wfs.size());
  if (n_chwfs > n_wfs && n_wfs > 0) n_chwfs = n_wfs;
  else n_wfs = n_chwfs;

  std::cout << "n_waveforms: " << n_waveforms << std::endl;
  std::cout << "n_chwfs: " << n_chwfs << std::endl;
  if (n_chwfs == 0) {
    std::cout << "No waveforms found for channel " << daphne_channel << std::endl;
    // print the available endpoints and channels with no duplicates
    std::cout << "Available endpoints and channels:" << std::endl;
    std::set<std::pair<int32_t, uint8_t>> unique_ep_ch;
    for (size_t i = 0; i < n_waveforms; ++i) {
      unique_ep_ch.insert({endpoints[i], channels[i]});
    }
    for (const auto& ep_ch : unique_ep_ch) {
      std::cout << "Endpoint: " << ep_ch.first << ", Channel: " << static_cast<int>(ep_ch.second) << std::endl;
    }
    wfs.clear();
    return;
  }

  // Read all adcs at once (much faster!)
  std::vector<uint16_t> all_adcs(n_waveforms * n_samples);
  std::cout << "reading adcs" << std::endl;
  ds_adcs.read(all_adcs.data(), H5::PredType::NATIVE_UINT16);
  std::cout << "done reading adcs" << std::endl;
  // Extract only the waveforms we want
  std::cout << "extracting wfs" << std::endl;
  wfs.resize(n_chwfs, std::vector<double>(n_samples));
  for (int i = 0; i < n_chwfs; ++i) {
    size_t row = idx_wfs[i];
    size_t offset = row * n_samples;
    for (size_t j = 0; j < n_samples; ++j) {
      wfs[i][j] = static_cast<double>(all_adcs[offset + j]);
    }
  }
  std::cout << "done extracting wfs" << std::endl;

}

void StructuredEthWaveformSetReader(const std::string fileName,
                                 std::vector<std::vector<double>>& wfs,
                                 const int& daphne_channel,
                                 int& n_wfs) {
  // Open the file
  H5::H5File file(fileName, H5F_ACC_RDONLY);
  int endpoint    = int(daphne_channel/100-1);
  uint8_t channel = uint8_t(daphne_channel%1000);
  std::cout << "ep " << endpoint << " ch " << int(channel) << std::endl;

  // Open datasets
  auto ds_adcs       = file.openDataSet("adcs");
  auto ds_endpoints  = file.openDataSet("endpoints");
  auto ds_channels   = file.openDataSet("channels");

  // Get dimensions of adcs: [n_waveforms, n_samples]
  H5::DataSpace dsp_adcs = ds_adcs.getSpace();
  hsize_t dims[2];
  dsp_adcs.getSimpleExtentDims(dims);
  size_t n_waveforms = dims[0];
  size_t n_samples   = dims[1];

  // Read endpoints and channels
  std::vector<int32_t> endpoints(n_waveforms);
  ds_endpoints.read(endpoints.data(), H5::PredType::NATIVE_INT32);

  std::vector<uint8_t> channels(n_waveforms);
  ds_channels.read(channels.data(), H5::PredType::NATIVE_UINT8);

  // Filter waveform indices
  std::vector<size_t> idx_wfs;
  for (size_t i = 0; i < n_waveforms; ++i) {
    if (endpoints[i] == endpoint && channels[i] == channel)
      idx_wfs.push_back(i);
  }

  int n_chwfs = static_cast<int>(idx_wfs.size());
  if (n_chwfs > n_wfs && n_wfs > 0) n_chwfs = n_wfs;
  else n_wfs = n_chwfs;

  std::cout << "n_waveforms: " << n_waveforms << std::endl;
  std::cout << "n_chwfs: " << n_chwfs << std::endl;
  if (n_chwfs == 0) {
    std::cout << "No waveforms found for channel " << daphne_channel << std::endl;
    // print the available endpoints and channels with no duplicates
    std::cout << "Available endpoints and channels:" << std::endl;
    std::set<std::pair<int32_t, uint8_t>> unique_ep_ch;
    for (size_t i = 0; i < n_waveforms; ++i) {
      unique_ep_ch.insert({endpoints[i], channels[i]});
    }
    for (const auto& ep_ch : unique_ep_ch) {
      std::cout << "Endpoint: " << ep_ch.first << ", Channel: " << static_cast<int>(ep_ch.second) << std::endl;
    }
    wfs.clear();
    return;
  }

  // Read all adcs at once (much faster!)
  std::vector<uint16_t> all_adcs(n_waveforms * n_samples);
  std::cout << "reading adcs" << std::endl;
  ds_adcs.read(all_adcs.data(), H5::PredType::NATIVE_UINT16);
  std::cout << "done reading adcs" << std::endl;
  // Extract only the waveforms we want
  std::cout << "extracting wfs" << std::endl;
  wfs.resize(n_chwfs, std::vector<double>(n_samples));
  for (int i = 0; i < n_chwfs; ++i) {
    size_t row = idx_wfs[i];
    size_t offset = row * n_samples;
    for (size_t j = 0; j < n_samples; ++j) {
      wfs[i][j] = static_cast<double>(all_adcs[offset + j]);
    }
  }
  std::cout << "done extracting wfs" << std::endl;

}

// ---------------------------------------------------------
// --- END METHODS TO READ WAVEFORMS -----------------------
// ---------------------------------------------------------


// ---------------------------------------------------------
// --- METHODS TO READ WHATEVER ----------------------------
// ---------------------------------------------------------

//**********************************************************
std::vector<std::pair<std::string, std::vector<double>>> read_vec_pair_CSV(const std::string& filename) {
//**********************************************************
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  std::vector<std::pair<std::string, std::vector<double>>> data;
  std::string line;

  // Read the header line
  if (std::getline(infile, line)) {
    std::stringstream ss(line);
    std::string header;
    while (std::getline(ss, header, ',')) {
      if (!header.empty()) {
        data.emplace_back(header, std::vector<double>{});
      }
    }
  }

  // Read the data lines
  while (std::getline(infile, line)) {
    std::stringstream ss(line);
    std::string value;
    size_t col_index = 0;

    while (std::getline(ss, value, ',') && col_index < data.size()) {
      try {
        if (!value.empty()) {
          data[col_index].second.push_back(std::stod(value)); // Convert to double
        }
      } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid number at column " << col_index + 1 << std::endl;
        data[col_index].second.push_back(0.0); // Default value for invalid numbers
      }
      col_index++;
    }
  }

  infile.close();
  return data;
}
#endif /* G_Read_hpp */
