//
//  G_Read.hpp
//
//  Created by Federico Galizzi on 28/09/23
//

#ifndef G_Read_hpp
#define G_Read_hpp





#ifndef hdf5torootclass_cxx
#include "../Class/ProtoduneHD/wffunctions2.h"
#include "../Class/ProtoduneHD/hdf5torootclass.h"
#endif // !hdf5torootclass_cxx
using namespace std;


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
void PDHD_ch_wfs(std::string fileName, vector<vector<double>>& y, int this_ch, int WFs){
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

       for (size_t ievt=0; ievt<nentries && wf_counter<WFs; ievt++){ // loop over entries in root file

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

// Read CSV file with #WF=WFs len-tick long
//**********************************************************
void CSV_double_WF_Binary(std::string filename, vector<vector<double>>& y, int& n_wf, int len){
//**********************************************************
  std::ifstream file(filename);

  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  std::string line;
  size_t row_count = 0;

  while (std::getline(file, line) && row_count < n_wf) {
    std::vector<double> row;
    std::stringstream ss(line);
    std::string value;
    size_t col_count = 0;

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

    if (row.size() != len) {
      std::cerr << "Warning: Row " << row_count + 1 << " does not have "
                << len << " elements. Filling with zeros." << std::endl;
      while (row.size() < len) {
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



#endif /* G_Read_hpp */
