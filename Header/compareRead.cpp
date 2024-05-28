// compare speed 

//*********************************************
template <typename T>
void SelCalib_WF(vector<vector<T>>& y, vector<vector<T>>& y2, int pre, T sat_low, T sat_up, T bsl){
//*********************************************
  T max_el, min_el;
  size_t len = y[0].size();
  size_t wfs = y.size();
    
  for (size_t i=0; i<wfs; i++) {
    max_el = *max_element( y[i].begin(), y[i].begin()+pre);
    min_el = *min_element( y[i].begin(), y[i].begin()+pre);
    
    if (max_el<bsl && min_el > -bsl) {
      max_el = *max_element( y[i].begin()+pre, y[i].end());
      min_el = *min_element( y[i].begin()+pre, y[i].end());
      
      if (max_el<sat_up && min_el > sat_low) {
        y2.push_back(y[i]);
      }
    }
  }  
  return;
}


//*********************************************
template <typename T>
void SelCalib_WF_resize(vector<vector<T>>& y, vector<vector<T>>& y2, int pre, T sat_low, T sat_up, T bsl){
//*********************************************
  T max_el, min_el;
  size_t len = y[0].size();
  size_t wfs = y.size();
  y2.resize(wfs, vector<T>(len) );
    
  for (size_t i=0; i<wfs; i++) {
    max_el = *max_element( y[i].begin(), y[i].begin()+pre);
    min_el = *min_element( y[i].begin(), y[i].begin()+pre);
    
    if (max_el<bsl && min_el > -bsl) {
      max_el = *max_element( y[i].begin()+pre, y[i].end());
      min_el = *min_element( y[i].begin()+pre, y[i].end());
      
      if (max_el<sat_up && min_el > sat_low) {
        for(size_t j=0; j<len; j++) y2[i][j] = y[i][j];
      }
    }
  }  
  return;
}


void compareRead(){
  std::string fileName = "Ch1.dat";
  std::vector<double> y, y2;
  std::vector<std::vector<double>> yr, yr2;

  int WFs = 10000;
  int len = 5000;

  y.resize(WFs*len, 0);
  y2.resize(WFs*len, 0);
  yr.resize(WFs, std::vector<double>(len));
  yr2.resize(WFs, std::vector<double>(len));


  auto start = std::chrono::high_resolution_clock::now();
  
  for(int i=0; i<10; i++) SelCalib_WF(yr, yr2, len/2, 500., -500., 200.); 
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;
  std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  
  for(int i=0; i<10; i++) SelCalib_WF_resize(yr, yr2, len/2, 500., -500., 200.);  

  end = std::chrono::high_resolution_clock::now();
  duration = end - start;
  std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;



}
