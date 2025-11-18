//
//  G_Utility.hpp
//
//  Created by Federico Galizzi on 03/10/23
//

#include "TComplex.h"
#include "TVirtualFFT.h"
#include <cstddef>
#ifndef G_Utility_hpp
  #include <stdio.h>
  #include <iostream>
  #include <fstream>
  #include <utility>
  #include <vector>

  #include <TF1.h>
  #include <TMatrixD.h>
  #include "TFitResult.h"
  #include <TH1.h>
  #include <TGraphErrors.h>
  #include <TGraph.h>
#define G_Utility_hpp


using namespace std;

//*********************************************
template <typename T>
  std::vector<T> linspace(T a, T b, size_t N) {
  T h = (b - a) / static_cast<T>(N-1);
  std::vector<T> xs(N);
  typename std::vector<T>::iterator x;
  T val;
  
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
    *x = val;
  return xs;     // delta function
 }

// Compute the integral of a TGraph object
//*********************************************
double g_integral(TGraph* g, double x0, double x1) {
//*********************************************
  auto fc_gintegral = [g](double *x, double* p) {
    return p[0]*g->Eval(x[0]);
  };
  TF1 f("f", fc_gintegral, x0, x1, 1);
  f.SetNpx(g->GetN());
  f.SetParameter(0, 1.);
  return f.Integral(x0, x1,1e-4);
}

// Scale all the values of a TGraph
void g_scale(TGraph* g, double c) {
  for (int i=0; i<g->GetN(); i++) {
    g->GetY()[i] *= c;
  }
  return;
}


// Find x where g(x)=y in range [x0;x1]
//*********************************************
inline double g_find_x(TGraph* g, double y, double x0, double x1, double
    epsilon = 1e-2) {
//*********************************************
  auto fc_g= [g](double *x, double* p) {
    return p[0]*g->Eval(x[0]);
  };
  TF1 f("f", fc_g, x0, x1, 1);
  f.SetNpx(g->GetN());
  f.SetParameter(0, 1.);

  return f.GetX(y, x0, x1, epsilon);
}


//*********************************************
void VecDouble_in_Binary(std::string fileName, std::vector<double>& vec){
//*********************************************
  double t;
  std::ofstream OutFile (fileName, ios::binary);
  for(size_t i = 0; i < vec.size(); i++) OutFile.write(reinterpret_cast<char*>( &vec[i] ), sizeof(t));
  
  std::cout << "Vector saved in ---> " << fileName << std::endl;
  OutFile.close();
}

// Error propagation given a TFitResultPtr and a TF1 object
//*********************************************
double error_propagation(TFitResultPtr fit_res, TF1* fit_func, int idx_par1,
                         int idx_par2, std::string operation){
//*********************************************
  double result = 0.;
  double p1 = fit_func->GetParameter(idx_par1);
  double p2 = fit_func->GetParameter(idx_par2);
  double e1 = fit_func->GetParError(idx_par1);
  double e2 = fit_func->GetParError(idx_par2);
  TMatrixD cov = fit_res->GetCovarianceMatrix();
  double cov_p1_p2 = cov(idx_par1, idx_par2);
  
  if (operation == "sum") 
    result = sqrt(e1 * e1 + e2 * e2 + 2 * cov_p1_p2);
  else if (operation == "sub") 
    result = sqrt(e1 * e1 + e2 * e2 - 2 * cov_p1_p2);
  else if (operation == "mul") 
    result = sqrt((e1 * e1) * (p2 * p2) + (e2 * e2) * (p1 * p1) + 2 * p1 * p2 * cov_p1_p2);
  else if (operation == "div") 
    result = sqrt((e1 * e1) / (p2 * p2) + (e2 * e2) * (p1 * p1) / (p2 * p2 * p2 * p2) 
                  - 2 * p1 * cov_p1_p2 / (p2 * p2 * p2));
  else
    throw std::invalid_argument("Invalid operation: must be 'sum', 'sub', 'mul', or 'div'.");

  return result;
}

// Error propagation given two parameters and their errors (no covariance)
//*********************************************
double error_propagation(double p1, double e1, double p2, double e2, std::string operation){
//*********************************************
  double result = 0.;
  if (operation == "sum") 
    result = sqrt(e1 * e1 + e2 * e2);
  else if (operation == "sub") 
    result = sqrt(e1 * e1 + e2 * e2);
  else if (operation == "mul") 
    result = sqrt((e1 * e1) * (p2 * p2) + (e2 * e2) * (p1 * p1));
  else if (operation == "div") 
    result = sqrt((e1 * e1) / (p2 * p2) + (e2 * e2) * (p1 * p1) / (p2 * p2 * p2 * p2));
  else if (operation == "sqrt_sum"){
    double f2 = p1*p1 + p2*p2;
    result    = sqrt( (e1*e1 * p1*p1)/f2 + (e2*e2 * p2*p2)/f2 );
  }
  else if (operation == "sqrt_sub"){
    double f2 = p1*p1 - p2*p2;
    result    = sqrt( (e1*e1 * p1*p1)/f2 + (e2*e2 * p2*p2)/f2 );
  }
  else
    throw std::invalid_argument("Invalid operation: must be 'sum', 'sub', 'mul', 'div', 'sqrt_sum' or 'sqrt_sub'.");

  return result;
}

//*********************************************
// Probably, this is overestimating the error
TGraphErrors* Build_CX_Graph_Cov (TF1* fgaus, TH1* hI, TFitResultPtr FitRes, double& avg_n_photons){
//*********************************************
  int npeaks = int(fgaus->GetXmax()/fgaus->GetParameter(1)+1); // only fitted peaks
  vector<double> Pi(npeaks, 0.0);
  vector<double> Err_Pi(npeaks, 0.0);
  vector<double> X(npeaks, 0.0);
  vector<double> Err_X(npeaks, 0.0);
  double G = fgaus->GetParameter(1);
  double mu_0 = fgaus->GetParameter(0);
 
  for (int peak=0; peak<npeaks; peak++) {
    TF1 mock_fgaus = *fgaus;
    for(int i = 0; i<fgaus->GetNumberFreeParameters(); i++){
      if(i>3 && i!=peak+4) mock_fgaus.SetParameter(i, 0.);
    }
 
    double sigma = sqrt(pow(fgaus->GetParameter(2),2)+peak*pow(fgaus->GetParameter(3),2));
    double ampl  = fgaus->GetParameter(4+peak);
    double area  = ampl*sigma*sqrt(2*TMath::Pi());

    Pi[peak] = area/(hI->GetBinWidth(5)*hI->GetEntries());
    Err_Pi[peak] = mock_fgaus.IntegralError(mu_0+(peak-4)*G, mu_0+(peak+4)*G, FitRes->GetParams(), 
                        FitRes->GetCovarianceMatrix().GetMatrixArray())/(hI->GetBinWidth(5)*hI->GetEntries());
    X[peak] = peak;
    
    if(peak==0) avg_n_photons = -log(Pi[peak]);
  }

  std::cout << "\n\nYOU ARE OVERESTIMATING THE ERROR??" << std::endl;
  TGraphErrors* g_CX = new TGraphErrors(Pi.size(), &X[0], &Pi[0], &Err_X[0], &Err_Pi[0]);
  return g_CX;
}

//***********************************************
TGraphErrors* Build_CX_Graph(TF1* fgaus, TH1* hI, double& avg_n_photons){
//*********************************************
  int npeaks = int(fgaus->GetXmax()/fgaus->GetParameter(1)+1); // only fitted peaks
  vector<double> Pi(npeaks, 0.0);
  vector<double> Err_Pi(npeaks, 0.0);
  vector<double> X(npeaks, 0.0);
  vector<double> Err_X(npeaks, 0.0);
  double G = fgaus->GetParameter(1);
  double mu_0 = fgaus->GetParameter(0);
 
  for (int peak=0; peak<npeaks; peak++) {
    double sigma0 = fgaus->GetParameter(2);
    double sigmac = fgaus->GetParameter(3);
    double sigma  = sqrt(pow(sigma0,2)+peak*pow(sigmac,2));
    double err_sigma0 = fgaus->GetParError(2);
    double err_sigmac = fgaus->GetParError(3);
    double err_sigma  = error_propagation(sigma0, err_sigma0, sqrt(peak)*sigmac, sqrt(peak)*err_sigmac, "sqrt_sum"); 
    double ampl  = fgaus->GetParameter(4+peak);
    double err_ampl = fgaus->GetParError(4+peak);
    double area  = ampl*sigma*sqrt(2*TMath::Pi());

    Pi[peak] = area/(hI->GetBinWidth(5)*hI->GetEntries());
    double err_stat = sqrt(area)/(hI->GetBinWidth(5)*hI->GetEntries());
    double err_fit  = error_propagation(ampl, err_ampl, sigma, err_sigma, "mul")*sqrt(2*TMath::Pi())/(hI->GetBinWidth(5)*hI->GetEntries());
    Err_Pi[peak] = sqrt( err_stat*err_stat + err_fit*err_fit );
    X[peak] = peak;
    
    if(peak==0) avg_n_photons = -log(Pi[peak]);
  }

  TGraphErrors* g_CX = new TGraphErrors(Pi.size(), &X[0], &Pi[0], &Err_X[0], &Err_Pi[0]);
  return g_CX;
}


// Given an average waveform, it prints the rise time 10%->90% and the fall
// time 90->10. It assumes a flat baseline before the pulse.
//*********************************************
void RiseFallTimeUndershoot(std::vector<double>& waveform, const double& tick_len, int& int_up){
//*********************************************
  double* wf = waveform.data();
  auto* g_wf = new TGraph(waveform.size(), wf);
  double x0, x1, xm, undershoot, r_time, f_time;

  if (int_up < 0 || int_up >= waveform.size()) int_up = waveform.size()-1;

  // Find the maximum amplitude
  double max_amplitude = *std::max_element(waveform.begin(), waveform.begin()+int_up);
  double min_amplitude = *std::min_element(waveform.begin()+int_up, waveform.end());
  undershoot = -min_amplitude/max_amplitude*100;

  // Calculate the threshold levels (10% and 90% of max amplitude)
  double threshold_10 = 0.1 * max_amplitude;
  double threshold_90 = 0.9 * max_amplitude;

  // Find the rise time
  auto iter_start = std::find_if(waveform.begin(), waveform.end(), [threshold_10](double value) {
    return value >= threshold_10; });

  auto iter_end = std::find_if(iter_start, waveform.end(), [threshold_90](double value) {
    return value >= threshold_90; });

  r_time = std::distance(iter_start, iter_end)*tick_len*1000;

  // Find the fall time
  iter_start = std::find_if(iter_end, waveform.end(), [max_amplitude](double value) {
    return value == max_amplitude; });

  iter_start = std::find_if(iter_start, waveform.end(), [threshold_90](double value){
    return value <= threshold_90; });

  iter_end = std::find_if(iter_start, waveform.end(), [threshold_10](double value){
    return value <= threshold_10; });

  f_time = std::distance(iter_start, iter_end)*tick_len*1000;
 
  std::cout << "Max " << max_amplitude << std::endl;
  std::cout << "Min " << min_amplitude << std::endl;
  std::cout << "Rise time 10%->90% [ns] " << r_time << "\nFall time 90%->10% [ns] " << f_time << std::endl;

  xm = g_find_x(g_wf, max_amplitude*0.98, 0., double(waveform.size()), 0.1);
  x0 = g_find_x(g_wf, threshold_10, 0., xm, 0.1);
  x1 = g_find_x(g_wf, threshold_90, x0, xm, 0.1);
  r_time = (x1-x0)*tick_len*1000;
  x0 = g_find_x(g_wf, threshold_90, xm, double(waveform.size()), 0.1);
  x1 = g_find_x(g_wf, threshold_10, x0, double(waveform.size()), 0.1);
  f_time = (x1-x0)*tick_len*1000;
  
  std::cout << "\n \nNew method " << std::endl;
  std::cout << "Rise time 10%->90% [ns] " << r_time << "\nFall time 90%->10% [ns] " << f_time << std::endl;
  std::cout << "Undershoot [%]: " << undershoot << std::endl;

}

//*********************************************
void min_max_element (vector<vector<double>>& y, double& ymin, double& ymax){
//*********************************************
  ymin = std::numeric_limits<double>::max();
  ymax = std::numeric_limits<double>::min();

  double min, max;

  for (auto vec : y){
    min = *min_element(std::begin(vec), std::end(vec));
    max = *max_element(std::begin(vec), std::end(vec));
    if (max > ymax) ymax = max;
    if (min < ymin) ymin = min;
  }

}

// Two functions to print the tuple's content easily
//*********************************************
template<typename Tuple, std::size_t... Is>
//*********************************************
void print_tuple_impl(const Tuple& t, std::index_sequence<Is...>) {
  ((std::cout << std::get<Is>(t) << "\t"), ...);
  std::cout << std::endl;
}
template<typename... Ts>
void print_tuple(const std::tuple<Ts...>& t) {
    print_tuple_impl(t, std::index_sequence_for<Ts...>{});
}

// Get the name of the current directory 
//*********************************************
std::string This_Directory_Name(){
//*********************************************
  try{
    std::filesystem::path currentPath = std::filesystem::current_path();
    std::string currentDirName = currentPath.filename().string();
  return currentDirName;
  } 
  catch (const std::filesystem::filesystem_error& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return "DIRECTORY NOT FOUND !!";
}

// Take the most frequent value of a vector
//*********************************************
template<typename T>
T Vector_MPV(std::vector<T> vec){
//*********************************************
  std::unordered_map<T, int> my_map;

  // Count the occurrences of each value in the vector
  for (T num : vec) my_map[num]++;
  
  // Find the most frequent value
  int maxFrequency = 0;
  T mostFrequentValue = vec[0]; // Initialize with the first value
  for (const auto& pair : my_map) {
      if (pair.second > maxFrequency) {
          maxFrequency = pair.second;
          mostFrequentValue = pair.first;
      }
  } 
  return mostFrequentValue; 
}

//Print the keys as header and the values in columns in a .csv file
//*********************************************
void print_vec_pair_csv(std::string filename, std::vector<std::pair<std::string,double>> my_map,
                        std::string comment=""){
//*********************************************
  std::vector<std::pair<std::string,double>>::iterator it;
  // Check if the file already exists
  std::ifstream infile(filename);
  bool file_exists = infile.good();
  infile.close();
  // Open the file in append mode
  std::ofstream outfile;
  outfile.open(filename, std::ios_base::app); // Append to file if exists
  // If the file doesn't exist, add the header row
  if (!file_exists) {
    for(it = my_map.begin(); it!=my_map.end(); it++) outfile << it->first << ",";
    outfile << std::endl;
  }
  // Append new data lines
  for(it = my_map.begin(); it!=my_map.end(); it++) outfile << it->second << ",";
  if (comment!="n") outfile << comment << ",";
  
  outfile << std::endl;
  // Close the file
  outfile.close();
  return;
}


//*********************************************
int extract_channel_from_filename(std::string filename) {
//*********************************************
  // Find the position of the first underscore and the dot
  std::size_t underscore_pos = filename.find('_');
  std::size_t dot_pos = filename.find('.');
  int ch = -999;
  // Extract the number part between the underscore and dot
  if (underscore_pos != std::string::npos && dot_pos != std::string::npos) {
    std::string number_str = filename.substr(underscore_pos + 1, dot_pos - underscore_pos - 1);
    // Convert the extracted string to an integer
    ch = std::stoi(number_str);
  }
  else std::cout << "Invalid channel name" << std::endl;

 return ch;
}

// Rotate vector in numpy.roll fashion
//*********************************************
template<typename T>
void vector_roll(std::vector<T>& vec, int roll){
//*********************************************
  if (roll > 0) std::rotate(vec.begin(), vec.begin()+vec.size()-roll, vec.end());
  if (roll < 0) std::rotate(vec.begin(), vec.begin()-roll, vec.end());
  return;
}

//*********************************************
void allign_wfs(vector<vector<double>>& waveforms, const int x_half_height) {
//*********************************************
 
  for (auto& waveform : waveforms){
    double max_amplitude = *std::max_element(waveform.begin(), waveform.end());
    double* wf = waveform.data();
    TGraph* g_wf = new TGraph(waveform.size(), wf);
    std::cout << "find" << std::endl;
    int x_half = (int) g_find_x(g_wf, 0.5*max_amplitude, 0., double(waveform.size()), 0.1);
    std::cout << "roll" << std::endl;
    vector_roll(waveform, x_half-x_half_height);
    std::cout << "end \n\n" << std::endl;
  }

}



// Compute the FFT of a real (time domain) signal
//*********************************************
void compute_r2c_fft(TComplex* c_fft, double* xt, int len){
//*********************************************
  double fft_re[len]; double fft_im[len];

  TVirtualFFT* fft_r2c = TVirtualFFT::FFT(1, &len, "M R2C");
  fft_r2c->SetPoints(xt);
  fft_r2c->Transform();
  fft_r2c->GetPointsComplex(fft_re, fft_im);
  
  for (int j=0; j<len*0.5+1; j++) c_fft[j] = TComplex(fft_re[j], fft_im[j]);
  return;
}

// Compute the FFT of a complex (frequency domain) signal
//*********************************************
void compute_c2r_fft(double* xt, TComplex* c_fft, int len){
//*********************************************
  double* xy;
  double fft_re[size_t(len*0.5+1)]; double fft_im[size_t(0.5*len+1)];
  for (int j=0; j<len*0.5+1; j++) {
    fft_re[j] = c_fft[j].Re();
    fft_im[j] = c_fft[j].Im();
  }

  TVirtualFFT* fft_c2r = TVirtualFFT::FFT(1, &len, "M C2R");
  fft_c2r->SetPointsComplex(fft_re, fft_im);
  fft_c2r->Transform();
  xy = fft_c2r->GetPointsReal();
  for (int j=0; j<len; j++) xt[j] = xy[j];
  return;
}

// Frequency domain convolution
//*********************************************
void freq_domain_convolution(TComplex* c_fft1, TComplex* c_fft2, TComplex* c_fft_out, const int& len){
//*********************************************
  for (int j=0; j<len*0.5+1; j++) {
    c_fft_out[j] = c_fft1[j]*c_fft2[j];
  }
  return;
}

//*********************************************
void Build_Matched_Filter(TComplex* G, vector<double> t_template){
//*********************************************
  int len = t_template.size();

  double xt[len];
  double G_re[len]; double G_im[len];
  TVirtualFFT* fft = TVirtualFFT::FFT(1, &len, "M R2C");
  
  for (int j=0; j<len; j++) xt[j] = t_template[len-j-1];

  fft = TVirtualFFT::FFT(1, &len, "M R2C");
  fft->SetPoints(xt);
  fft->Transform();
  fft->GetPointsComplex(G_re, G_im);


  for (int j=0; j<len*0.5+1; j++) G[j] = TComplex(G_re[j], G_im[j]);
  
}

//*********************************************
void Build_Wiener_Filter(TComplex* G, vector<double>& t_template, double n2){
//*********************************************
  size_t len = t_template.size();
  std::vector<double> xs_v(len, 0.0); xs_v[1] = 1.;
  TComplex templ[len]; compute_r2c_fft(templ, t_template.data(), len);
  TComplex xS[len];    compute_r2c_fft(xS, xs_v.data(), len);
  for (size_t i=0; i<len; i++){
    G[i] = TComplex::Conjugate(templ[i])*xS[i].Rho2() / (templ[i].Rho2()*xS[i].Rho2() + n2);
  }
}

//*********************************************
template<typename T>
void vectorVector_to_vector(vector<T>& vec1, vector<vector<T>>& vec2) {
//*********************************************
  vec1.clear();
  for (const auto& vec : vec2) {
    for (const auto& val : vec) {
      vec1.push_back(val);
    }
  }
}


#endif /* G_Utility_hpp */
