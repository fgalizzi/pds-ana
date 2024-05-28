//
//  G_Utility.hpp
//
//  Created by Federico Galizzi on 03/10/23
//

#ifndef G_Utility_hpp
#define G_Utility_hpp

#include <stdio.h>

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
  for(int i = 0; i < vec.size(); i++) OutFile.write(reinterpret_cast<char*>( &vec[i] ), sizeof(t));
  
  std::cout << "Vector saved in ---> " << fileName << std::endl;
  OutFile.close();
}


//*********************************************
TGraph* Build_CX_Graph (TF1* fgaus, TH1* hI){
//*********************************************
  double norm, mean, sigma;
  int npeaks = fgaus->GetNumberFreeParameters()-3;
  double Pi[npeaks];
  double G = fgaus->GetParameter(1);
 
  TF1* f[npeaks];
 
  for (int i=0; i<npeaks; i++) {
    f[i] = new TF1("ff","[0]*TMath::Gaus(x[0],[1],[2])");
    mean = fgaus->GetParameter(0)+G*i;
    norm = fgaus->GetParameter(4+i);
    sigma = sqrt(pow(fgaus->GetParameter(2),2)+i*pow(fgaus->GetParameter(3),2));
    f[i]->SetParameters(norm, mean, sigma);
    Pi[i] = f[i]->Integral(mean-4*sigma, mean+4*sigma, 1e-5)/(hI->GetBinWidth(5)*hI->GetEntries());
  }

  auto* g_CX = new TGraph(npeaks, Pi);
  return g_CX;
}


// Given an average waveform, it prints the rise time 10%->90% and the fall
// time 90->10. It assumes a flat baseline before the pulse.
//*********************************************
void RiseFallTimeUndershoot(std::vector<double>& waveform, const double tick_len, double& r_time,
    double& f_time, double& undershoot) {
//*********************************************
 
  double* wf = waveform.data();
  auto* g_wf = new TGraph(waveform.size(), wf);
  double x0, x1, xm;

  // Find the maximum amplitude
  double max_amplitude = *std::max_element(waveform.begin(), waveform.end());
  double min_amplitude = *std::min_element(waveform.begin(), waveform.end());
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

#endif /* G_Utility_hpp */
