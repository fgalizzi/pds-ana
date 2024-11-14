#include "../classe.hpp"
#include "Rtypes.h"
#include "TComplex.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <TF1.h>
#include <cstddef>
#include <vector>

void double_expo(double* x, const double* p, size_t len){
  double a_f = p[0];
  double t_f = p[1];
  double a_s = p[2];
  double t_s = p[3];
  double sig = p[4];
  double t_0 = p[5];
  double con = p[6];

  double tf_inv = 1./t_f;
  double ts_inv = 1./t_s;
  double sig_inv = 1./sig;
  double sig_tf = sig*tf_inv;
  double sig_ts = sig*ts_inv;
  double sig2_tf = sig*sig/(2*t_f*t_f);
  double sig2_ts = sig*sig/(2*t_s*t_s);
  double sqrt2_inv = 1.0 / TMath::Sqrt2();

  for(size_t i=0; i<len; i++){
    double td = t_0-double(i);
    double res = a_f*exp(td*tf_inv + sig2_tf)*TMath::Erfc((td*sig_inv+sig_tf)*sqrt2_inv)*0.5 +
                 a_s*exp(td*ts_inv + sig2_ts)*TMath::Erfc((td*sig_inv+sig_ts)*sqrt2_inv)*0.5 +
                 con;
    // if (isnan(res)==true || res<1.e-40) res = 0.;
    
    x[i] = res;
  }

}


//*********************************************
void Build_FFT(TComplex* G, double* xt, int len){
//*********************************************
  double G_re[len]; double G_im[len];
  TVirtualFFT* fft = TVirtualFFT::FFT(1, &len, "M R2C");
  
  fft->SetPoints(xt);
  fft->Transform();
  fft->GetPointsComplex(G_re, G_im);

  for (int j=0; j<len*0.5+1; j++) G[j] = TComplex(G_re[j], G_im[j]);
  return;
}

double conv_templ_dexp(const double* x, const double* p,
                       TComplex* templ_fft, size_t len){
  
  TComplex dexp_fft[len];
  double   dexp_td[len];
  double* xy;
  TComplex xY[len]; double xY_re[len]; double xY_im[len];
  int len_ = int(len);
  
  double_expo(&dexp_td[0], p, len);
  Build_FFT(&dexp_fft[0], &dexp_td[0], len_);

  for (int j=0; j<len*0.5+1; j++) {
    xY[j] = templ_fft[j]*dexp_fft[j]; 
    xY_re[j] = xY[j].Re(); xY_im[j] = xY[j].Im();
  }

  TVirtualFFT* fft = TVirtualFFT::FFT(1, &len_, "M C2R");

  fft->SetPointsComplex(xY_re, xY_im);
  fft->Transform();
  xy = fft->GetPointsReal();


  return xy[int(x[0])];
}


void cla::Convolution(){

  templ_f = "/Users/federico/PhD/PDE/Templates/Template.dat";
  muon_f = "/Users/federico/PhD/PDE/Muon.dat";
  size_t   nsample = memorydepth;
  TComplex templ_fft[nsample];
  double   templ_td[nsample];
  
  std::vector<std::vector<double>> templ_v, avg_muon_v;
  std::vector<double> templ(memorydepth, 0.0); 
  std::vector<double> avg_muon(memorydepth, 0.0); 
  std::vector<double> sin_muon(memorydepth, 0.0); 
  std::vector<double> time(memorydepth, 0.0); 
  std::vector<double> e_x(memorydepth, 1.0); 
  std::vector<double> e_y(memorydepth, 10.0); 
  // Create the template vector
  CompleteWF_Binary(templ_f, templ_v, 1, memorydepth);
  CompleteWF_Binary(muon_f, avg_muon_v, 1, memorydepth);
  
  for(size_t i=0; i<memorydepth; i++){
    templ_td[i] = templ_v[0][i];
    avg_muon[i] = avg_muon_v[0][i];
    time[i] = double(i);
  }




  rotate(avg_muon.begin(), avg_muon.begin()+avg_muon.size()-70, avg_muon.end());
  Build_FFT(&templ_fft[0], &templ_td[0], memorydepth);

  double params[7] = {a_fast, tau_fast/tick_len, a_slow, tau_slow/tick_len, sigma/tick_len, t_0/tick_len}; 
  
  TF1* fitFunc = new TF1("fitFunc",
                         [&](double* x, double* p) {
                            return conv_templ_dexp(x, p, &templ_fft[0], nsample);
                         },
                         1400., 2200., 7);
 
  std::cout << fitFunc->GetNumberFitPoints() << std::endl;
  fitFunc->SetParameters(a_fast, tau_fast/tick_len, a_slow, tau_slow/tick_len, sigma/tick_len, t_0/tick_len);
  fitFunc->SetParNames("A_{s}", "#tau_{s}", "A_{t}", "#tau_{t}", "#sigma", "t_{0}", "c");
  fitFunc->SetNpx(2000);
  fitFunc->SetLineColor(kBlack);

  
  auto start = std::chrono::high_resolution_clock::now();
  for(size_t i=0; i<1e4; i++){
    fitFunc->Eval(50.);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;
  std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
  // return;
  

  TGraphErrors* g_muon = new TGraphErrors(memorydepth, &time[0], &avg_muon[0],
                                         &e_x[0], &e_y[0]);
  std::cout << fitFunc->GetNumberFitPoints() << std::endl;
  // TFitter::SetMaxIterations(10);
  // g_muon->Fit(fitFunc, "RSN");
  for(size_t i=0; i<memorydepth; i++) sin_muon[i]=fitFunc->Eval(i);
  TGraph* g_sint = new TGraph(memorydepth, &time[0], &sin_muon[0]);
  
  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,800);
  c2->cd();
  g_sint->SetLineColor(kRed);
  g_sint->Draw();
  g_muon->Draw("same");
  c2->Modified(); c2->Update();


  return;
}
