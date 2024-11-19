#include "../classe.hpp"
#include "Fit/FitResult.h"
#include "Rtypes.h"
#include "TComplex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TVirtualFFT.h"

#include <TF1.h>
#include <cstddef>
#include <iostream>
#include <vector>

#include <Math/Functor.h>
#include "Fit/Fitter.h"

//*********************************************
void double_expo(double* y, const double* p, size_t len, double tick_len){
//*********************************************
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
    double td = t_0-double(i)*tick_len;
    double res = a_f*exp(td*tf_inv + sig2_tf)*TMath::Erfc((td*sig_inv+sig_tf)*sqrt2_inv)*0.5 +
                 a_s*exp(td*ts_inv + sig2_ts)*TMath::Erfc((td*sig_inv+sig_ts)*sqrt2_inv)*0.5 +
                 con;
    // if (isnan(res)==true || res<1.e-40) res = 0.;
    
    y[i] = res;
  }

}


//*********************************************
void Build_FFT(TComplex* G, double* xt, int len){
//*********************************************
  double G_re[len]; double G_im[len];

  TVirtualFFT* fft_r2c = TVirtualFFT::FFT(1, &len, "M R2C");
  fft_r2c->SetPoints(xt);
  fft_r2c->Transform();
  fft_r2c->GetPointsComplex(G_re, G_im);
  
  for (int j=0; j<len*0.5+1; j++) G[j] = TComplex(G_re[j], G_im[j]);
  return;
}

//*********************************************
double* conv_templ_dexp(const double* p,
                       TComplex* templ_fft, size_t len, double tick_len){
//*********************************************
  TComplex dexp_fft[len];
  double   dexp_td[len];
  TComplex xY[len]; double xY_re[len]; double xY_im[len];
  int len_ = int(len);
  double norm = 1./double(len);
  double* xy;

  double_expo(&dexp_td[0], p, len, tick_len);
  Build_FFT(&dexp_fft[0], &dexp_td[0], len_);

  for (int j=0; j<len*0.5+1; j++) {
    xY[j] = templ_fft[j]*dexp_fft[j]; 
    xY_re[j] = xY[j].Re(); xY_im[j] = xY[j].Im();
  }

  TVirtualFFT* fft = TVirtualFFT::FFT(1, &len_, "M C2R");

  fft->SetPointsComplex(xY_re, xY_im);
  fft->Transform();
  xy = fft->GetPointsReal();

  for (int j=0; j<len; j++) xy[j] *= norm; 

  return xy;
}


void cla::Convolution(){

  templ_f = "/Users/federico/PhD/PDE/Templates/Template.dat";
  muon_f = "/Users/federico/PhD/PDE/Muon.dat";
  size_t   nsample = memorydepth;
  TComplex templ_fft[nsample], dexp_fft[nsample];
  double   templ_td[nsample], dexp_td[nsample];
  
  std::vector<std::vector<double>> templ_v, avg_muon_v;
  std::vector<double> avg_muon(memorydepth, 0.0); 
  std::vector<double> sin_muon(memorydepth, 0.0); 
  std::vector<double> time(memorydepth, 0.0); 
  std::vector<double> e_x(memorydepth, 0.001); 
  std::vector<double> e_y(memorydepth, 0.5); 
  
  // Create the template vector
  CompleteWF_Binary(templ_f, templ_v, 1, memorydepth);
  CompleteWF_Binary(muon_f, avg_muon_v, 1, memorydepth);
  
  for(size_t i=0; i<memorydepth; i++){
    templ_td[i] = templ_v[0][i];
    avg_muon[i] = avg_muon_v[0][i];
    time[i] = double(i)*tick_len;
  }


  //rotate vector to adjust the convolution offset
  rotate(avg_muon.begin(), avg_muon.begin()+avg_muon.size()-roll, avg_muon.end());
  Build_FFT(&templ_fft[0], &templ_td[0], memorydepth);

  double params[7] = {a_fast, tau_fast, a_slow, tau_slow, sigma, t_0}; 
  double* par = &params[0];
  double* xy = conv_templ_dexp(par, &templ_fft[0], nsample, tick_len);
  TComplex xY[nsample]; double xY_re[nsample]; double xY_im[nsample];

  TGraphErrors* g_muon = new TGraphErrors(memorydepth, &time[0], &avg_muon[0],
                                         &e_x[0], &e_y[0]);

  auto chi2Function = [&](const double* par){
    int np = g_muon->GetN();
    double chi2 = 0.;
    double* x = g_muon->GetX();
    double* y = g_muon->GetY();
    double norm = 1./double(memorydepth);

    double_expo(&dexp_td[0], par, nsample, tick_len);
    Build_FFT(&dexp_fft[0], &dexp_td[0], memorydepth);

    for (int j=0; j<memorydepth*0.5+1; j++) {
      xY[j] = templ_fft[j]*dexp_fft[j]; 
      xY_re[j] = xY[j].Re(); xY_im[j] = xY[j].Im();
    }

    TVirtualFFT* fft = TVirtualFFT::FFT(1, &memorydepth, "M C2R");

    fft->SetPointsComplex(xY_re, xY_im);
    fft->Transform();
    xy = fft->GetPointsReal();
    for (int j=0; j<memorydepth; j++) xy[j] *= norm; 

    int int_fit_l = int(fit_l/tick_len);
    int int_fit_u = min(int(fit_u/tick_len), np);

    for (int i=int_fit_l; i<int_fit_u; i++){
      double d = y[i]-xy[i];
      double err = g_muon->GetErrorY(i);
      chi2 += d*d/(err*err);
    }
    return chi2/double(int_fit_u-int_fit_l);
  };

  ROOT::Math::Functor fnc(chi2Function, 7);
  ROOT::Fit::Fitter fitter;
  fitter.SetFCN(fnc, par);
  fitter.Config().ParSettings(0).SetName("A_{s}");
  fitter.Config().ParSettings(1).SetName("#tau_{s}");
  fitter.Config().ParSettings(2).SetName("A_{t}");
  fitter.Config().ParSettings(3).SetName("#tau_{t}");
  fitter.Config().ParSettings(4).SetName("#sigma");
  fitter.Config().ParSettings(5).SetName("t_{0}");
  fitter.Config().ParSettings(6).SetName("c");

  bool ok;
  if(!no_fit){
    ok = fitter.FitFCN();
    const ROOT::Fit::FitResult &result = fitter.Result();
    result.Print(std::cout);
  }

  for(size_t i=0; i<memorydepth; i++) sin_muon[i]=xy[i];
  TGraph* g_sint = new TGraph(memorydepth, &time[0], &sin_muon[0]);
  
  
  double A_f   = par[0];
  double tau_f = par[1];
  double A_s   = par[2];
  double tau_s = par[3];
  std::cout <<"\n \nFast/Slow intensities" << std::endl;
  std::cout <<"Is = " << A_f*tau_f/(A_f*tau_f+A_s*tau_s) << std::endl;  
  std::cout <<"It = " << A_s*tau_s/(A_f*tau_f+A_s*tau_s) << std::endl;  





  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,800);
  c2->cd();
  g_muon->GetXaxis()->SetTitle("Time #mus");
  g_muon->GetYaxis()->SetTitle("Amplitude [a.u.]");
  // g_muon->Set;

  g_sint->SetLineWidth(2);
  g_sint->SetLineColor(kRed);
  g_muon->Draw();
  g_sint->Draw("same");
  c2->Modified(); c2->Update();

  std::cout << par[0]<< std::endl;
  return;
}