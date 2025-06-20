#include "../classe.hpp"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"
#include <cstddef>
#include <vector>


//*********************************************
double complicated_double_expo(double* x, double* p){
//*********************************************
  double amp = p[0];
  double f_f = p[1];
  double t_f = p[2];
  double t_s = p[3];
  double sig = p[4];
  double con = p[5];
  double t_0 = p[6];

  double tf_inv = 1./t_f;
  double ts_inv = 1./t_s;
  double sig_inv = 1./sig;
  double sig_tf = sig*tf_inv;
  double sig_ts = sig*ts_inv;
  double sig2_tf = sig*sig/(2*t_f*t_f);
  double sig2_ts = sig*sig/(2*t_s*t_s);
  double sqrt2_inv = 1.0 / TMath::Sqrt2();

  double result;
  
  double td = t_0-x[0];
  result  = amp * f_f * exp(td*tf_inv + sig2_tf)*TMath::Erfc((td*sig_inv+sig_tf)*sqrt2_inv)*0.5 +
                 amp*(1-f_f)*exp(td*ts_inv + sig2_ts)*TMath::Erfc((td*sig_inv+sig_ts)*sqrt2_inv)*0.5 +
                 con;

  return result;
}

//*********************************************
double* conv_templfft_arraytd(TComplex* templ_fft, double* array_td, size_t len){
//*********************************************
  TComplex array_fft[len];
  TComplex xY[len];
  int len_ = int(len);
  double norm = 1./double(len);
  double* xy;

  compute_r2c_fft(array_fft, array_td, len_);
  freq_domain_convolution(templ_fft, array_fft, xY, len);
  compute_c2r_fft(xy, xY, len);

  for (int j=0; j<len_; j++) xy[j] *= norm; 

  return xy;
}


void cla::Deconvolution(){
  size_t   nsample = memorydepth;
  TComplex templ_fft[nsample], dexp_fft[nsample];
  double   templ_td[nsample], dexp_td[nsample];
  TComplex G[nsample];
  
  std::vector<std::vector<double>> templ_v, avg_muon_v, deco_avg_muon_v;
  std::vector<double> avg_muon_original(memorydepth, 0.0); 
  std::vector<double> sin_muon(memorydepth, 0.0); 
  std::vector<double> time(memorydepth, 0.0); 
  std::vector<double> e_x(memorydepth, tick_len); 
  std::vector<double> e_y(memorydepth, yerr); 

  CompleteWF_Binary(templ_f, templ_v, 1, memorydepth);
  CompleteWF_Binary(muon_f, avg_muon_v, 1, memorydepth);
  Build_Wiener_Filter(&G[0], templ_v[0], n2_);
  FilterAllWF(avg_muon_v, deco_avg_muon_v, &G[0]);
  std::vector<double> deco_avgmuon = deco_avg_muon_v[0];
  // Put the maximum at tick 200
  int max_index = std::distance(deco_avgmuon.begin(), std::max_element(deco_avgmuon.begin(), deco_avgmuon.end())); 
  vector_roll(deco_avgmuon, -max_index+200);
  SubBaseline(deco_avgmuon, 300, 0, deco_avgmuon.size()-300);

  for(int i=0; i<memorydepth; i++) time[i] = double(i)*tick_len;
  TGraphErrors* g_muon = new TGraphErrors(memorydepth, &time[0], &deco_avgmuon[0], 
                                          &e_x[0], &e_y[0]);
  
  double params[6] = {amp, f_fast, tau_fast, tau_slow, sigma, 0.};
  double* par = &params[0];

  TF1 *f1 = new TF1("f1", complicated_double_expo, fit_l, fit_u, 7);
  f1->SetParNames("A", "Fast fraction", "#tau_{fast}", "#tau_{slow}", "#sigma", "c", "t_{0}");
  f1->SetParameters(amp, f_fast, tau_fast, tau_slow, sigma, 0., tick_len*double(roll));
  f1->SetNpx(2000);

  // Scan different t0s
  int best_fit_roll = 0;
  double min_chi2 = 1.e10;
  double fit_amp, fit_f_fast, fit_tau_fast, fit_tau_slow, fit_sigma, fit_t0, fit_c;
  double err_fit_amp, err_fit_f_fast, err_fit_tau_fast, err_fit_tau_slow, err_fit_sigma, err_fit_t0, err_fit_c;
  double best_params[7] = {0.};
  if (!nofit){
    for (size_t idx_fit=0; idx_fit<30; idx_fit++){
      std::vector<double> rolled_avgmuon = deco_avgmuon;
      int fit_roll = idx_fit-15;
      if (idx_fit!=0) vector_roll(rolled_avgmuon, fit_roll); // to start from your best guess
      f1->SetParameters(amp, f_fast, tau_fast, tau_slow, sigma, 0., tick_len*double(fit_roll+roll));
      f1->FixParameter(6, tick_len*double(fit_roll+roll));

      std::cout << "---- Fit " << idx_fit << " ---------" << std::endl;
      TFitResultPtr fit_res = g_muon->Fit(f1, "SR");
      std::cout << "\n" << std::endl;

      if (fit_res->MinFcnValue()<min_chi2 && fit_res==0){
        min_chi2 = fit_res->MinFcnValue();
        best_fit_roll = fit_roll;
        for (size_t i=0; i<7; i++) best_params[i] = fit_res->GetParams()[i];
        fit_amp = fit_res->GetParams()[0]; err_fit_amp = fit_res->GetErrors()[0];
        fit_f_fast = fit_res->GetParams()[1]; err_fit_f_fast = fit_res->GetErrors()[1];
        fit_tau_fast = fit_res->GetParams()[2]; err_fit_tau_fast = fit_res->GetErrors()[2];
        fit_tau_slow = fit_res->GetParams()[3]; err_fit_tau_slow = fit_res->GetErrors()[3];
        fit_sigma = fit_res->GetParams()[4]; err_fit_sigma = fit_res->GetErrors()[4];
        fit_c = fit_res->GetParams()[5]; err_fit_c = fit_res->GetErrors()[5];
        fit_t0 = fit_res->GetParams()[6]; err_fit_t0 = fit_res->GetErrors()[6];
      }


      std::cout << "\n" << std::endl;
    }
    par  = &best_params[0];
    f1->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5], par[6]);
  }
      


  std::cout << "Fit results: " << std::endl;
  std::cout << "Min chi2  = " << min_chi2 << std::endl;
  std::cout << "Ampl      = " << fit_amp << "+-" << err_fit_amp << std::endl;
  std::cout << "Fast frac = " << fit_f_fast << "+-" << err_fit_f_fast << std::endl;
  std::cout << "tau_fast  = " << fit_tau_fast << "+-" << err_fit_tau_fast << std::endl;
  std::cout << "tau_slow  = " << fit_tau_slow << "+-" << err_fit_tau_slow << std::endl;
  std::cout << "sigma     = " << fit_sigma << "+-" << err_fit_sigma << std::endl;
  std::cout << "t0        = " << fit_t0 << "+-" << err_fit_t0 << std::endl;

  TFitResultPtr fit_res = nullptr;
  TCanvas* c_deco = new TCanvas("c_deco","c_deco",0,0,800,600);
  c_deco->cd();
  g_muon->Draw("AP");
  
  if (!nofit) g_muon->Fit(f1, "R");
  else {
    f1->SetLineColor(kBlue);
    f1->Draw("same");
  }
  
  c_deco->Modified(); c_deco->Update();
}
