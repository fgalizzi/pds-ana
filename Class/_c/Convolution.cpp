#include "../classe.hpp"

// ****************************************************************
// Macro to perform the convolution of the spe response with the
// liquid argon scintillation time profile to fit the muon signal.
// How to use:
// 1- set nofit=1 to see how good are the fit parameters you set
//    Parameters to set are: a_fast, tau_fast, a_slow, tau_slow,
//    roll, fit_l and fit_u.
//    The roll parameter is used to adjust the convolution offset
//    The fit_l and fit_u parameters are used to set the fit range
//    The a_fast, tau_fast, a_slow, tau_slow are the amplitudes 
//    and the tau of fast and slow components.
// 2- set nofit=0 to perform the fit
// ****************************************************************

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////
string pde_result_file = "/eos/home-f/fegalizz/PDE_MiB/PDE_Results/Conv_results/res.csv";
///////////////////////////////////////////////////////////////////

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

  double tf_inv = 1./t_f;
  double ts_inv = 1./t_s;

  for(size_t i=0; i<len; i++){
    double td = -double(i)*tick_len;
    if (i==0) y[i] = a_f + a_s;
    else      y[i] = a_f*exp(td*tf_inv) + a_s*exp(td*ts_inv);
    // if (i==0) y[i] = a_f + a_s;
    // else      y[i] = a_f*exp(td*tf_inv) + a_s*exp(td*ts_inv);
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
double* conv_templ_dexp(const double* p, TComplex* templ_fft, size_t len, double tick_len){
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


//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Convolution(){

  size_t   nsample = memorydepth;
  TComplex templ_fft[nsample], dexp_fft[nsample];
  double   templ_td[nsample], dexp_td[nsample];
  
  std::vector<std::vector<double>> templ_v, avg_muon_v;
  std::vector<double> avg_muon_original(memorydepth, 0.0); 
  std::vector<double> avg_muon(memorydepth, 0.0); 
  std::vector<double> sin_muon(memorydepth, 0.0); 
  std::vector<double> time(memorydepth, 0.0); 
  std::vector<double> e_x(memorydepth, tick_len); 
  std::vector<double> e_y(memorydepth, yerr); 
 

  // Create the template and average muon vectors
  CompleteWF_Binary(templ_f, templ_v, 1, memorydepth);
  CompleteWF_Binary(muon_f, avg_muon_v, 1, memorydepth);
  
  for(size_t i=0; i<memorydepth; i++){
    templ_td[i] = templ_v[0][i];
    avg_muon_original[i] = avg_muon_v[0][i];
    time[i] = double(i)*tick_len;
  }


  // Rotate vector to adjust the convolution offset. Need a good guess!
  vector_roll(avg_muon_original, roll);
  avg_muon = avg_muon_original;
  Build_FFT(&templ_fft[0], &templ_td[0], memorydepth);

  double params[4] = {a_fast, tau_fast, a_slow, tau_slow}; 
  double* par = &params[0];
  double* xy = conv_templ_dexp(par, &templ_fft[0], nsample, tick_len);
  TComplex xY[nsample]; double xY_re[nsample]; double xY_im[nsample];


  // ---- FIT -------------------------------------------------------- 
  int best_fit_roll = 0;
  if(!nofit){
    bool fit_check;
    double best_params[4];
    double min_chi2 = 1.e10;
    double fit_a_fast, err_a_fast;  
    double fit_tau_fast, err_tau_fast;
    double fit_a_slow, err_a_slow;  
    double fit_tau_slow, err_tau_slow;

 
    // Scan different t0s
    for (size_t idx_fit=0; idx_fit<30; idx_fit++){
      avg_muon = avg_muon_original;
      int fit_roll = idx_fit-15;
      if (idx_fit!=0) vector_roll(avg_muon, fit_roll); // to start from your best guess
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

        return chi2;
      };

      int int_fit_l = int(fit_l/tick_len);
      int int_fit_u = int(fit_u/tick_len);
      ROOT::Math::Functor fnc(chi2Function, 4);
      ROOT::Fit::Fitter fitter;
      fitter.SetFCN(fnc, par);
      fitter.Config().ParSettings(0).SetName("A_{f}");
      fitter.Config().ParSettings(1).SetName("#tau_{f}");
      fitter.Config().ParSettings(2).SetName("A_{s}");
      fitter.Config().ParSettings(3).SetName("#tau_{s}");

      std::cout << "---- Fit " << idx_fit << " ---------" << std::endl;
      fit_check = fitter.FitFCN();
      fitter.SetNumberOfFitPoints(static_cast<size_t>(int_fit_u-int_fit_l));
      auto result = fitter.Result();
      result.Print(std::cout);
      std::cout << "\n" << std::endl;

      if (result.MinFcnValue()<min_chi2 && fit_check==true){
        std::cout << result.GetParams()[0] << " " << result.GetParams()[0] << std::endl;
        min_chi2 = result.MinFcnValue();
        best_fit_roll = fit_roll;
        for (size_t i=0; i<4; i++) best_params[i] = result.GetParams()[i];
        fit_a_fast   = result.GetParams()[0]; err_a_fast   = result.GetErrors()[0];
        fit_tau_fast = result.GetParams()[1]; err_tau_fast = result.GetErrors()[1];
        fit_a_slow   = result.GetParams()[2]; err_a_slow   = result.GetErrors()[2];
        fit_tau_slow = result.GetParams()[3]; err_tau_slow = result.GetErrors()[3];
      }
      
      par = &best_params[0];
      xy  = conv_templ_dexp(par, &templ_fft[0], nsample, tick_len);
      std::cout << "\n\n\n" << std::endl;
    }
 

    // Update class members with the best fit parameters
    a_fast   = fit_a_fast;  
    tau_fast = fit_tau_fast;
    a_slow   = fit_a_slow;  
    tau_slow = fit_tau_slow;
   

    if (print == true){
      vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 
      double date, electronic;
      string comment;

      std::cout << "Data of the run YYMMDD format (e.g. 240228)" << std::endl;
      std::cin >> date;
      std::cout << "Old/new electronic (old=0, new=1)" << std::endl;
      std::cin >> electronic;
      std::cout << "Add a comment (type n if no comment); e.g. \"am or pm\", bias..." << std::endl;
      std::cin >> comment;

      feature_value.push_back({"Date", date});
      feature_value.push_back({"Electronic (old=0, new=1)", electronic});
      feature_value.push_back({"A fast", fit_a_fast});
      feature_value.push_back({"A fast err", err_a_fast});
      feature_value.push_back({"Tau fast [mus]", fit_tau_fast});
      feature_value.push_back({"Tau fast err [mus]", err_tau_fast});
      feature_value.push_back({"A slow", fit_a_slow});
      feature_value.push_back({"A slow err", err_a_slow});
      feature_value.push_back({"Tau slow [mus]", fit_tau_slow});
      feature_value.push_back({"Tau slow err [mus]", err_tau_slow});
      feature_value.push_back({"I fast", fit_a_fast*fit_tau_fast/(fit_a_fast*fit_tau_fast+fit_a_slow*fit_tau_slow)});
      feature_value.push_back({"I slow ", fit_a_slow*fit_tau_slow/(fit_a_fast*fit_tau_fast+fit_a_slow*fit_tau_slow)});
      feature_value.push_back({"I slow pure", fit_a_slow*1.5/(fit_a_fast*fit_tau_fast+fit_a_slow*1.5)});

      print_vec_pair_csv(pde_result_file, feature_value, comment);
    }


    std::cout << "The best parameters are " << std::endl;
    std::cout << "a_fast   = " << par[0] << std::endl;
    std::cout << "tau_fast = " << par[1] << std::endl;
    std::cout << "a_slow   = " << par[2] << std::endl;
    std::cout << "tau_slow = " << par[3] << std::endl;
  
    std::cout <<"\nFast/Slow intensities" << std::endl;
    std::cout <<"I_fast = " << fit_a_fast*fit_tau_fast/(fit_a_fast*fit_tau_fast+fit_a_slow*fit_tau_slow) << std::endl;  
    std::cout <<"I_slow = " << fit_a_slow*fit_tau_slow/(fit_a_fast*fit_tau_fast+fit_a_slow*fit_tau_slow) << std::endl;  
    std::cout <<"I_slow pure = " << fit_a_slow*1.5/(fit_a_fast*fit_tau_fast+fit_a_slow*1.5) << std::endl;  
  }


  // ---- PLOTS --------------------------------------------------------
  avg_muon = avg_muon_original;
  vector_roll(avg_muon, best_fit_roll);
  TGraphErrors* g_muon = new TGraphErrors(memorydepth, &time[0], &avg_muon[0],
                                          &e_x[0], &e_y[0]);

  for(size_t i=0; i<memorydepth; i++) sin_muon[i]=xy[i];
  TGraph* g_sint = new TGraph(memorydepth, &time[0], &sin_muon[0]);

  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,800);
  c2->cd();
  g_muon->GetXaxis()->SetTitle("Time [#mus]");
  g_muon->GetYaxis()->SetTitle("Amplitude [a.u.]");

  g_sint->SetLineWidth(2);
  g_sint->SetLineColor(kRed);
  g_muon->Draw();
  g_sint->Draw("same");
  c2->Modified(); c2->Update();

  return;
}
