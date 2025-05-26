#include "../classe.hpp"

// ****************************************************************
// Macro to perform the convolution of the spe response with the
// liquid argon scintillation time profile to fit the muon signal.
// How to use:
// 1- set nofit=1 to see how good are the fit parameters you set
//    Parameters to set are: amp, f_fast, tau_fast, tau_slow, roll,
//    yerr, fit_l and fit_u.
//    The roll parameter is used to adjust the convolution offset
//    The fit_l and fit_u parameters are used to set the fit range
//    Then amp=amplitude, f_fast=fast fraction, tau_fast, tau_slow
//    are the the tau of fast and slow components.
//    The yerr parameter is used to set the error on the muon wf.
// 2- set nofit=0 to perform the fit
// ****************************************************************

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////
string pde_result_folder = "/eos/home-f/fegalizz/PDE_MiB/PDE_Results/Conv_results/";
string pde_result_file = pde_result_folder+"res.csv";
///////////////////////////////////////////////////////////////////

#include "Fit/FitResult.h"
#include "Rtypes.h"
#include "TComplex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TVirtualFFT.h"
#include "TLegend.h"

#include <TF1.h>
#include <cstddef>
#include <iostream>
#include <vector>

#include <Math/Functor.h>
#include "Fit/Fitter.h"

//*********************************************
void conv_templfft_arraytd(double* conv, TComplex* templ_fft, double* array_td, size_t len){
//*********************************************
  TComplex array_fft[len];
  TComplex xY[len];
  int len_ = int(len);
  double norm = 1./double(len);
  
  compute_r2c_fft(&array_fft[0], &array_td[0], len_);
  freq_domain_convolution(&templ_fft[0], &array_fft[0], xY, len);
  compute_c2r_fft(&conv[0], xY, len);

  for (int j=0; j<len_; j++) conv[j] *= norm;

<<<<<<< HEAD
  return;
=======
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
>>>>>>> deff1318b99aedb8334620d3240944927f952bce
}


//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Convolution(){

  size_t   nsample = memorydepth;
  double norm = 1./double(memorydepth);
  TComplex templ_fft[nsample], dexp_fft[nsample];
  double   templ_td[nsample], dexp_td[nsample], convolution[nsample];
  
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
  compute_r2c_fft(&templ_fft[0], &templ_td[0], memorydepth);

  double params[4] = {amp, f_fast, tau_fast, tau_slow}; 
  double* par = &params[0];
  double_expo_timedomain(&dexp_td[0], par, memorydepth, tick_len);
  conv_templfft_arraytd(&convolution[0], &templ_fft[0], &dexp_td[0], nsample);
  TComplex xY[nsample]; double xY_re[nsample]; double xY_im[nsample];


  // ---- FIT -------------------------------------------------------- 
  int best_fit_roll = 0;
  double min_chi2 = 1.e10;
  double NDf = 0.;
  int int_fit_l = int(fit_l/tick_len);
  int int_fit_u = min(int(fit_u/tick_len), memorydepth);
  bool fit_check;

  if(!nofit){
    double best_params[4];
    double fit_amp, err_amp;  
    double fit_f_fast, err_f_fast;  
    double fit_tau_fast, err_tau_fast;
    double fit_tau_slow, err_tau_slow;
 
    // Scan different t0s
    for (size_t idx_fit=0; idx_fit<30; idx_fit++){
      avg_muon = avg_muon_original;
      int fit_roll = idx_fit-15;
      if (idx_fit!=0) vector_roll(avg_muon, fit_roll); // to start from your best guess
      TGraphErrors* g_muon = new TGraphErrors(memorydepth, &time[0], &avg_muon[0], 
                                              &e_x[0], &e_y[0]);

      auto chi2Function = [&](const double* par){
        double chi2 = 0.;
        double* y = g_muon->GetY();

        compute_c2r_fft(&convolution[0], xY, memorydepth);
        double_expo_timedomain(&dexp_td[0], par, memorydepth, tick_len);
        conv_templfft_arraytd(&convolution[0], &templ_fft[0], &dexp_td[0], nsample);

        for (int i=int_fit_l; i<int_fit_u; i++){
          double d = y[i]-convolution[i];
          chi2 += d*d/(yerr*yerr);
        }

        return chi2;
      };

      ROOT::Math::Functor fnc(chi2Function, 4);
      ROOT::Fit::Fitter fitter;
      fitter.SetFCN(fnc, par);
      fitter.Config().ParSettings(0).SetName("Amplitude");
      fitter.Config().ParSettings(1).SetName("Fast fraction");
      fitter.Config().ParSettings(2).SetName("#tau_{fast}");
      fitter.Config().ParSettings(3).SetName("#tau_{slow}");

      std::cout << "---- Fit " << idx_fit << " ---------" << std::endl;
      fit_check = fitter.FitFCN();
      fitter.SetNumberOfFitPoints(static_cast<size_t>(int_fit_u-int_fit_l));
      auto result = fitter.Result();
      result.Print(std::cout);
      std::cout << "\n" << std::endl;

      if (result.MinFcnValue()<min_chi2 && fit_check==true){
        std::cout << result.GetParams()[0] << " " << result.GetParams()[0] << std::endl;
        min_chi2 = result.MinFcnValue();
        NDf = result.Ndf();
        best_fit_roll = fit_roll;
        for (size_t i=0; i<4; i++) best_params[i] = result.GetParams()[i];
        fit_amp      = result.GetParams()[0]; err_amp      = result.GetErrors()[0];
        fit_f_fast   = result.GetParams()[1]; err_f_fast   = result.GetErrors()[1];
        fit_tau_fast = result.GetParams()[2]; err_tau_fast = result.GetErrors()[2];
        fit_tau_slow = result.GetParams()[3]; err_tau_slow = result.GetErrors()[3];
      }
      
      par  = &best_params[0];

      double_expo_timedomain(&dexp_td[0], par, memorydepth, tick_len);
      conv_templfft_arraytd(&convolution[0], &templ_fft[0], &dexp_td[0], nsample);
      std::cout << "\n" << std::endl;
    }
 

    // Update class members with the best fit parameters
    amp      = fit_amp;   
    f_fast   = fit_f_fast;  
    tau_fast = fit_tau_fast;
    tau_slow = fit_tau_slow;
    par[0] = fit_amp;
   

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
      feature_value.push_back({"Fit low [mus]", fit_l});
      feature_value.push_back({"Fit up [mus]", fit_u});
      feature_value.push_back({"Amplitude", fit_amp});
      feature_value.push_back({"Amplitude err", err_amp});
      feature_value.push_back({"Fast fraction", fit_f_fast});
      feature_value.push_back({"Fast fraction err", err_f_fast});
      feature_value.push_back({"Tau fast [mus]", fit_tau_fast});
      feature_value.push_back({"Tau fast err [mus]", err_tau_fast});
      feature_value.push_back({"Tau slow [mus]", fit_tau_slow});
      feature_value.push_back({"Tau slow err [mus]", err_tau_slow});
      feature_value.push_back({"Y error", yerr});
      feature_value.push_back({"Min chi2", min_chi2});
      feature_value.push_back({"NDf", NDf});

      print_vec_pair_csv(pde_result_file, feature_value, comment);
      update_thisfile(pde_result_folder+Form("Convolution_ana_parameters_%i_el_%i", int(date), int(electronic)));
    }


    cout << "\n---------------------------------------------------  \n" << endl;
    cout << "The best fit gives: " << std::endl;
    cout << "Min chi2  = " << min_chi2 << std::endl;
    cout << "NDf       = " << NDf << std::endl;
    cout << "Chi2/NDf  = " << min_chi2/NDf << std::endl;
    cout << "Ampl      = " << par[0] << std::endl;
    cout << "Fast frac = " << par[1] << std::endl;
    cout << "tau_fast  = " << par[2] << std::endl;
    cout << "tau_slow  = " << par[3] << std::endl;
    cout << "\n---------------------------------------------------  \n" << endl;
  }


  // ---- PLOTS --------------------------------------------------------
  avg_muon = avg_muon_original;
  vector_roll(avg_muon, best_fit_roll);
  TGraphErrors* g_muon = new TGraphErrors(memorydepth, &time[0], &avg_muon[0],
                                          &e_x[0], &e_y[0]);
<<<<<<< HEAD
  
  for(int i=0; i<memorydepth; i++) sin_muon[i]=convolution[i];
=======

  for(size_t i=0; i<memorydepth; i++) sin_muon[i]=xy[i];
>>>>>>> deff1318b99aedb8334620d3240944927f952bce
  TGraph* g_sint = new TGraph(memorydepth, &time[0], &sin_muon[0]);
  g_sint->SetLineWidth(2);
  g_sint->SetLineColor(kRed);

  auto legend = new TLegend(0.6,0.7,0.9,0.9);
  if(!nofit){
    double last_fit_params[4] = {par[0], par[1], par[2], par[3]};
    double* pp = &last_fit_params[0];

    auto chi2Function = [&](const double* pp){
      double chi2 = 0.;
      double* y = g_muon->GetY();

      compute_c2r_fft(&convolution[0], xY, memorydepth);
      double_expo_timedomain(&dexp_td[0], pp, memorydepth, tick_len);
      conv_templfft_arraytd(&convolution[0], &templ_fft[0], &dexp_td[0], nsample);

      for (int i=int_fit_l; i<int_fit_u; i++){
        double d = y[i]-convolution[i];
        chi2 += d*d/(yerr*yerr);
      }

      return chi2;
    };

    ROOT::Math::Functor fnc(chi2Function, 4);
    ROOT::Fit::Fitter fitter;
    fitter.SetFCN(fnc, pp);
    fitter.Config().ParSettings(0).SetName("Amplitude");
    fitter.Config().ParSettings(1).SetName("Fast fraction");
    fitter.Config().ParSettings(2).SetName("#tau_{fast}");
    fitter.Config().ParSettings(3).SetName("#tau_{slow}");

    std::cout << "---- Final Fit ------------------------" << std::endl;
    fit_check = fitter.FitFCN();
    fitter.SetNumberOfFitPoints(static_cast<size_t>(int_fit_u-int_fit_l));
    auto result = fitter.Result();
    result.Print(std::cout);

    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    std::string chi2_entry = Form("#chi2/n_{dof} = %.2f", result.MinFcnValue()/result.Ndf());
    std::string amp_entry = Form("A = %.2f #pm %.2f", result.GetParams()[0], result.GetErrors()[0]);
    std::string ffast_entry = Form("f_{fast} = %.3f #pm %.3f", result.GetParams()[1], result.GetErrors()[1]);
    std::string taufast_entry = Form("#tau_{fast} = %.3f #pm %.3f", result.GetParams()[2], result.GetErrors()[2]);
    std::string tauslow_entry = Form("#tau_{slow} = %.3f #pm %.3f", result.GetParams()[3], result.GetErrors()[3]);
    legend->AddEntry(g_muon, "Average waveform", "p");
    legend->AddEntry(g_sint, "Fit", "l");
    legend->AddEntry((TObject*)0, chi2_entry.c_str(), "");
    legend->AddEntry((TObject*)0, amp_entry.c_str(), "");
    legend->AddEntry((TObject*)0, ffast_entry.c_str(), "");
    legend->AddEntry((TObject*)0, taufast_entry.c_str(), "");
    legend->AddEntry((TObject*)0, tauslow_entry.c_str(), "");
    
  }


  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,800);
  c2->cd();
  g_muon->GetXaxis()->SetTitle("Time [#mus]");
  g_muon->GetYaxis()->SetTitle("Amplitude [ADC]");

  g_muon->Draw();
  if (!nofit) legend->Draw("same");
  g_sint->Draw("same");
  c2->Modified(); c2->Update();

  return;
}
