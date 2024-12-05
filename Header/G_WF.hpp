//
//  G_WF.hpp
//
//  Created by Federico Galizzi on 02/10/23
//

#ifndef G_WF_hpp
#define G_WF_hpp

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <vector>
// #include "G_Read.hpp"
#include "G_Utility.hpp"

#include <TH2.h>
#include <TComplex.h>
#include <TVirtualFFT.h>
#include <TCanvas.h>


using namespace std;


//*********************************************
template <typename T>
void avgWF (const vector<vector<T>>& y, vector<T>& avg_wf){
//*********************************************
  size_t len = y[0].size();
  
  if (avg_wf.size() < 5) avg_wf.resize(len, 0);
  
  for (auto& wf : y) {
    for (unsigned int j = 0; j < len; j++) avg_wf[j] += wf[j];
  }

  double w = 1./y.size();
  for (int i = 0; i < len; i++) avg_wf[i] *= w;
}

// Subtract baseline to all the len-long waveform in all_wf. Baseline
// computed in pre-trigger
//*********************************************
void SubBaseline(std::vector<std::vector<double>>& all_wf, int pre, bool invert){
//*********************************************
  if (pre<20) {
    std::cout << "Too few ticks to subtract the baseline \t \t " ;
    exit(0);
  }
   
  int WFs = all_wf.size();
  int len = all_wf[0].size();
  double baseline = 0.;
  
  for(int i=0; i<WFs; i++){
    baseline = 0.;
    for (int j = 0; j<pre; j++) baseline += all_wf[i][j];

    baseline /= (double) pre;
    if(invert==false) for (int j=0; j<len; j++) all_wf[i][j] -= baseline;
    else for (int j=0; j<len; j++) all_wf[i][j] = -all_wf[i][j]+baseline;
  }
}


// Subtract baseline using the most provable value
//********************************************
void SubBaseline2(std::vector<std::vector<double>>& all_wf, double rms, bool invert){
//*********************************************
  double t, norm;
  double low = -3*rms;
  double up  =  3*rms;
  double min_sampling = double(all_wf[0].size()*0.05);
  
  for (auto& wf : all_wf){
    t = Vector_MPV(wf);
    if(invert==false) for (auto& v : wf) v -= t;
    else for (auto& v : wf) v = -v+t;
  }

  for (auto& wf : all_wf){
    norm = 0.;
    t = 0;
    
    for (auto& v : wf){
      if(v > low && v < up){
        t += v; norm += 1;
      }
    }

    if(norm > min_sampling){
      t = t/norm;
      for (auto& v : wf) v -= t; 
    }
  }
}


// With the entire set of WFs (all_wf) it  build the calibration histogram integrating [I_low;I_up]
//*********************************************
TH1D* BuildRawChargeHisto(std::vector<std::vector<double>>& all_wf , std::vector<double>&int_wf, 
    int I_low, int I_up, int nbins){
//*********************************************
  int len = all_wf[0].size();
  int_wf.erase(int_wf.begin(), int_wf.end());
  double t = 0.;

  for(auto& wf : all_wf){
    t = 0;
    for(int i=I_low; i<=I_up; i++) t += wf[i];
    int_wf.push_back(t);
  }

  //Find the best x-axis range of the histogram
  // Duplicate to preserve the int_wf order
  std::vector int_dual = int_wf;
  std::sort(int_dual.begin(), int_dual.end());
  double x_low = int_dual[int(int_dual.size()*0.01)];
  double x_up  = int_dual[int(int_dual.size()*0.99)];

  TH1D* hI  = new TH1D("hI" ,"hI", nbins, x_low, x_up);
  
  for (auto& val : int_wf) hI->Fill(val);
 
  return hI;
}

// With the entire set of WFs (all_wf) it  build the calibration histogram integrating [I_low;I_up]
//*********************************************
TH1D* BuildRawChargeHisto(std::vector<std::vector<double>>& all_wf , std::vector<double>& int_wf,
    int I_low, int I_up, double hmin, double hmax, int nbins){
//*********************************************
  int len = all_wf[0].size();
  int_wf.erase(int_wf.begin(), int_wf.end());
  double t = 0.;
    
  for(auto& wf : all_wf){
    t = 0;
    for(int i=I_low; i<=I_up; i++) t += wf[i];
    int_wf.push_back(t);
  }

  TH1D* hI  = new TH1D("hI" ,"hI", nbins, hmin, hmax);
  
  for (auto& val : int_wf) hI->Fill(val);
  
  return hI;
 
}

// With the non-saturating WFs (ns_wf) it  build the F_prompt histogram integrating [I_low;I_up] and [I_low;I_prompt]
//*********************************************
TH1D* BuildFpromptHisto(std::vector<std::vector<double>>& ns_wf, std::vector<std::vector<double>>& mu_wf, 
    int I_low, int I_up, int I_pr, double f_pr){
//*********************************************
  double t, prompt_integral, integral;
  int len = ns_wf[0].size();
  std::vector<double> f_wf;
  
  for(auto& wf: ns_wf){
    prompt_integral= accumulate(wf.begin()+I_low, wf.begin()+I_pr, 0.);
    integral= accumulate(wf.begin()+I_low, wf.begin()+I_up, 0.);
    t = prompt_integral/integral;

    f_wf.push_back(t);
    if (t < f_pr) mu_wf.push_back(wf);
  }

  TH1D* hI  = new TH1D("hI" ,"hI", 1000, 0., 1.);
  for (auto& val : f_wf) hI->Fill(val);
  
  return hI;
}

/*
// With the non-saturating WFs (ns_wf) build the F_prompt histogram
// integrating [I_low;I_up] and [I_low;I_prompt]
// *********************************************
TH1D* AllFpromptHisto(std::vector<std::vector<double>>& ns_wf, int I_low, int I_up, int I_pr){
// *********************************************
  double t, prompt_integral, integral;
  int len = ns_wf[0].size();
  std::vector<double> f_wf;
  TH1D* hwf = new TH1D("hwf","hwf",len,0,len);
  
  for(auto& wf: ns_wf){
    prompt_integral= accumulate(wf.begin()+I_low, wf.begin()+I_pr, 0.);
    integral= accumulate(wf.begin()+I_low, wf.begin()+I_up, 0.);
    t = prompt_integral/integral;
    f_wf.push_back(t);
  }

  TH1D* h_p  = new TH1D("h_prompt" ,"h_prompt", 500, 0., 1.);
  for (auto& val : f_wf) h_p->Fill(val);
  
  return h_p;
}
*/

// Build the F_prompt vs Charge 2d-histogram integrating [I_low;I_up] and [I_low;I_prompt]
//*********************************************
TH2D* BuildChargeFpromptHisto(vector<vector<double>>& wfs,
                              int I_low, int I_up, int I_pr){
//*********************************************
  int len = wfs[0].size();
  size_t nwfs = wfs.size();
  std::vector<double> f_wf, i_wf;
  
  for(auto& wf: wfs){
    double prompt_integral= accumulate(wf.begin()+I_low, wf.begin()+I_pr, 0.);
    double integral= accumulate(wf.begin()+I_low, wf.begin()+I_up, 0.);
    f_wf.push_back(prompt_integral/integral);
    i_wf.push_back(integral);
  }
  
  double ymin = *min_element(std::begin(i_wf), std::end(i_wf));
  double ymax = *max_element(std::begin(i_wf), std::end(i_wf));
  
  TH2D* h2_charge_fprompt = new TH2D("h2_charge_fprompt",
                                     Form("%s;%s;%s", "fprompt_vs_charge", "Charge [ADC #times ticks]", "F prompt"),
                                     500, ymin, ymax, 500, 0., 1.);
  for (size_t i=0; i<nwfs; i++) h2_charge_fprompt->Fill(i_wf[i], f_wf[i]);
  
  return h2_charge_fprompt;
}

// With the non-saturating WFs (ns_wf) it  build the F_prompt histogram integrating [I_low;I_up] and [I_low;I_prompt]
//*********************************************
TH2D* BuildChargeFpromptHisto(vector<vector<double>>& ns_wf, vector<vector<double>>& mu_wf, 
                              int I_low, int I_up, int I_pr, double f_pr){
//*********************************************
  double tx, ty, prompt_integral, integral;
  int len = ns_wf[0].size();
  size_t wfs = ns_wf.size();
  std::vector<double> f_wf, i_wf;
  
  for(auto& wf: ns_wf){
    prompt_integral= accumulate(wf.begin()+I_low, wf.begin()+I_pr, 0.);
    integral= accumulate(wf.begin()+I_low, wf.begin()+I_up, 0.);
    f_wf.push_back(prompt_integral/integral);
    i_wf.push_back(integral);
    
    if (prompt_integral/integral < f_pr) mu_wf.push_back(wf);
  }
  
  tx = *min_element(std::begin(i_wf), std::end(i_wf));
  ty = *max_element(std::begin(i_wf), std::end(i_wf));
  
  TH2D* hI  = new TH2D("hI" ,"hI", 500, tx, ty, 500, 0., 1.);
  for (size_t i=0; i<wfs; i++) hI->Fill(i_wf[i], f_wf[i]);
  
  return hI;
}
// Select WFs on their integral basis, store them in sel_wf and compute their average spe_wf
//*********************************************
void Avg_Sel_WF (std::vector<std::vector<double>>& all_wf,
    std::vector<std::vector<double>>& sel_wf, 
    std::vector<double>& spe_wf, const std::vector<double>& int_wf,
    double I_low, double I_up){
//*********************************************
  int nspe_wf=0;
  int len = all_wf[0].size();
  
  for (int i = 0; i < int_wf.size(); i++) {
    if (int_wf[i] > I_low && int_wf[i] < I_up) {
      nspe_wf += 1;
      sel_wf.push_back(all_wf[i]);
    }
  }
  
  avgWF(sel_wf, spe_wf);
  std::cout << "N_sel " << nspe_wf << std::endl;

}

// Same, without storing the selected WFs
//*********************************************
void Avg_Sel_WF (std::vector<std::vector<double>>& all_wf, std::vector<double>& spe_wf, 
    const std::vector<double>& int_wf, double I_low, double I_up){
//*********************************************
  int nspe_wf=0;
  size_t len = all_wf[0].size();
  spe_wf.erase(spe_wf.begin(), spe_wf.end());
  spe_wf.resize(len, 0.);
  
  for (int i = 0; i < int_wf.size(); i++) {
    if (int_wf[i] > I_low && int_wf[i] < I_up) {
      nspe_wf += 1;
      for (size_t j=0; j<len; j++) spe_wf[j] += all_wf[i][j];
    }
  }
  
  double norm = double(1./nspe_wf);
  for (size_t j=0; j<len; j++) spe_wf[j] *= norm;
  cout << "\nAvg_sel_WF() WFs candidates " << nspe_wf << endl;
}

// Build a TH2D with all the PSD of the waveform and compute the average
//*********************************************
TGraph* build_avg_spectral_density(int nsample, double t1, double t0,
    std::vector<std::vector<double>>& wf) {
//*********************************************
  double dt = (t1-t0)/nsample;
  const int nsample_ = nsample;
  double nwindow = wf.size();
  double    xn[nsample_];
  TComplex  xN[nsample_];
  double    xN_re[nsample_];
  double    xN_im[nsample_];
  double scale = 1./nwindow;
  double c_scale = 1./nsample;
  

  int nsample_fft = 0.5*nsample;
  TGraph* g_avg_spectral_density = new TGraph(nsample_fft);
  for (int j=0; j<nsample_fft; j++)
    g_avg_spectral_density->SetPoint(j, j/t1, 0.);
  double ymin = 1e-6;
  double ymax = 9e+1;
  int    nbinsy = 100;
  
  std::vector<double> ybin_exponent = linspace(TMath::Log10(ymin), TMath::Log10(ymax), nbinsy);
  std::vector<double> ybins(100, 0);
  for (int iy=0; iy<100; iy++) ybins[iy] = TMath::Power(10., ybin_exponent[iy]);

  TH2D* h2_spectral_density = new TH2D("h2_spectral_density",
      Form("%s;%s;%s", "Noise spectral density", "Frequency [MHz]", "Amplitude [A.U.]"),
      nsample_fft,
      0.,
      (double) 0.5*nsample/t1,
      nbinsy-1, &ybins.at(0));
  h2_spectral_density->GetXaxis()->CenterTitle();
  h2_spectral_density->GetYaxis()->CenterTitle();
  
  //  FFT
  TVirtualFFT* fft = TVirtualFFT::FFT(1, &nsample, "M R2C");

  for (int iw=0; iw<nwindow; iw++) {
    for (int ip=0; ip<nsample_; ip++) xn[ip] = (double) wf[iw][ip];
    
    fft->SetPoints(xn);
    fft->Transform();
    fft->GetPointsComplex(xN_re, xN_im);

    for (int j=0; j<nsample_fft; j++) {
      xN[j] = TComplex(xN_re[j], xN_im[j])*c_scale;
      h2_spectral_density->Fill(j/t1, xN[j].Rho2());
      g_avg_spectral_density->GetY()[j] += (xN[j].Rho2()*scale);
    }
  }

  TCanvas* cNoise = new TCanvas("cNoise", "Noise spectral density", 100, 100, 800, 600);
  cNoise->SetLogy(1);
  cNoise->SetLogx(1);
  cNoise->SetTicks(1, 1);
  cNoise->SetGrid(1, 1);
  h2_spectral_density->Draw("colz");
  g_avg_spectral_density->SetLineColor(kGray+2);
  g_avg_spectral_density->SetLineWidth(2);
  g_avg_spectral_density->Draw("l");
  
  return g_avg_spectral_density;
}

// Build a TH2D with all the PSD of the waveform and compute the average
//*********************************************
TGraph* build_avg_spectral_density(int nsample, double t1, double t0,
    std::vector<std::vector<double>>& wf, double res) {
//*********************************************
  double dt = (t1-t0)/nsample;
  const int nsample_ = nsample;
  double nwindow = wf.size();
  double    xn[nsample_];
  TComplex  xN[nsample_];
  double    xN_re[nsample_];
  double    xN_im[nsample_];
  double scale = 1./nwindow;
  double c_scale = 1./nsample;
  double t;
  

  int nsample_fft = 0.5*nsample+1;
  TGraph* g_avg_spectral_density = new TGraph(nsample_fft);
  for (int j=0; j<nsample_fft; j++)
    g_avg_spectral_density->SetPoint(j, j/t1, 0.);
  double ymin = -100;
  double ymax = -20;
  int    nbinsy = 100;

  TH2D* h2_spectral_density = new TH2D("h2_spectral_density",
      Form("%s;%s;%s", "Noise spectral density", "Frequency [MHz]", "Power Spectral Density [db]"),
      nsample_fft,
      0.,
      (double) 0.5*nsample/t1,
      nbinsy-1, ymin, ymax);
  
  //  FFT
  TVirtualFFT* fft = TVirtualFFT::FFT(1, &nsample, "M R2C");

  for (int iw=0; iw<nwindow; iw++) {
    for (int ip=0; ip<nsample_; ip++) xn[ip] = (double) wf[iw][ip];
    
    fft->SetPoints(xn);
    fft->Transform();
    fft->GetPointsComplex(xN_re, xN_im);

    for (int j=0; j<nsample_fft; j++) {
      xN[j] = TComplex(xN_re[j], xN_im[j]);
      t = 10*TMath::Log10(xN[j].Rho2()*c_scale/(pow(2,res*2)));
      //t = 20*TMath::Log10(xN[j].Rho2()*c_scale/(pow(2,res)));
      h2_spectral_density->Fill(j/t1, t);
      g_avg_spectral_density->GetY()[j] += (t*scale);
    }
  }

  TCanvas* cNoise = new TCanvas("cNoise", "Power Spectral Density", 100, 100, 800, 600);
  cNoise->SetLogz(1);
  cNoise->SetLogx(1);
  cNoise->SetTicks(1, 1);
  cNoise->SetGrid(1, 1);
  h2_spectral_density->Draw("colz");
  g_avg_spectral_density->SetLineColor(kGray+2);
  g_avg_spectral_density->SetLineWidth(2);
  g_avg_spectral_density->Draw("l");
  
  return g_avg_spectral_density;
}

// As build_avg_spectral_density but without TH2D stuff
//*********************************************
TGraph* build_ch_fft(int nsample, double t1, double t0,
    std::vector<std::vector<double>>& wf, double res) {
//*********************************************
  double dt = (t1-t0)/nsample;
  const int nsample_ = nsample;
  double nwindow = wf.size();
  double    xn[nsample_];
  TComplex  xN[nsample_];
  double    xN_re[nsample_];
  double    xN_im[nsample_];
  double scale = 1./nwindow;
  double c_scale = 1./nsample;
  double t;
  

  int nsample_fft = 0.5*nsample;
  TGraph* g_avg_spectral_density = new TGraph(nsample_fft);
  for (int j=0; j<nsample_fft; j++)
    g_avg_spectral_density->SetPoint(j, j/t1, 0.);
  double ymin = -100;
  double ymax = -20;
  int    nbinsy = 100;

    //  FFT
  TVirtualFFT* fft = TVirtualFFT::FFT(1, &nsample, "M R2C");

  for (int iw=0; iw<nwindow; iw++) {
    for (int ip=0; ip<nsample_; ip++) xn[ip] = (double) wf[iw][ip];
    
    fft->SetPoints(xn);
    fft->Transform();
    fft->GetPointsComplex(xN_re, xN_im);

    for (int j=0; j<nsample_fft; j++) {
      xN[j] = TComplex(xN_re[j], xN_im[j]);
      t = 10*TMath::Log10(xN[j].Rho2()*c_scale/(pow(2,res*2)));
      //t = 20*TMath::Log10(xN[j].Rho2()*c_scale/(pow(2,res)));
      g_avg_spectral_density->GetY()[j] += (t*scale);
    }
  }

    return g_avg_spectral_density;
}
//*********************************************
template <typename T>
void SubVec_to_WFs(vector<vector<T>>& y, vector<T>& sub){
//*********************************************
  size_t len = sub.size();
  for(auto& wf : y)
    for (size_t j = 0; j < len; j++) wf[j] -= sub[j];
}


// Pick up the first saturating WF 
//*********************************************
template <typename T>
void Sat_WF(vector<vector<T>>& y, vector<vector<T>>& y2, T sat_up){
//*********************************************
  T max_el;
  size_t len = y[0].size();

  for (auto wf : y) {
    max_el = *max_element( wf.begin(), wf.end());
    
    if (max_el>sat_up) {
      y2.push_back(wf);
      return;
    }
  }
    
  std::cout << "\n \n !!! \n No saturating WF \n Threshold was set to " << sat_up << "\n \n" ;
}

// Select the waveforms where the prepulse-trigger ticks are within the
// [-bsl;+bsl], the pulse amplitude in [range_low;range_up] and without lower
// saturation (sat_low)
//*********************************************
template <typename T>
void SelPDE_WF(vector<vector<T>>& y, vector<vector<T>>& y2, int pre, int int_prompt,
               T sat_low, T range_low, T range_up, T bsl, T rms){
//*********************************************
  T max_el, min_el, t;
  size_t len = y[0].size();
  size_t wfs = y.size();
  int bsl_counter = 0;
  int amp_counter = 0;
  int sel_counter = 0;

  vector<double> first_avg(len, 0.0);
  vector<bool> preliminary_selection(wfs, false);
  vector<bool> final_selection(wfs, false);

  y2.reserve(wfs);
    
  for (size_t i=0; i<wfs; i++) {
    max_el = *max_element( y[i].begin(), y[i].begin()+pre);
    min_el = *min_element( y[i].begin(), y[i].begin()+pre);
   
    // Select wfs with no pulses in the pre-trigger
    if (max_el<bsl && min_el > -bsl) {
      bsl_counter++; 
      max_el = *max_element( y[i].begin()+pre, y[i].begin()+int_prompt);
      min_el = *min_element( y[i].begin()+pre, y[i].begin()+int_prompt);
      // Select wfs with amplitude in given range and with no lower saturation
      if (max_el<range_up && max_el>range_low && min_el > sat_low){
        amp_counter++;
        t = max_el;
        max_el = *max_element( y[i].begin()+int_prompt, y[i].end());
        min_el = *min_element( y[i].begin()+int_prompt, y[i].end());
        //
        if (max_el < t && min_el > sat_low){
          sel_counter++;
          for(size_t j=0; j<len; j++) first_avg[j] += y[i][j];
          preliminary_selection[i]=true;
        } 
      }     
    }
  }  
 
  double norm = 1./ *max_element(first_avg.begin(), first_avg.end());
  for (auto& e : first_avg) e*=norm;

  for (size_t i=0; i<wfs; i++) {
    if(preliminary_selection[i]==true){
      vector wf = y[i];
      norm = 1./ *max_element(wf.begin(), wf.end());
      double norm_bsl = rms*norm*5.;
      for (auto& e : wf) e *= norm;

      bool select_shape = true;
      for (size_t j=0; j<len; j++){
        if (wf[j] > first_avg[j]+norm_bsl || wf[j] < first_avg[j]-norm_bsl){
          select_shape = false;
          j = len+1;
        }
      }
      final_selection[i]=select_shape;
    }

    if(final_selection[i]==true) y2.push_back(y[i]);
  }
  
  cout << "WFs selected on the baseline "  << bsl_counter << "/" << wfs << endl;
  cout << "WFs selected on the amplitude " << amp_counter << "/" << wfs << endl;
  cout << "WFs preliminary selected " << sel_counter << "/" << wfs << endl;
  cout << "WFs selected " << y2.size() << "/" << wfs << endl;

  return;
}

// Return the waveform where the prepulse-trigger ticks are within the
// sat_low-sat_up range
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

// Select waveforms used for calibration on their integral (from spe_low to inf)
// and shape: it must be equal to the one of an average wf up to 5*bsl fluctuations
//*********************************************
template <typename T>
void SelTemplate_WF(const vector<vector<T>>& calib_wfs, vector<vector<T>>& template_wfs,
                    const vector<T>& int_wfs, const T spe_low, const T bsl,
                    const int int_low, const int int_up){
//*********************************************
  T max_el, min_el, norm, norm_bsl;
  size_t len = calib_wfs[0].size();
  size_t wfs = calib_wfs.size();

  vector<T> calib_avg_wf;
  vector<bool> selection_flag;
  
  avgWF(calib_wfs, calib_avg_wf);
  max_el = *max_element( calib_avg_wf.begin()+int_low, calib_avg_wf.begin()+int_up);
  norm = 1/max_el;
  for (auto& e : calib_avg_wf) e *= norm;

  for (size_t i=0; i<wfs; i++){
    if(int_wfs[i]>spe_low){
      auto wf = calib_wfs[i];
      max_el = *max_element(wf.begin()+int_low, wf.begin()+int_up);
      norm = 1/max_el;
      norm_bsl = bsl*norm*5.;
      for (auto& e : wf) e *= norm;
      bool flag = true;
      for (size_t j=0; j<len; j++){
        if (wf[j] > calib_avg_wf[j]+norm_bsl ||
            wf[j] < calib_avg_wf[j]-norm_bsl){
          flag = false;
          continue;
        }
      }
      selection_flag.push_back(flag);
    }
    else{
      selection_flag.push_back(false);
    }
  }

  if(selection_flag.size()!=wfs){
    std::cout << "Problem in selecting template wfs" << std::endl;
    return;
  }

  for (size_t i=0; i<wfs; i++){
    if (selection_flag[i]==true) template_wfs.push_back(calib_wfs[i]);
  }

  return;
}

// Display num waveforms belonging to a single vector
//*********************************************
template <typename T>
void DisplayWFs (const vector<vector<T>>& y, T tt, int num){
//*********************************************
  size_t len = y[0].size();

  for (int i = 0; i<num; i++) {
    TH1D *waa = new TH1D("Waveform", "Waveform",len,0,len*tt);
    for(size_t bin = 0; bin < len; bin++){
      waa->SetBinContent(bin+1, y[i][bin]);}
    waa->Draw();gPad->Update();gPad->WaitPrimitive();
    delete waa;
  }

}

// Display num waveforms belonging to two vectors
//*********************************************
template <typename T>
void DisplayWFs (const vector<vector<T>>& y, const vector<vector<T>>& y2, T tt, int num){
//*********************************************
  int len = y[0].size();

  for (int i = 0; i<num; i++) {
    TH1D *h1 = new TH1D("h1", "h1",len,0,len*tt);
    TH1D *h2 = new TH1D("h2", "h2",len,0,len*tt);
    h1->SetLineColor(2);
    h2->SetLineColor(4);
    for(int bin = 0; bin < len; bin++){
      h1->SetBinContent(bin+1, y[i][bin]);
      h2->SetBinContent(bin+1, y2[i][bin]);
    }

    h1->Draw();h2->Draw("same");gPad->Update();gPad->WaitPrimitive();
    delete h1; delete h2;
  }
}

//
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

// Filter all the wfs according to the G filter
//*********************************************
void FilterAllWF(const vector<vector<double>>& all_wf, vector<vector<double>>& filt_wf, TComplex G[]){
//*********************************************
  int len = all_wf[0].size();
  filt_wf.resize(all_wf.size(), vector<double>(len));
  double xv[len];
  double* xy;
  double c_scale = 1./len;
  TComplex xV[len]; double xV_re[len]; double xV_im[len];
  TComplex xY[len]; double xY_re[len]; double xY_im[len];
  TVirtualFFT* fft = TVirtualFFT::FFT(1, &len, "M R2C");
  
  for (int i=0; i<all_wf.size(); i++) {
    for (int j=0; j<len; j++) xv[j] = all_wf[i][j];
    
    fft = TVirtualFFT::FFT(1, &len, "M R2C");
    fft->SetPoints(xv);
    fft->Transform();
    fft->GetPointsComplex(xV_re, xV_im);
    
    for (int j=0; j<len*0.5+1; j++) {
      xV[j] = TComplex(xV_re[j], xV_im[j]);
      xY[j] = G[j]*xV[j]; 
      xY_re[j] = xY[j].Re(); xY_im[j] = xY[j].Im();
    }
    
    fft = TVirtualFFT::FFT(1, &len, "M C2R");
    fft->SetPointsComplex(xY_re, xY_im);
    fft->Transform();
    xy = fft->GetPointsReal();
    
    for (int j=0; j<len; j++) filt_wf[i][j] = xy[j]*c_scale;
     
  }
}


 
//*********************************************
template <typename T>
void MovingAverageWF (std::vector<T> in, std::vector<T>& out, int w){
//*********************************************
  T sum = 0.;
 
  if (out.size()>0) out.erase(out.begin(), out.end());
  if (w <= 0) std::cout << "Change window \n" ;  
  
  for (int i=0; i<in.size(); i++) {
    sum += in[i];

    if (i >= w){
      sum -= in[i-w];
      out.push_back(sum/w);
    } else out.push_back( sum/((T)i+1.) );
 
  } 
}

//*********************************************
template <typename T>
void MovingAverageWF (vector<vector<T>> in, vector<vector<T>>& out, int w){
//*********************************************
  T sum = 0.;
  size_t wfs = in.size();
  size_t len = in[0].size();

  out.resize(wfs, vector<T>(len));
 
  if (w <= 0) std::cout << "Change window \n" ;  
 
  for(size_t i=0; i<wfs; i++){
    sum = 0;

    for (int j=0; j<len; j++) {
      sum += in[i][j];

      if (j >= w){
        sum -= in[i][j-w];
        out[i][j] = sum/w;
      } else out[i][j] = sum/((T)j+1.);
    }
  }   
}

//*********************************************
template <typename T>
void AllignWFs (vector<vector<T>> wfs){
//*********************************************
  double thr = *max_element(wfs[0].begin(), wfs[0].end)*0.5; 
  auto iter_start = std::find_if(wfs[0].begin(), wfs[0].end(), [thr](double value) {
    return value >= thr; });
  auto iter_ref = std::find_if(wfs[0].begin(), wfs[0].end(), [thr](double value) {
    return value >= thr; });

  for(auto wf : wfs){
    thr = *max_element(wf.begin(), wf.end())*0.5;
    iter_start = std::find_if(wf.begin(), wf.end(), [thr](double value) {
      return value >= thr; });

  }

}

//*********************************************
vector<double> TriggerTime(vector<double>& waveform){
//*********************************************
  vector<double> trgs;
  for(size_t i=0; i<waveform.size(); i++) if (waveform[i] > 0.5){
    trgs.push_back(i);
    while(waveform[i]>0.5) i++;
  }
  return trgs;
}


#endif /* G_WF_hpp */
