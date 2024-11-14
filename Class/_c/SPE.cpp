#include "../classe.hpp"
// Define apply_filter, spe_low, spe_up before running this macro!
// You can do it with LED_Analysis.cpp or Filt_Analysis.cpp and updating these
// variables with LoadFitParameters(fgaus)

void cla::SPE() {
  double y0, y1, r_time, f_time, undershoot;
  vector<double> x, avg, int_wf, t_templ;
  vector<vector<double>> sel_wf, spe_wf, avg_wf, filt_wf;
  size_t nsample = memorydepth;
  TComplex G[nsample];
    
  for (int i = 0; i < memorydepth; i++) x.push_back( (double) i*tick_len);
  
  // Read and subtract the baseline
  read();

  SelCalib_WF(wfs, sel_wf, prepulse_ticks, sat_low, sat_up, bsl);
  TH1D* hI = nullptr;
 
  //Build the charge histo according to moving average and filter option
  if(mov_win == true) {
    MovingAverageWF(sel_wf, sel_wf, win);
    hI = BuildRawChargeHisto(sel_wf, int_wf, int_low+win, int_up+win, nbins);
  }
  else{
    if(apply_filter == false) hI = BuildRawChargeHisto(sel_wf, int_wf, int_low, int_up, nbins);
    if(apply_filter == true){
      CompleteWF_Binary(templ_f, t_templ, memorydepth); // t_templ = time domain template
      Build_Matched_Filter(&G[0], t_templ);
      FilterAllWF(sel_wf, filt_wf, G);
      hI = BuildRawChargeHisto(filt_wf, int_wf, 0, 1, nbins);
    }
  }

  //Using sel_wf to preserve the raw information -> spe_wf are the single pe 
  //cadidates, avg their average 
  Avg_Sel_WF (sel_wf, spe_wf, avg, int_wf, spe_low, spe_up);
  if (spe_wf.size() == 0){
    class_skip = 1;
    return;
  }
  avg_wf.push_back(avg);
  RiseFallTimeUndershoot(avg, tick_len, r_time, f_time, undershoot);
  spe_ampl  = *std::max_element(avg.begin(), avg.end());
  spe_under = *std::min_element(avg.begin(), avg.end());


  if(print==true){
    VecDouble_in_Binary("Template.dat", avg);
    print = false;
  }
  //Draw the spe FFT
  TGraph* gAvg= build_avg_spectral_density(memorydepth, tick_len*memorydepth, tick_len, avg_wf, res);
  gAvg->SetLineColor(kGray+2); gAvg->SetLineWidth(2);
  gAvg->Draw("l");
  
  TCanvas *c_Spe_FFT = new TCanvas("c_Spe_FFT","c_Spe_FFT",20,20,1000,900);
  c_Spe_FFT->SetTicks(1, 1);      c_Spe_FFT->SetGrid(1, 1);
  c_Spe_FFT->cd();
  gAvg->Draw();
  c_Spe_FFT->Modified(); c_Spe_FFT->Update();

  //Draw the SPE wf
  TGraph *g1 = new TGraph(avg.size(), &x[0], &avg[0]);
  g1->GetXaxis()->SetTitle("Time [#mus]");
  g1->GetYaxis()->SetTitle("Amplitude [ADC]");
  g1->SetLineColor(2); g1->SetLineWidth(2);
  g1->SetTitle("Average spe-waveform");
  
  TCanvas *c_Spe = new TCanvas("c_Spe","c_Spe",20,20,1000,900);
  c_Spe->cd();
  g1->Draw("AL");
  c_Spe->Modified(); c_Spe->Update();
  
  //TERMINAL OUTPUTS + charghe histo for spe waveforms
  double overall_integral = 0;
  for(auto e : avg ) overall_integral += e;
  double spe_avg_gain = 0;
  for(int i=int_low; i<=int_up; i++) spe_avg_gain += avg[i];
  

  std::cout << "\n\nSPE integral [0;memorydepth] " << overall_integral << std::endl;
  std::cout << "SPE integral [int_low;int_up] " << spe_avg_gain << std::endl;
  

  std::cout << "\n\nColdbox table Ampl - under - r_time - f_time - under" << std::endl;
  std::cout << spe_ampl << "\t" << undershoot*spe_ampl/100. << "\t" << r_time << "\t" << f_time << "\t" <<
    undershoot << "\n\n" << std::endl; 

  if (fgaus != nullptr){
    spe_ampl_correction = compute_spe_correction(fgaus);
    spe_ampl /= spe_ampl_correction;
    std::cout << "\nSpe ampl corrected = " << spe_ampl << std::endl;
  }



}
