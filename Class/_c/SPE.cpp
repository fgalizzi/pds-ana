
void cla::SPE() {
  double y0, y1, r_time, f_time, undershoot;
  vector<double> x, avg, int_wf;
  vector<vector<double>> y2, y3, avg_wf;
    
  for (size_t i = 0; i < memorydepth; i++) x.push_back( (double) i*tick_len);
  
  // Read and subtract the baseline
  read();

  SelCalib_WF(wfs, y2, prepulse_ticks, sat_low, sat_up, bsl);
  TH1D* hI = BuildRawChargeHisto(y2, int_wf, int_low, int_up, nbins);
  if(mov_win == true) MovingAverageWF(y2, y2, win);
  Avg_Sel_WF (y2, y3, avg, int_wf, spe_low, spe_up);
 
  avg_wf.push_back(avg);
  RiseFallTimeUndershoot(avg, tick_len, r_time, f_time, undershoot);
  spe_ampl = *std::max_element(avg.begin(), avg.end());
  //MovingAverageWF(avg_wf, avg_wf, 100);

  if(print==true) VecDouble_in_Binary("Template.dat", avg);

  double sum = 0;
  for(auto e : avg ) sum+=e;
  std::cout << "SPE integral [0;memorydepth] " << sum << std::endl;

  TGraph* gAvg= build_avg_spectral_density(memorydepth, tick_len*memorydepth, tick_len, avg_wf, res);
  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,900);
  c2->cd();
  gAvg->Draw();
  c2->Modified();
  c2->Update();

  TGraph *g1 = new TGraph(avg.size(), &x[0], &avg[0]);
  g1->GetXaxis()->SetTitle("Time [#mus]");
  g1->GetYaxis()->SetTitle("ADC counts");
  g1->SetLineColor(2);
  TCanvas *c1 = new TCanvas("c1","c1",20,20,1000,900);
  c1->cd();
  g1->Draw();
  c1->Modified();
  c1->Update();

  std::cout << "\n\nColdbox table Ampl - under - r_time - f_time - under" << std::endl;
  std::cout << spe_ampl << "\t" << undershoot*spe_ampl/100. << "\t" << r_time << "\t" << f_time << "\t" <<
    undershoot << "\n\n" << std::endl; 

}
