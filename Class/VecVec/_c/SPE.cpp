
void cla::SPE() {
  double y0, y1;
  vector<double> x, y2, y3, avg_wf, int_wf;
    
  for (size_t i = 0; i < MEMORYDEPTH; i++) x.push_back( (double) i*TICK_LEN);
  
  // Read and subtract the baseline
  read();

  SelCalib_WF(WFS, y2, MEMORYDEPTH, PREPULSE_TICKS, SAT_LOW, SAT_UP, BSL);
  TH1D* hI = BuildRawChargeHisto(y2, int_wf, MEMORYDEPTH, INT_LOW, INT_UP);
  Avg_Sel_WF (y2, y3, avg_wf, int_wf, SPE_LOW, SPE_UP);
 
  Print_RiseFallTime(avg_wf, TICK_LEN);
  //MovingAverageWF(avg_wf, avg_wf, 100);

  if(PRINT==true) VecDouble_in_Binary("SPE.dat", avg_wf);

  TGraph* gAvg= build_avg_spectral_density(MEMORYDEPTH, TICK_LEN*MEMORYDEPTH, TICK_LEN, avg_wf, RES);
  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,900);
  c2->cd();
  gAvg->Draw();
  c2->Modified();
  c2->Update();

  TGraph *g1 = new TGraph(avg_wf.size(), &x[0], &avg_wf[0]);
  g1->GetXaxis()->SetTitle("Time [#mus]");
  g1->GetYaxis()->SetTitle("ADC counts");
  g1->SetLineColor(2);
  TCanvas *c1 = new TCanvas("c1","c1",20,20,1000,900);
  c1->cd();
  g1->Draw();
  c1->Modified();
  c1->Update();

}
