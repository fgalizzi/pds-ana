
void cla::Noise_PSD(){
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  gStyle->SetPalette(kSunset);
  
  double t;
  vector<double> x, avg, int_wf;
  vector<vector<double>> noise, noise2, avg_wf;
    
  for (size_t i = 0; i < memorydepth; i++) x.push_back( (double) i);

  // Read and subtract the baseline
  read();
  
  SelCalib_WF(wfs, noise, prepulse_ticks, -bsl, bsl, bsl);
  TH1D* hI = BuildRawChargeHisto(noise, int_wf, int_low, int_up, nbins);
  Avg_Sel_WF (noise, noise2, avg, int_wf, mu0_low, mu0_up);

  avg_wf.push_back(avg);
  
  TGraph* gAvg= build_avg_spectral_density(memorydepth, tick_len*memorydepth, tick_len, avg_wf, res);
  gAvg->SetLineColor(2);
  TGraph* gNoise_spectral_density = build_avg_spectral_density(memorydepth,
      tick_len*memorydepth, tick_len, noise2, res);


  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,900);
  c2->cd();
  //gNoise_spectral_density->Draw();
  gAvg->Draw();
  c2->Modified();
  c2->Update();
  
  TGraph *g1 = new TGraph(avg.size(), &x[0], &avg[0]);
  g1->GetXaxis()->SetTitle("Time [#mus]");
  g1->GetYaxis()->SetTitle("ADC counts");
  g1->SetLineColor(2);

  TCanvas *c3 = new TCanvas("c3","c3",30,30,1000,900);
  c3->cd();
  g1->Draw();
  c3->Modified();
  c3->Update();

  if(print==true){
    std::ofstream OutFile  ("FFT_ch.dat", ios::binary);
    std::ofstream OutFile2 ("FFT_ch_avg.dat", ios::binary);
    for(int i = 0; i < gNoise_spectral_density->GetN(); i++){
    t = gNoise_spectral_density->GetPointY(i);
    OutFile.write(reinterpret_cast<char*> (&t), sizeof(t));
    t = gAvg->GetPointY(i);
    OutFile2.write(reinterpret_cast<char*> (&t), sizeof(t));
    t = gNoise_spectral_density->GetPointX(i);
    OutFile.write(reinterpret_cast<char*> (&t), sizeof(t));
    OutFile2.write(reinterpret_cast<char*> (&t), sizeof(t));
    }   
    std::cout << "Vector saved in ---> " << std::endl;
    OutFile.close();
    OutFile2.close();
  }
  
}
 
