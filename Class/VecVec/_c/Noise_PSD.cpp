
void cla::Noise_PSD(){
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  gStyle->SetPalette(kSunset);
  
  double t;
  std::vector<Double_t> x, avg_wf, noise2, noise, y3, int_wf;
    
  for (size_t i = 0; i < MEMORYDEPTH; i++) x.push_back( (double) i);

  // Read and subtract the baseline
  read();

  std::cout << WFS.size() << std::endl;
  //SelCalib_WF(WFS, noise2, MEMORYDEPTH, PREPULSE_TICKS, SAT_LOW, BSL, BSL);
  TH1D* hI = BuildRawChargeHisto(WFS, int_wf, MEMORYDEPTH, INT_LOW, INT_UP);
  Avg_Sel_WF (WFS, y3, avg_wf, int_wf, MU0_LOW, MU0_UP);

  //std::cout << "N wf for FFT : " << nnn << std::endl;

  TGraph* gAvg= build_avg_spectral_density(MEMORYDEPTH, TICK_LEN*MEMORYDEPTH, TICK_LEN, avg_wf, RES);
  gAvg->SetLineColor(2);
  TGraph* gNoise_spectral_density = build_avg_spectral_density(MEMORYDEPTH, TICK_LEN*MEMORYDEPTH, TICK_LEN, y3, RES);

  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,900);
  c2->cd();
  gNoise_spectral_density->Draw();
  gAvg->Draw("same");
  c2->Modified();
  c2->Update();

  TGraph *g1 = new TGraph(avg_wf.size(), &x[0], &avg_wf[0]);
  g1->GetXaxis()->SetTitle("Time [#mus]");
  g1->GetYaxis()->SetTitle("ADC counts");
  g1->SetLineColor(2);

  TCanvas *c3 = new TCanvas("c3","c3",30,30,1000,900);
  c3->cd();
  g1->Draw();
  c3->Modified();
  c3->Update();
  

  if(PRINT==true){
    std::ofstream OutFile ("FFT_file.dat", ios::binary);
    for(int i = 0; i < gNoise_spectral_density->GetN(); i++){
    t = gNoise_spectral_density->GetPointY(i);
    OutFile.write(reinterpret_cast<char*> (&t), sizeof(t));
    t = gNoise_spectral_density->GetPointX(i);
    OutFile.write(reinterpret_cast<char*> (&t), sizeof(t));
    }   
    std::cout << "Vector saved in ---> " << std::endl;
    OutFile.close();
  }
  
}
 
