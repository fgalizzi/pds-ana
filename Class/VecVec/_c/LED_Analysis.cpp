//*** MAIN ************************************
void cla::LED_Analysis(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  std::vector<double> sel_wf, int_wf;

  // Read and subtract the baseline
  read();

  //SubBaseline(WFS, MEMORYDEPTH, PREPULSE_TICKS);
  SelCalib_WF(WFS, sel_wf, MEMORYDEPTH, PREPULSE_TICKS, SAT_LOW, SAT_UP, BSL);  
  
  //TH1D* hI = BuildRawChargeHisto(sel_wf, int_wf, MEMORYDEPTH, INT_LOW, INT_UP);
  TH1D* hI = BuildRawChargeHisto(sel_wf, int_wf, MEMORYDEPTH, INT_LOW, INT_UP,
      HMIN, HMAX, NBINS);
  TH1D* hFind = new TH1D(*hI);

  // Go for peaks: create an instance of TSpectrum
  hFind->Rebin(8);
  TSpectrum *s = new TSpectrum(5);
  s->Background(hFind, 50, "same");
  int npeaks = s->Search(hFind,1,"",0.2);
  Double_t *xpeaks;
  xpeaks = s->GetPositionX();
  //hFind->Draw();s->Draw();gPad->Update();gPad->WaitPrimitive();

  std::vector<double> xp;
  for(int i=0; i<npeaks; i++) xp.push_back(xpeaks[i]);
  std::sort(xp.begin(), xp.end());
  std::cout << "Peaks' positions:" << std::endl;
  for(auto xx:xp) std::cout << xx << std::endl;

  // Parameters for the fit function
  double par[20] = {0}; //Norm_Consts, sigma_0 , sigma_cell, G
  
  // Fit on the first peak to estimate sigma_0
  TF1* f0 = new TF1("f0", "gaus", -xp[1]*0.3, xp[1]*0.3);
  f0->SetParameters(hI->GetBinContent(hI->GetMaximumBin()), xp[0], xp[1]*0.5);
  hI->Fit("f0", "R");
  
  // Initial guess using TSpectrum
  par[0] = xp[0];               //peak 0
  par[1] = xp[1]-xp[0];         //peak 1
  par[2] = f0->GetParameter(2); //sigma_0
  par[3] = par[2]*0.1;          //sigma_cell
  
  // Set the initial guess manually 
  if(MANUAL==true){
    par[0] = 0.;                //peak 0
    par[1] = (SPE_LOW+SPE_UP)/2;//peak 1
    par[2] = (S0_LOW+S0_UP)/2;  //sigma_0
    par[3] = (SC_LOW+SC_UP)/2;  //sigma_cell
  }
 
  if(hI->GetBinContent(hI->GetMaximumBin())<70) hI->Rebin(2);
  if(npeaks<NMAXPEAKS) npeaks=NMAXPEAKS;
  for(int i = 0 ; i < npeaks ; i++){
    par[i+4] = hI->GetBinContent(hI->GetMaximumBin())*0.8;//peaky[i];
  }
  

  // Gaussian sumatory
  TF1* fgaus = new TF1("fgaus", fRandomName, -par[2]*2., par[1]*(npeaks-1), npeaks+4);
  if(MANUAL==true){
    fgaus = new TF1("fgaus", fRandomName, FIT_LOW, FIT_UP, npeaks+4);
  }
  fRandomName_set(fgaus);
  fgaus->SetParameters(par);
  fgaus->SetParLimits(0, -par[1]*0.2, par[1]*0.2);
  fgaus->SetParLimits(1, par[1]*0.4, par[1]*1.9);
  fgaus->SetParLimits(2, par[2]*0.3, par[2]*2.);
  fgaus->SetParLimits(3, 0., par[2]*2.);
  
  // Set limits manually
  if(MANUAL==true){
    fgaus->SetParLimits(0,MU0_LOW,MU0_UP);
    fgaus->SetParLimits(1,SPE_LOW,SPE_UP);
    fgaus->SetParLimits(2,S0_LOW,S0_UP);
    fgaus->SetParLimits(3,SC_LOW,SC_UP);
  }
  
  for(int i = 0 ; i < npeaks ; i++) fgaus->SetParLimits(i+4,0,2700);

  TCanvas *c3 = new TCanvas("c3","c3",0,0,500,450);
  c3->cd();
  hFind->Draw();
  c3->Modified();
  c3->Update();
 
  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,900);
  c1->cd();

  // Fit histogram
  hI->Draw();
  cout << "Fit ... " << endl;
  hI->Fit("fgaus", "R");
  cout << "... end fit. " << endl;
  
  c1->Modified();
  c1->Update();
 

  return; 
  //fit CX
  auto g_CX = Build_CX_Graph(fgaus, hI);
  TF1* f_CX = new TF1("f_CX", fCX, -0.5, 5.5, 2);
  fCX_set(f_CX);
  
  g_CX->Fit("f_CX", "R");
  
  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,900);
  c2->cd();
  g_CX->Draw();
  c2->Modified();
  c2->Update();
  
}
