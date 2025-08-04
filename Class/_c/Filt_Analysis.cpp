#include "../classe.hpp"

// ****************************************************************
// Description
// ****************************************************************

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Filt_Analysis(){

  vector<vector<double>> sel_wf, filt_wf;
  vector<double> int_wf, t_templ, f_noise;
  size_t nsample = memorydepth;
  TComplex G[nsample];
  string filt =  "matched";

  // Read and subtract the baseline
  read(); // wfs to analyze
  CompleteWF_Binary(templ_f, t_templ, memorydepth); // t_templ = time domain template
  SelCalib_WF(wfs, sel_wf, prepulse_ticks, sat_low, sat_up, bsl);  

  // Build the filter G in freq. domain
  if (filt == "matched"){
    Build_Matched_Filter(&G[0], t_templ);
  } /*else if (filt = "wiener") {
    CompleteWF_Binary(noise_f, f_noise, memorydepth); // f_noise = noise FFT
    Build_Wiener_Filter(G, t_templ, f_noise); 
  } */
  else {
    std::cout << "\n\n This filter is not implemented: " << filt << std::endl;
    return;
  }

  FilterAllWF(sel_wf, filt_wf, G);
  
  if (ndisplay>0) DisplayWFs(filt_wf, 1., ndisplay);

  TH1D* hI = new TH1D();

  hI = BuildRawChargeHisto(filt_wf, int_wf, 0, 1, nbins);
  
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
  if(manual==true){
    par[0] = 0.;                //peak 0
    par[1] = double(spe_low+spe_up)/2;//peak 1
    par[2] = double(s0_low+s0_up)/2;  //sigma_0
    par[3] = double(sc_low+sc_up)/2;  //sigma_cell
  }
 
  if(hI->GetBinContent(hI->GetMaximumBin())<70) hI->Rebin(2);
  if(npeaks<nmaxpeaks) npeaks=nmaxpeaks;
  for(int i = 0 ; i < npeaks ; i++){
    par[i+4] = hI->GetBinContent(hI->GetMaximumBin())*0.8;//peaky[i];
  }
  

  // Gaussian sumatory
  fgaus = new TF1("fgaus", fRandomName, -par[2]*2., par[1]*(npeaks-1), npeaks+4);
  if(manual==true){
    fgaus = new TF1("fgaus", fRandomName, fit_low, fit_up, npeaks+4);
  }
  fRandomName_set(fgaus);
  fgaus->SetParameters(par);
  fgaus->SetParLimits(0, -par[1]*0.2, par[1]*0.2);
  fgaus->SetParLimits(1, par[1]*0.4, par[1]*1.9);
  fgaus->SetParLimits(2, par[2]*0.3, par[2]*2.);
  fgaus->SetParLimits(3, 0., par[2]*2.);
  
  // Set limits manually
  if(manual==true){
    fgaus->SetParLimits(0,mu0_low,mu0_up);
    fgaus->SetParLimits(1,spe_low,spe_up);
    fgaus->SetParLimits(2,s0_low,s0_up);
    fgaus->SetParLimits(3,sc_low,sc_up);
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
 
  double SNR, q1, q1q2, s0;
  SNR = fgaus->GetParameter(1)/fgaus->GetParameter(2);
  q1  = fgaus->GetParameter(1)+fgaus->GetParameter(0);
  q1q2= fgaus->GetParameter(1);
  pedestal = q1q2;
  spe_charge = q1;
  s0  = fgaus->GetParameter(2);
  //std::cout << "\n\nColdbox table SNR - q1 - q1q2 - s0" << std::endl;
  //std::cout << SNR << "\t" << q1 << "\t" << q1q2 << "\t" << s0 << "\n\n" << std::endl; 
  std::cout << "\n\nColdbox June: Gain - S0 - SNR" << std::endl;
  std::cout << q1q2 << "\t" << s0 << "\t" << SNR << "\n\n" << std::endl; 
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
