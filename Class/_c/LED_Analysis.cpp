//*** MAIN ************************************
void cla::LED_Analysis(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  vector<vector<double>> sel_wf; // Wavefroms I select for the analysis
  vector<double> int_wf;         // Integrals of the selected waveforms

  // Read and subtract the baseline
  read();

  // Select the with a baseline within the [-bsl; bsl] range and 
  // fully contained in the [sat_low; sat_up] range
  SelCalib_WF(wfs, sel_wf, prepulse_ticks, sat_low, sat_up, bsl);  
  
  TH1D* hI = new TH1D();

  if (manual == 0) {
    if(mov_win == true){
      MovingAverageWF(sel_wf, sel_wf, win);   
      hI = BuildRawChargeHisto(sel_wf, int_wf, int_low+win, int_up+win, nbins);
    } else hI = BuildRawChargeHisto(sel_wf, int_wf, int_low, int_up, nbins);
  }
  else {
    if(mov_win == true){
      MovingAverageWF(sel_wf, sel_wf, win);   
      hI = BuildRawChargeHisto(sel_wf, int_wf, int_low+win, int_up+win,
        hmin, hmax, nbins);
    } else hI = BuildRawChargeHisto(sel_wf, int_wf, int_low, int_up,
        hmin, hmax, nbins);
  }
  
  TH1D* hFind = new TH1D(*hI); //Delcared to plot it later
 
  // See private_methods.hpp
  // "fgaus" is the function to fit the chargehistogram, it is initialized 
  // by looking for peaks in the histo or manually by setting manual=true
  fgaus = set_charge_fit_function(hI, hFind);
  
  TCanvas *c3 = new TCanvas("c3","c3",0,0,500,450);
  c3->cd();
  hFind->Draw();
  c3->Modified();
  c3->Update();
 

  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,900);
  c1->cd();
  hI->Draw();// Fit histogram
  cout << "Fit ... " << endl;
  hI->Fit("fgaus", "R");
  cout << "... end fit. " << endl;
  
  c1->Modified();
  c1->Update();
 
  double SNR, q1, q1q2, s0;
  SNR = fgaus->GetParameter(1)/fgaus->GetParameter(2);
  q1  = fgaus->GetParameter(1)+fgaus->GetParameter(0);
  q1q2= fgaus->GetParameter(1);
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
