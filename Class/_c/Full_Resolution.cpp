#include "../classe.hpp"

TH2D* X_Histo(vector<vector<double>>& all_wf, int I_low, int I_up,
             int L_up, int window, int x_bins, int y_bins){
  int len = all_wf[0].size();
  int n_wf = all_wf.size();
  std::vector<double> charge, pp, d_charge, d_pp;
  charge.resize(n_wf, 0.);
  pp.resize(n_wf, 0.);
  vector<double> t_wf;
  
  int i=0;
  for(auto wf : all_wf) {
    auto in  = wf.begin()+I_low;
    auto fin = wf.begin()+I_up;
    charge[i] = std::reduce(in, fin, 0.0);
    i++;
  }

  MovingAverageWF(all_wf, all_wf, window); 
  
  i = 0;
  for(auto wf : all_wf){
    auto in  = wf.begin()+I_low+window;
    auto fin = wf.begin()+I_up+window;
    double max = *std::max_element(in, fin);
    double min = *std::min_element(in, fin-I_up+L_up);
    pp[i] = 30*(max-min);
    
    i++;
  }
  
  //Dual vectors
  d_charge = charge; d_pp = pp;
  sort(d_charge.begin(), d_charge.end());
  sort(d_pp.begin(), d_pp.end());
  double x_low = d_charge[size_t(d_charge.size()*0.01)];
  double x_up  = d_charge[size_t(d_charge.size()*0.99)];
  double y_low = d_pp[size_t(d_charge.size()*0.01)];
  double y_up  = d_pp[size_t(d_charge.size()*0.99)];
  TH2D* h2 = new TH2D("h2", Form("%s,%s,%s", "X-Histo", "Charge", "Peak-Peak"), 
                      x_bins, x_low, x_up, int(y_up-y_low), -1.1*y_low, 1.1*y_up);
  //h1 = new TH1D("h", "h", 500, 10*x_low*y_low, x_up*y_up);
  for (size_t i = 0; i<d_charge.size(); i++) h2->Fill(charge[i], pp[i]);
  TF1* f = new TF1("ff", "pol1", x_low, x_up);
  h2->Fit(f, "R");
  h2->Draw(); gPad->Update(); gPad->WaitPrimitive();
  double theta = atan(f->GetParameter(1));
  double c_t = cos(theta);
  double s_t = sin(theta);

  for (size_t i = 0; i<d_charge.size(); i++){
    charge[i] = c_t*charge[i]+s_t*pp[i];
    pp[i]     = -s_t*charge[i]+c_t*pp[i];
  }
 

  d_charge = charge; d_pp = pp;
  sort(d_charge.begin(), d_charge.end());
  sort(d_pp.begin(), d_pp.end());
  x_low = d_charge[size_t(d_charge.size()*0.01)];
  x_up  = d_charge[size_t(d_charge.size()*0.99)];
  y_low = d_pp[size_t(d_charge.size()*0.01)];
  y_up  = d_pp[size_t(d_charge.size()*0.99)];
  h2 = new TH2D("h2", Form("%s,%s,%s", "X-Histo", "Charge", "Peak-Peak"), 
                      x_bins, x_low, x_up, int(y_up-y_low), -1.1*y_low, 1.1*y_up);

  for (size_t i = 0; i<d_charge.size(); i++) h2->Fill(charge[i], pp[i]);
  f = new TF1("f2", "pol1", x_low, x_up);
  h2->Fit(f, "R");
  return h2;
}


//*** MAIN ************************************
void cla::Full_Resolution(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  vector<vector<double>> sel_wf;
  vector<double> int_wf;

  // Read and subtract the baseline
  read();

  //SubBaseline(wfs, memorydepth, prepulse_ticks);
  SelCalib_WF(wfs, sel_wf, prepulse_ticks, sat_low, sat_up, bsl);  
  
  TH2D* h2 = X_Histo(sel_wf, int_low, int_up, int_up+100, win, 200, 200);
  h2->Draw("COLZ");




















  //TH1D* hI = BuildRawChargeHisto(sel_wf, int_wf, memorydepth, int_low, int_up);
  
/*  TH1D* hFind = new TH1D(*hI);

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
    par[1] = (spe_low+spe_up)/2;//peak 1
    par[2] = (s0_low+s0_up)/2;  //sigma_0
    par[3] = (sc_low+sc_up)/2;  //sigma_cell
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
  */
}
