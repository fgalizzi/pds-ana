#include "classe.hpp"

using namespace std;

//*********************************************
TF1* cla::set_charge_fit_function(TH1D* hI, TH1D* hFind, bool avoid_auto_peak){
//*********************************************
  if(hFind==nullptr) TH1D* hFind = new TH1D(*hI);
  // Go for peaks: create an instance of TSpectrum
  hFind->Rebin(8);
  TSpectrum *s = new TSpectrum(5);
  s->Background(hFind, 50, "same");
  int npeaks = s->Search(hFind,1,"",0.2);
  Double_t *xpeaks;
  xpeaks = s->GetPositionX();

  std::vector<double> xp;
  for(int i=0; i<npeaks; i++) xp.push_back(xpeaks[i]);
  std::sort(xp.begin(), xp.end());
  
  if (xp.size()<2 || avoid_auto_peak==true){
    cout << "\nWARNING: using spe_low - spe_up\n" << endl; 
    xp = {0, (spe_low+spe_up)/2.};
  }
  
  // Parameters for the fit function
  double par[20] = {0}; //Norm_Consts, sigma_0 , sigma_cell, G
  
  // Fit on the first peak to estimate sigma_0
  TF1* f0 = new TF1("f0", "gaus", -xp[1]*0.35, xp[1]*0.2);
  f0->SetParameters(hI->GetBinContent(hI->GetMaximumBin()), xp[0], xp[1]*0.5);
  hI->Fit("f0", "R");

  // Initial guess using TSpectrum
  par[0] = xp[0];               //peak 0
  par[1] = xp[1]-xp[0];         //peak 1
  par[2] = f0->GetParameter(2); //sigma_0
  par[3] = par[2]*0.4;          //sigma_cell
  
  // Set the initial guess manually 
  if(manual==true){
    par[0] = 0.;                //peak 0
    par[1] = (spe_low+spe_up)/2;//peak 1
    par[2] = (s0_low+s0_up)/2;  //sigma_0
    par[3] = (sc_low+sc_up)/2;  //sigma_cell
  }
 // if(hI->GetBinContent(hI->GetMaximumBin())<70) hI->Rebin(2);
  if(npeaks<nmaxpeaks) npeaks=nmaxpeaks;
  for(int i = 0 ; i < npeaks ; i++){
    par[i+4] = hI->GetBinContent(hI->GetMaximumBin())*0.8;//peaky[i];
  }

  // Gaussian sumatory
  fMultiGauss f_multi_gauss;
  f_multi_gauss.npeaks = npeaks;
  TF1* ff = new TF1("fgaus", f_multi_gauss, -par[2]*2., par[1]*(npeaks-1), npeaks+4);
  fMultiGauss_set(ff);
  ff->SetParameters(par);
  ff->SetParLimits(0, -par[1]*0.2, par[1]*0.2);
  ff->SetParLimits(1, par[1]*0.4, par[1]*1.9);
  ff->SetParLimits(2, par[2]*0.3, par[2]*2.);
  ff->SetParLimits(3, 0., par[2]*2.);
  
  // Set limits manually
  if(manual==true){
    ff->SetParLimits(0,mu0_low,mu0_up);
    ff->SetParLimits(1,spe_low,spe_up);
    ff->SetParLimits(2,s0_low,s0_up);
    ff->SetParLimits(3,sc_low,sc_up);
  }
  
  for(int i = 0 ; i < npeaks; i++) ff->SetParLimits(i+4,0,5000);
  return ff;
}

// Compute the correction factor for spe amplitude estimate
// considering that there could be 0pes and 2pes wfs among the
// selected spe candidate 
double cla::compute_spe_correction(TF1* f){
  double mu0 = f->GetParameter(0);
  double s0 = f->GetParameter(2);
  double sc = f->GetParameter(3);
  double s1 = sqrt(s0*s0+sc*sc);
  double s2 = sqrt(s0*s0+2*sc*sc);
  
  TF1* gaus0 = new TF1("gaus0", "gaus", -10.*s0, 10.*s0);
  gaus0->SetParameters(f->GetParameter(4), mu0, s0);
  TF1* gaus2 = new TF1("gaus2", "gaus", mu0+spe_charge-10.*s2, mu0+spe_charge+10.*s2);
  gaus2->SetParameters(f->GetParameter(6), mu0+2*spe_charge, s2);
 
  double ll = mu0+spe_charge-s1;
  double ul = mu0+spe_charge+s1;
  double f_corr = f->Integral(ll,ul)-gaus0->Integral(ll,ul)+gaus2->Integral(ll,ul);
  f_corr /= f->Integral(ll,ul); 

  return f_corr;
}

//*********************************************
void cla::self_histos(TH1D* h_all, TH1D* h_trg, std::vector<double>& int_wf){
//*********************************************
  int len = wfs[0].size();
  int count = 0;
  std::vector<int> trg_bool;
  double spe_norm = 1./spe_charge;
  double t;
  std::vector<double> trgs;
  int got_ya;
   
  TH1D* hI = new TH1D();
  
  if(apply_filter == 0){
    hI = BuildRawChargeHisto(wfs, int_wf, int_low, int_up, nbins);
  } 
  else{
    hI = BuildRawChargeHisto(wfs, int_wf, 0, 1, nbins);
  }
 
  // Spectra of true positive
  for(auto wf : trg_wf){
    trgs = TriggerTime(wf);
    got_ya = 0; //false
    for (auto tr : trgs){
      if (tr > pretrg && tr < afttrg) got_ya = 1; //true
    }
    trg_bool.push_back(got_ya);
  }

  for (size_t i=0; i<int_wf.size(); i++){
    t = (int_wf[i]-pedestal)*spe_norm;
    h_all->Fill(t);
    if (trg_bool[i] == 1) h_trg->Fill(t);
  }
  std::cout << "#Self-Trigger in coincedence with LED " << count << std::endl;
 
}


void cla::remove_saturated_wfs(){
  std::vector<std::vector<double>> wfs_new;
  for(auto wf : wfs){
    bool sat = false;
    for(auto sample : wf){
      if(sample >= pow(2,res)-1 || sample == 0){
        sat = true;
        break;
      }
    }
    if(!sat) wfs_new.push_back(wf);
  }
  wfs = wfs_new;
  n_wf = wfs.size();
}

void cla::remove_too_little_wfs(double peak_to_peak){
  std::vector<std::vector<double>> wfs_new;
  for(auto wf : wfs){
    double max_el = *max_element(wf.begin(), wf.end());
    double min_el = *min_element(wf.begin(), wf.end());
    if((max_el - min_el) >= peak_to_peak){
      wfs_new.push_back(wf);
    }
  }
  wfs = wfs_new;
  n_wf = wfs.size();
}
