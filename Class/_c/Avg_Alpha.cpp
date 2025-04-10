#include "../classe.hpp"
#include <cstdio>
#include <numeric>
#include <string>

// ****************************************************************
// This macro is used to create an average alpha-waveform and save it
// in a binary file, setting print=true, once you are happy with the
// selection cuts. The idea is to select the alpha waveforms on an
// f_prompt basis and remove the ones that have light in coincidence.
//
// Parameters you can tune:
// - int_low, int_prompt, int_up <-> integration limits
// - f_prompt (integral[int_low;int_prompt]/integral[int_low;int_up])
//      We select waveforms with f_prompt lower than this value.
// - sat_low <-> discard wfs with samples below this level
// - amp_low, amp_up <-> the max element in the [int_low;int_prompt]
//      range must fall in [amp_low;amp_up]
// - bsl <-> accept only wfs with sample whithin [-bsl;+bsl] in the
//      pre trigger [0;prepulse_ticks]
// - rms <-> the selected wfs can differ from the average only by random
//      fluctuation, we must ensure there are no pulses in the tail of
//      our signal.
// ****************************************************************

////////////////////////////////////////////////////////
/////// HARD CODE //////////////////////////////////////

std::string alpha_files_path = "/eos/home-f/fegalizz/PDE_MiB/PDE_Results/Alpha_files/";
double nl_sigma = 2.5; // Lower selection = mu - nl_sigma*sigma
double nr_sigma = 0.5; // Upper selection = mu + nr_sigma*sigma


////////////////////////////////////////////////////////

Double_t fAlphaPeak(Double_t *x,Double_t *par){
  Double_t fitval = par[0]/2/par[1]* TMath::Exp( (x[0]- par[2])/par[1] + par[3]*par[3]/2/par[1]/par[1]) * TMath::Erfc(1/sqrt(2)* ((x[0]- par[2])/par[3] + par[3]/par[1]) );
  return fitval;
}

void FitAlpha_Spectrum(TH1D* h_charge, double& mu, double& sigma){
  Double_t fAlphaPeak(Double_t *x,Double_t *par);
  mu = h_charge->GetXaxis()->GetBinCenter(h_charge->GetMaximumBin());
  float xmax = h_charge->GetXaxis()->GetXmax();
  TF1* fgaus = new TF1("fgaus","gaus",mu,xmax);
  h_charge->Fit("fgaus","QNR");
  mu = fgaus->GetParameter(1); 
  sigma = fgaus->GetParameter(2);
  float xmin = mu - 2*sigma;
  xmax = mu + 8*sigma;
  int npar = 4;
  TF1 *func = new TF1("fAlphaPeak", fAlphaPeak, xmin, xmax, npar);
  func->SetParNames ("Amplitude", "Tau", "mu", "sigma");
  func->SetParameters(3e8, 2e6, mu, sigma);

  // faccio un primo fit per avere una stima corretta di mu e sigma
  h_charge->Fit("fAlphaPeak","QNR");

  // ora correggo il range e faccio il fit vero e proprio
  xmin = func->GetParameter("mu") - 3.*func->GetParameter("sigma");
  xmax = func->GetParameter("mu") + 8.*func->GetParameter("sigma");
  func->SetRange(xmin, xmax);
  h_charge->Fit("fAlphaPeak","RMES");
  mu = func->GetParameter("mu"); sigma = func->GetParameter("sigma");
}

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Avg_Alpha(){
  std::vector<double> avg_alpha, int_wf;
  std::vector<std::vector<double>> alpha_wfs, sel_wfs;

  // Read, subtract the baseline, select wfs according to bsl in the
  // pre-trigger, build and fit the alpha spectrum
  read();
  SelCalib_WF(wfs, sel_wfs, prepulse_ticks, -1e6, 1e6, bsl);  
  h_charge = BuildRawChargeHisto(sel_wfs, int_wf, int_low, int_up, nbins);
  double mu, sigma;
  FitAlpha_Spectrum(h_charge, mu, sigma);

  //-----------------------------------------------------------------
  //----- SELECTION LOOPS -------------------------------------------
  double max_el, min_el, t, integral;
  size_t len = sel_wfs[0].size();
  size_t nwfs = sel_wfs.size();
  int bsl_counter = 0;
  int int_counter = 0;
  int fpr_counter = 0;
  int sel_counter = 0;

  vector<double> first_avg(len, 0.0);
  vector<bool> preliminary_selection(nwfs, false);
  vector<bool> final_selection(nwfs, false);

  alpha_wfs.reserve(nwfs);
  std::cout << "Selection range " << mu-nl_sigma*sigma << " " << mu+nr_sigma*sigma << std::endl;
  std::cout << "Mu " << mu << std::endl;
  std::cout << "Sigma " << sigma << std::endl;
    
  for (size_t i=0; i<nwfs; i++) {
    if (int_wf[i] < mu+nl_sigma*sigma && int_wf[i] > mu-nr_sigma*sigma) {
      int_counter++;
      t = *max_element(sel_wfs[i].begin()+int_low, sel_wfs[i].begin()+int_up);
      max_el = *max_element(sel_wfs[i].begin()+int_up, sel_wfs[i].end());

      if (max_el < 0.5*t) {
        sel_counter++;
        for(size_t j=0; j<len; j++) first_avg[j] += sel_wfs[i][j];
        preliminary_selection[i]=true;
      }
    }
  }  
 
  double norm = 1./ *max_element(first_avg.begin(), first_avg.end());
  for (auto& e : first_avg) e*=norm;

  for (size_t i=0; i<nwfs; i++) {
    if(preliminary_selection[i]==true){
      vector wf = sel_wfs[i];
      norm = 1./ *max_element(wf.begin(), wf.end());
      double norm_bsl = rms*norm*5.;
      for (auto& e : wf) e *= norm;

      bool select_shape = true;
      for (size_t j=0; j<len; j++){
        if (wf[j] > first_avg[j]+norm_bsl || wf[j] < first_avg[j]-norm_bsl){
          select_shape = false;
          j = len+1;
        }
      }
      final_selection[i]=select_shape;
    }

    if(final_selection[i]==true) alpha_wfs.push_back(sel_wfs[i]);
  }
  
  cout << "WFs selected on the baseline "  << nwfs << "/" << wfs.size() << endl;
  cout << "WFs selected on the integral " << int_counter << "/" << wfs.size() << endl;
  cout << "WFs preliminary selected " << sel_counter << "/" << wfs.size() << endl;
  cout << "WFs selected " << alpha_wfs.size() << "/" << wfs.size() << endl;

 
  //----- SELECTION LOOPS END ---------------------------------------
  //-----------------------------------------------------------------




  cout << "\nAlpha canidates: " << alpha_wfs.size() << endl;
  avgWF(alpha_wfs, avg_alpha);
  RiseFallTimeUndershoot(avg_alpha, tick_len, int_up);

  if(print==true) {
    std::cout << wf_file << std::endl;
    string outfile_name;
    cout << "Name of the alpha file (without .dat)" << endl;
    cin >> outfile_name;
    outfile_name = alpha_files_path+outfile_name+".dat";
    VecDouble_in_Binary(outfile_name, avg_alpha);
    print = false;
  }




  //-----------------------------------------------------------------
  //----- PLOTS -----------------------------------------------------
  if (plot == true){
    double ymin, ymax;
    // All wfs persistence ------------------------------------------ 
    min_max_element(sel_wfs, ymin, ymax);
    TH2D* h2_sel_pers = new TH2D("h2_sel_pers", Form("%s;%s;%s", "Persistence", "Ticks", "ADC Counts"),
        memorydepth/2, 0., memorydepth, 
        120, ymin, ymax);

    for (size_t i=0; i< sel_wfs.size(); i++)
      for (int j=0; j<memorydepth; j=j+2) h2_sel_pers->Fill(j, sel_wfs[i][j]);

    TCanvas *c_sel_pers = new TCanvas("c_sel_pers","c_sel_pers",0,0,500,450);
    c_sel_pers->cd();
    h2_sel_pers->Draw("COLZ");
    c_sel_pers->SetLogz();
    c_sel_pers->Modified(); c_sel_pers->Update();
    
    // Alpha wfs persistence ----------------------------------------
    min_max_element(alpha_wfs, ymin, ymax);
    TH2D* h2_alpha_pers = new TH2D("h2_alpha_pers", Form("%s;%s;%s", "Persistence", "Ticks", "ADC Counts"),
        memorydepth/2, 0., memorydepth, 
        120, ymin, ymax);
    
    for (size_t i=0; i< alpha_wfs.size(); i++) 
      for (int j=0; j<memorydepth; j=j+2) h2_alpha_pers->Fill(j, alpha_wfs[i][j]);
    
    TCanvas *c_alpha_pers = new TCanvas("c_alpha_pers","c_alpha_pers",500,0,500,450);  
    c_alpha_pers->cd();
    h2_alpha_pers->Draw("COLZ");
    c_alpha_pers->SetLogz();
    c_alpha_pers->Modified(); c_alpha_pers->Update();

    // Alpha average ------------------------------------------------
    TGraph *g_alpha_avg = new TGraph(avg_alpha.size(), &avg_alpha[0]);
    //g_alpha_avg->GetXaxis()->SetTitle("Time [#mus]");
    g_alpha_avg->GetYaxis()->SetTitle("ADC counts");
    g_alpha_avg->SetLineColor(2);

    TCanvas *c_alpha_avg = new TCanvas("c_alpha_avg","c_alpha_avg",0,550,500,450);
    c_alpha_avg->cd();
    g_alpha_avg->Draw();
    c_alpha_avg->Modified(); c_alpha_avg->Update();

    // Charge histogram ---------------------------------------------
    TCanvas *c_charge = new TCanvas("c_charge","c_charge",550,550,500,450);
    auto ll = new TLine(); auto rl = new TLine();
    c_charge->cd();
    h_charge->Draw();
    ll->DrawLine(mu-nl_sigma*sigma, 0., mu-nl_sigma*sigma, h_charge->GetMaximum());
    rl->DrawLine(mu+nr_sigma*sigma, 0., mu+nr_sigma*sigma, h_charge->GetMaximum());
    c_charge->Modified(); c_charge->Update();

  }
}
