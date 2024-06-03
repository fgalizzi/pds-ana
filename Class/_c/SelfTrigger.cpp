#ifndef sf_bins
  int sf_bins = 100;
#endif

#ifndef sf_hmin
  double sf_hmin = -2.;
#endif

#ifndef sf_hmax
  double sf_hmax = 7.;
#endif

// Thr is in [N p.e.] unit
#ifndef thr
  double thr = 8./6.55;
#endif

#ifndef esteban
  bool esteban = true;
#endif

std::vector<double> TriggerTime(std::vector<double>& waveform){
  std::vector<double> trgs;
  for(size_t i=0; i<waveform.size(); i++) if (waveform[i] > 0.5){
    trgs.push_back(i);
    if(esteban == true) i += 258;
  }
  return trgs;
}

//*********************************************
void SelfHistos(std::vector<std::vector<double>>& all_wf,
    std::vector<std::vector<double>>& trg_wf, TH1D* h_all, TH1D* h_trg,
    std::vector<double>& int_wf, int I_low, int I_up, double spe_charge, double pedestal,
    double pretrg, double afttrg){
//*********************************************
  int len = all_wf[0].size();
  int count = 0;
  std::vector<int> trg_bool;
  TH1D* hwf = new TH1D("hwf","hwf",len,0,len);
  double spe_norm = 1./spe_charge;
  double t;
  std::vector<double> trgs;
  bool got_ya;
   
  // All wf spectra
  for(auto wf : all_wf){
    for(int i=0; i<len; i++) hwf->SetBinContent(i, wf[i]);
    ComputeIntegral(hwf, int_wf, I_low, I_up);
    hwf->Reset();
  }
  
  
  // Spectra of true positive
  for(auto wf : trg_wf){
    trgs = TriggerTime(wf);
    got_ya = 0;
    for (auto tr : trgs){
      if (tr > pretrg && tr < afttrg) got_ya = 1;
    }
    trg_bool.push_back(got_ya);
  }

  for (size_t i=0; i<int_wf.size(); i++){
    t = (int_wf[i]-pedestal)*spe_norm;
    h_all->Fill(t);
    if (trg_bool[i] == 1) h_trg->Fill(t);
  }
  std::cout << "#Self-Trigger in coincedence with LED " <<  count << std::endl;
 
}



//*** MAIN ************************************
void cla::SelfTrigger(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  vector<vector<double>> trg_wf;
  vector<double> int_wf, n_pe, eff_pe;
  double t;

  // Read and subtract the baseline
  read();
  
  CompleteWF_Binary_Swap(trg_f, trg_wf, n_wf, memorydepth);
  //if(invert == false) SubBaseline(trg_wf, prepulse_ticks);
  //if(invert == true ) SubBaseline_Invert(trg_wf, prepulse_ticks);
  
  //DisplayWFs(trg_wf, 1., 10);

  TH1D* hAll  = new TH1D("hAll" ,"hAll", sf_bins, sf_hmin, sf_hmax);
  TH1D* hTrg  = new TH1D("hTrg" ,"hTrg", sf_bins, sf_hmin, sf_hmax);
  hAll->GetXaxis()->SetTitle("N p.e.");
  hAll->GetYaxis()->SetTitle("Counts");
  hTrg->SetLineColor(kRed);
  
  SelfHistos(wfs, trg_wf, hAll, hTrg, int_wf, int_low, int_up, spe_charge, pedestal, pretrg, afttrg);

  TH1D* hTP = new TH1D("hTP", "", sf_bins, sf_hmin, sf_hmax);
  TH1D* hFP = new TH1D("hFP", "", sf_bins, sf_hmin, sf_hmax);
  
  
  for (int i=1; i<hAll->GetNbinsX(); i++){
    t = hAll->GetBinCenter(i);
    if (t<thr){
      hFP->SetBinContent(i, hTrg->GetBinContent(i));
    }
    else{
      hTP->SetBinContent(i, hTrg->GetBinContent(i));
    }
  }
 
  TEfficiency* eTP = 0;
  TEfficiency* eFP = 0;
  TEfficiency* eTrg= 0;
  
  eTP  = new TEfficiency(*hTP, *hAll); eTP->SetLineColor(kBlue);
  eFP  = new TEfficiency(*hFP, *hAll); eFP->SetLineColor(kBlack);
  eTrg = new TEfficiency(*hTrg,*hAll); eTrg->SetLineColor(kBlack);

  eTP->SetLineWidth(0.); eTP->SetMarkerStyle(20); eTP->SetMarkerColor(kOrange+8);
  eFP->SetLineWidth(0.); eFP->SetMarkerStyle(20); eFP->SetMarkerColor(kAzure+6);

  TCanvas *c_tr = new TCanvas("c_tr","c_tr",0,0,1000,900);
  c_tr->cd();
  hAll->Draw();
  hTrg->Draw("SAME");
  c_tr->Modified();
  c_tr->Update();

  TF1 *f1 = new TF1("f1","1/(1+exp(([0]-x)/[1]))",thr-0.6, thr+0.6);
  f1->SetParameters(thr, 0.2);
  f1->SetParNames("t_{0}", "#tau");
  f1->SetNpx(2000);
  
  eTrg->SetTitle(";N pe; Probability");
  TCanvas *c_eff = new TCanvas("c_eff","c_eff",10,10,1000,900);
  c_eff->cd();
  eTrg->Draw();
  eTrg->Fit(f1, "R");
  c_eff->Modified();
  c_eff->Update();
  
  auto l = new TLine(thr,0.,thr,1.);
  l->SetLineWidth(4);
  l->SetLineStyle(2);
  l->SetLineColor(kGreen-3);
  auto legend = new TLegend(0.7,0.45,0.88,0.55);
  legend->SetBorderSize(0);
  legend->AddEntry(eTP, "True Positive", "lep");
  legend->AddEntry(eFP, "False Positive", "lep");
  legend->AddEntry(l, "Trigger Threshold", "l");
  eTP->SetTitle(";N pe; Probability");

  TCanvas *c_pos = new TCanvas("c_pos","c_pos",20,20,1000,900);
  c_pos->cd();
  eTP->Draw();
  eFP->Draw("SAME");
  l->Draw("SAME");
  legend->Draw("SAME");
  c_pos->Modified();
  c_pos->Update();

  t = f1->Eval(thr);
  std::cout << "\n\nSigmoid t0 - sigma - eff @ thr" << std::endl;
  std::cout << f1->GetParameter(0) << "\t" << f1->GetParameter(1) << "\t" << t << std::endl;  

  if(print== true){
    TString output_name = "Eff_";
    output_name += Form("%.2f", thr);
    TFile* out = new TFile(output_name+".root", "recreate");
    auto eff_dir = out->mkdir("eff");
    eff_dir->cd();
    eTrg->Write();
    out->Close();
  }
}

/*
  for (int i=1; i<hAll->GetNbinsX(); i++){
    n_pe.push_back(hAll->GetBinCenter(i));
    t = hAll->GetBinContent(i);
    if (t==0) eff_pe.push_back(0);
    else {
      t = hTrg->GetBinContent(i)/t;
      eff_pe.push_back(t);
    }
  }
  std::cout << eff_pe.size() << " " << n_pe.size() << std::endl;  
  TGraph *g1 = new TGraph(eff_pe.size(), &n_pe[0], &eff_pe[0]);
  g1->GetXaxis()->SetTitle("N p.e.");
  g1->GetYaxis()->SetTitle("Efficiency");
  */
 
