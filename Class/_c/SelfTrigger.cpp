#ifndef sf_bins
  int sf_bins = 100;
#endif

#ifndef sf_hmin
  double sf_hmin = -2.;
#endif

#ifndef sf_hmax
  double sf_hmax = 7.;
#endif

std::vector<double> TriggerTime(std::vector<double>& waveform){
  std::vector<double> trgs;
  for(size_t i=0; i<waveform.size(); i++) if (waveform[i] > 0.5){
    trgs.push_back(i);
    while(waveform[i]>0.5) i++;
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
  int got_ya;
   
  // All wf spectra
  for(auto wf : all_wf){
    for(int i=0; i<len; i++) hwf->SetBinContent(i, wf[i]);
    ComputeIntegral(hwf, int_wf, I_low, I_up);
    hwf->Reset();
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
  std::cout << "#Self-Trigger in coincedence with LED " <<  count << std::endl;
 
}



//*** MAIN ************************************
void cla::SelfTrigger(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  vector<vector<double>> trg_wf;
  vector<double> int_wf, n_pe, eff_pe, fps, tps;
  double t, thr = -1e6;
  bool stop_search = false;
  std::vector<double> thrs = {1,2,3,4,5}; // N pe

  // Read and subtract the baseline
  read();
  
  CompleteWF_Binary_Swap(trg_f, trg_wf, n_wf, memorydepth);
  //if(invert == false) SubBaseline(trg_wf, prepulse_ticks);
  //if(invert == true ) SubBaseline_Invert(trg_wf, prepulse_ticks);
  
  //DisplayWFs(trg_wf, 1., 10);

  TH1D* hAll  = new TH1D("hAll" ,"hAll", sf_bins, sf_hmin, sf_hmax);
  TH1D* hTrg  = new TH1D("hTrg" ,"hTrg", sf_bins, sf_hmin, sf_hmax);
  TH1D* hAcc  = new TH1D("hAcc" ,"hAcc", sf_bins, sf_hmin, sf_hmax);
  hAll->GetXaxis()->SetTitle("N p.e."); hAcc->GetXaxis()->SetTitle("N p.e.");
  hAll->GetYaxis()->SetTitle("Counts"); hAcc->GetYaxis()->SetTitle("Counts");
  hTrg->SetLineColor(kRed);
  
  SelfHistos(wfs, trg_wf, hAll, hTrg, int_wf, int_low, int_up, spe_charge, pedestal, pretrg, afttrg);

  
  // Efficiency as Trg/All + fit for effective threshold
  TEfficiency* eTrg= 0;
  eTrg = new TEfficiency(*hTrg,*hAll); eTrg->SetLineColor(kBlack);
  
  for (int i=1; i<hAll->GetNbinsX(); i++){
    hAcc->SetBinContent(i,hAll->Integral(1,i)-hTrg->Integral(1,i)+hTrg->Integral(i,hAll->GetNbinsX()));
    if (eTrg->GetEfficiency(i)>0.45 && stop_search==0){
      thr = hAll->GetBinCenter(i);
      stop_search = 1;
    }
  }

  TF1* f1 = new TF1("f1","[2]/(1+exp(([0]-x)/[1]))",thr-3, thr+3);
  f1->SetParameters(thr, 0.2, 0.8);
  if(manual == true) {
    f1 = new TF1("f1","[2]/(1+exp(([0]-x)/[1]))",fit_low, fit_up);
    f1->SetParameters(t_0, 0.2, 0.8);
  }
  f1->SetParNames("t_{0}", "#tau", "#epsilon_{MAX}");
  f1->SetParLimits(2, 0, 1);
  f1->SetNpx(2000);

  eTrg->Fit(f1, "R");
  thr = f1->GetParameter(0);


  // T = True ; P = Positive
  TH1D* hTP = new TH1D("hTP", "", sf_bins, sf_hmin, sf_hmax);
  TH1D* hFP = new TH1D("hFP", "", sf_bins, sf_hmin, sf_hmax);
 
  for (auto thr_pe : thrs){
    hTP->Reset(); hFP->Reset();
    for (int i=1; i<hAll->GetNbinsX(); i++){
      t = hAll->GetBinCenter(i);
      if (t<thr_pe){
        hFP->SetBinContent(i, hTrg->GetBinContent(i));
      }
      else{
        hTP->SetBinContent(i, hTrg->GetBinContent(i));
      }
    }
    fps.push_back( hTrg->Integral(1, hTrg->FindBin(thr_pe) ) /
        hAll->Integral(1, hAll->FindBin(thr_pe) ) );
    tps.push_back( hTrg->Integral(hTrg->FindBin(thr_pe), hAll->GetNbinsX()) /
        hAll->Integral(hAll->FindBin(thr_pe), hAll->GetNbinsX()) );
  }
  
 
  TEfficiency* eTP = 0;
  TEfficiency* eFP = 0;
  
  eTP  = new TEfficiency(*hTP, *hAll); eTP->SetLineColor(kBlue);
  eFP  = new TEfficiency(*hFP, *hAll); eFP->SetLineColor(kBlack);

  eTP->SetLineWidth(0.); eTP->SetMarkerStyle(20); eTP->SetMarkerColor(kOrange+8);
  eFP->SetLineWidth(0.); eFP->SetMarkerStyle(20); eFP->SetMarkerColor(kAzure+6);

  TCanvas *c_tr = new TCanvas("c_tr","c_tr",0,0,1000,900);
  c_tr->cd();
  hAll->Draw();
  hTrg->Draw("SAME");
  c_tr->Modified();
  c_tr->Update();
  
  
  eTrg->SetTitle(";N pe; Probability");
  TCanvas *c_eff = new TCanvas("c_eff","c_eff",10,10,1000,900);
  c_eff->cd();
  eTrg->Draw();
  c_eff->Modified();
  c_eff->Update();
  
  /*auto l = new TLine(thr,0.,thr,1.);
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
  */
  /*TCanvas *c_acc = new TCanvas("c_acc","c_acc",30,30,1000,900);
  c_acc->cd();
  hAcc->Draw();
  c_acc->Modified();
  c_acc->Update();*/

  std::cout << "\n\nAccucarcy thr - err - Sigmoid thr - err - Tau - Err - Eff max - Err - Fps and Tps" << std::endl;
  std::cout << hAcc->GetBinCenter(hAcc->GetMaximumBin()) << "\t" << hAcc->GetBinWidth(2) << "\t"
    << thr << "\t" << f1->GetParError(0) << "\t"
    << f1->GetParameter(1) << "\t" << f1->GetParError(1) << "\t"
    << f1->GetParameter(2) << "\t" << f1->GetParError(2) << "\t"; 
  
  for (size_t i=0; i<thrs.size(); i++) std::cout << fps[i] << "\t" << tps[i] << "\t";

  std::cout<< "\n\n" << std::endl;

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
