#include "../classe.hpp"

// ****************************************************************
// Description
// ****************************************************************

#ifndef sf_bins
  int sf_bins = 100;
#endif

#ifndef sf_hmin
  double sf_hmin = -2.;
#endif

#ifndef sf_hmax
  double sf_hmax = 7.;
#endif


//*********************************************
void SelfHistos_(std::vector<std::vector<double>>& all_wf,
    std::vector<std::vector<double>>& trg_wf, TH1D* h_all, TH1D* h_trg,
    std::vector<double>& int_wf, int I_low, int I_up, double spe_charge, double pedestal,
    double pretrg, double afttrg, bool apply_filter, int nbins){
//*********************************************
  int len = all_wf[0].size();
  int count = 0;
  std::vector<int> trg_bool;
  double spe_norm = 1./spe_charge;
  double t;
  std::vector<double> trgs;
  int got_ya;
   
  TH1D* hI = new TH1D();
  
  if(apply_filter == 0){
    hI = BuildRawChargeHisto(all_wf, int_wf, I_low, I_up, nbins);
  } 
  else{
    hI = BuildRawChargeHisto(all_wf, int_wf, 0, 1, nbins);
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
//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Self_Trigger(){
//*********************************************
  vector<vector<double>> trg_wf, filt_wf;
  vector<double> int_wf, t_templ, n_pe, eff_pe, fps, tps;
  double t, max_eff, thr = -1e6;
  bool stop_search = false;
  size_t nsample = memorydepth;
  TComplex G[nsample];
  std::vector<double> thrs = {1,2,3,4,5}; // N pe

  // Read and subtract the baseline
  read();
  
  CompleteWF_Binary_Swap(trg_f, trg_wf, n_wf, memorydepth);
  if (display == true) DisplayWFs (wfs, trg_wf, tick_len, 10);

  TH1D* hAll  = new TH1D("hAll" ,"hAll", sf_bins, sf_hmin, sf_hmax);
  TH1D* hTrg  = new TH1D("hTrg" ,"hTrg", sf_bins, sf_hmin, sf_hmax);
  TH1D* hAcc  = new TH1D("hAcc" ,"hAcc", sf_bins, sf_hmin, sf_hmax);
  hAll->GetXaxis()->SetTitle("N p.e."); hAcc->GetXaxis()->SetTitle("N p.e.");
  hAll->GetYaxis()->SetTitle("Counts"); hAcc->GetYaxis()->SetTitle("Counts");
  hTrg->SetLineColor(kRed);
 

  if(apply_filter == 1){
    CompleteWF_Binary(templ_f, t_templ, memorydepth); // t_templ = time domain template
    Build_Matched_Filter(&G[0], t_templ);
    FilterAllWF(wfs, filt_wf, G);
    SelfHistos_(filt_wf, trg_wf, hAll, hTrg, int_wf, int_low, int_up, spe_charge,
             pedestal, pretrg, afttrg, apply_filter, sf_bins);
  }
  if(apply_filter == 0)SelfHistos_(wfs, trg_wf, hAll, hTrg, int_wf, int_low, int_up, spe_charge,
             pedestal, pretrg, afttrg, apply_filter, sf_bins);

  
  // Efficiency as Trg/All + fit for effective threshold
  TEfficiency* eTrg= 0;
  eTrg = new TEfficiency(*hTrg,*hAll); eTrg->SetLineColor(kBlack);
  max_eff = 0;
  for (int i=1; i<hAll->GetNbinsX(); i++){
    if(max_eff < eTrg->GetEfficiency(i)) max_eff = eTrg->GetEfficiency(i);
  }
  
  for (int i=1; i<hAll->GetNbinsX(); i++){
    hAcc->SetBinContent(i,hAll->Integral(1,i)-hTrg->Integral(1,i)+hTrg->Integral(i,hAll->GetNbinsX()));
    if (eTrg->GetEfficiency(i)>max_eff*0.5 && stop_search==0){
      thr = hAll->GetBinCenter(i);
      stop_search = 1;
    }
  }

  TF1* f1 = new TF1("f1","[2]/(1+exp(([0]-x)/[1]))", sf_hmin, sf_hmax);
  f1->SetParameters(thr, 0.2, max_eff);
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
  eTrg->SetTitle(This_Directory_Name().c_str());
  eTrg->SetName(This_Directory_Name().c_str());
  eTrg->Draw();
  c_eff->Modified();
  c_eff->Update();

  std::cout << "\n\nAccucarcy thr - err - Sigmoid thr - err - Tau - Err - Eff max - Err - Fps and Tps - pretrg - afttrg - n_trg" << std::endl;
  std::cout << hAcc->GetBinCenter(hAcc->GetMaximumBin()) << "\t" << hAcc->GetBinWidth(2) << "\t"
    << thr << "\t" << f1->GetParError(0) << "\t"
    << f1->GetParameter(1) << "\t" << f1->GetParError(1) << "\t"
    << f1->GetParameter(2) << "\t" << f1->GetParError(2) << "\t"; 
  
  for (size_t i=0; i<thrs.size(); i++) std::cout << fps[i] << "\t" << tps[i] << "\t";

  std::cout<< pretrg << "\t" << afttrg << "\t" << hTrg->GetEntries() <<  "\n\n" << std::endl;

  if(print== true){
    TString output_name = "../Self_FOM.root";
    TFile* out = new TFile(output_name, "update");
    auto eff_dir = out->mkdir("eff");
    out->cd("eff");
    eTrg->Write(("eff_"+wf_file).c_str());
    auto his_dir = out->mkdir("his");
    out->cd("his");
    hAll->Write(("hAll_"+wf_file).c_str());hTrg->Write(("hTrg_"+wf_file).c_str());
    out->Close();
  }
}
