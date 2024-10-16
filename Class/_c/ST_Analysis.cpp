#include "../classe.hpp"

// Self-Trigger Analysis
// External assumption: spe_charge, pedestal (pay attention to use the same integration
//                      window). To re-normalize the charge spectrum in p.e. units
// Build the trigger-time distribution, auto-set the acceptance window (coincidence 
// LED-selftrigger) according to the Non-LED uniform distribution, build the LED and the
// self-trigger histos, compute the efficiency / effective threshold / ... / with a 
// TEfficiency(LED_histo, self_histo) and print everything in a file.csv  

void cla::ST_Analysis(){
  vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 
  vector<vector<double>> filt_wf;
  vector<double> int_wf;
  vector<double> trgs;
  double t;
  int dark_trg_count = 0;

  feature_value.push_back({"Run", (double)calibration_run});
  feature_value.push_back({"Channel", (double)channel});

  
  CompleteWF_Binary_Swap(trg_f, trg_wf, n_wf, memorydepth);

  TH1D* h_Trg_Time  = new TH1D("h_Trg_Time", Form("%s;%s;%s","h_Trg_Time","Ticks","Counts"), memorydepth, 0, memorydepth);

  for (size_t i=0; i<trg_wf.size(); i++){
    trgs = TriggerTime(trg_wf[i]);
    for (auto tr : trgs){
      if (tr<prepulse_ticks) dark_trg_count++;
      if (tr >= 1) h_Trg_Time->Fill(tr); //Exclude triggers before the opening window
    }
  }
  auto x_h_Trg_Time_max = h_Trg_Time->GetMaximumBin();
  auto y_h_Trg_Time_max = h_Trg_Time->GetBinContent(x_h_Trg_Time_max);
  
  TF1* fc = new TF1("fc", "pol0", 2., prepulse_ticks);
  fc->SetNpx(2000); h_Trg_Time->Fit(fc, "R");
  double avg_trg_counts = fc->GetParameter(0);
  double thr_counts = avg_trg_counts+sqrt(avg_trg_counts);
 
  TGraph g_Trg_Time(h_Trg_Time);

  for (auto i = x_h_Trg_Time_max; i>1; i--){
    pretrg = i;
    if (h_Trg_Time->GetBinContent(i-1)<thr_counts) i=1;
  }
  for (auto i = x_h_Trg_Time_max; i<memorydepth; i++){
    afttrg = i;
    if (h_Trg_Time->GetBinContent(i+1)<thr_counts) i=memorydepth;
  }

  if(pretrg < prepulse_ticks){
    std::cout << "\n\npretrg<prepulse_ticks !!\n\n" << std::endl;
    pretrg = prepulse_ticks+1;
    return;
  }

  TF1* f1 = new TF1("f1", Form("%f+gaus",avg_trg_counts), pretrg, afttrg);
  if (manual == true ) f1 = new TF1("f1",Form("%f+gaus", avg_trg_counts), fit_low, fit_up);
  t = h_Trg_Time->GetBinContent(x_h_Trg_Time_max)-avg_trg_counts; // Gaussian height
  f1->SetParameters(t, h_Trg_Time->GetMaximumBin(), 1.4);
  f1->SetNpx(2000);
  h_Trg_Time->Draw(); f1->Draw("SAME");

  // Set t = to the uniform + gauss/2 height
  t = y_h_Trg_Time_max;
  std::cout << " " << t << std::endl;
  // Compute discrete FWHM
  int fwhm = 0;
  for(int i=0; i<h_Trg_Time->GetEntries(); i++){
    if (h_Trg_Time->GetBinContent(i)>t/2) {
      while (h_Trg_Time->GetBinContent(i)>t/2) {
       fwhm++; i++; 
      }
      continue;
    }
  }

  t = avg_trg_counts+f1->GetParameter(0)*0.5;
  double fwhm_continuous = g_find_x(&g_Trg_Time, t, f1->GetParameter(1), double(afttrg), 1e-3) -
                           g_find_x(&g_Trg_Time, t, double(pretrg), f1->GetParameter(1), 1e-3);
  
  if (display == true) DisplayWFs(trg_wf, 1., 10);
  
  TCanvas *c_trg = new TCanvas("c_tr","c_tr",20,20,1000,900);
  c_trg->cd();
  h_Trg_Time->SetTitle(This_Directory_Name().c_str());
  h_Trg_Time->SetName(This_Directory_Name().c_str());
  h_Trg_Time->Draw();
  h_Trg_Time->Fit(f1, "R");
  c_trg->Modified();c_trg->Update();
 
  t = dark_trg_count/(tick_len*(prepulse_ticks-1)*trg_wf.size())*1.e6;
  double coincidences = dark_trg_count*(afttrg-pretrg) / prepulse_ticks;

  vector<double> t_templ, n_pe, eff_pe, fps, tps;
  double max_eff, thr = -1e6;
  bool stop_search = false;
  size_t nsample = memorydepth;
  TComplex G[nsample];
  std::vector<double> thrs = {1,2,3,4,5}; // N pe

  // Read and subtract the baseline
  read();
  
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
    std::cout << "\n\nHAVE TO MODIFY\n\n" << std::endl;
    self_histos(hAll, hTrg, int_wf);
  }
  if(apply_filter == 0) self_histos(hAll, hTrg, int_wf);

  
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

  TF1* f_sigmoid = new TF1("f_sigmoid","[2]/(1+exp(([0]-x)/[1]))", sf_hmin, sf_hmax);
  f_sigmoid->SetParameters(thr, 0.2, max_eff);
  f_sigmoid->SetParNames("t_{0}", "#tau", "#epsilon_{MAX}");
  f_sigmoid->SetParLimits(2, 0, 1);
  f_sigmoid->SetNpx(2000);

  eTrg->Fit(f_sigmoid, "R");
  thr = f_sigmoid->GetParameter(0);


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

  feature_value.push_back({"Accuracy thr [pe]", hAcc->GetBinCenter(hAcc->GetMaximumBin())});
  feature_value.push_back({"Acc_thr err", hAcc->GetBinWidth(2)});
  feature_value.push_back({"Sigmoid thr [pe]", thr});
  feature_value.push_back({"Sig_thr err", f_sigmoid->GetParError(0)});
  feature_value.push_back({"Tau [pe]", f_sigmoid->GetParameter(1)});
  feature_value.push_back({"Tau err", f_sigmoid->GetParError(1)});
  feature_value.push_back({"Eff_max", f_sigmoid->GetParameter(2)});
  feature_value.push_back({"Eff_max err", f_sigmoid->GetParError(2)});

  for (size_t i=0; i<thrs.size(); i++) {
    feature_value.push_back({Form("fp%zu",i+1), fps[i]});
    feature_value.push_back({Form("tp%zu",i+1), tps[i]});
  }
  
  feature_value.push_back({"pretrg [ticks]", pretrg});
  feature_value.push_back({"afttrg [ticks]", afttrg});
  feature_value.push_back({"N trg", hTrg->GetEntries()});

  t = dark_trg_count/(tick_len*(prepulse_ticks-1)*trg_wf.size())*1.e6;
  coincidences = dark_trg_count*(afttrg-pretrg) / prepulse_ticks;
  
  feature_value.push_back({"T0 trg [ticks]", f1->GetParameter(1)});
  feature_value.push_back({"Jitter [ticks]", f1->GetParameter(2)});
  feature_value.push_back({"Jitter err", f1->GetParError(2)});
  feature_value.push_back({"FWHM_d", fwhm});
  feature_value.push_back({"FWHM_c", fwhm_continuous});
  feature_value.push_back({"Non-LED rate [Hz]", t});
  feature_value.push_back({"Coincidences", coincidences});



  if(print==true) print_vec_pair_csv(Form("ST_results_run_%zu.csv",calibration_run), feature_value);


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


if(print== true){
    TString output_name = "../Self_FOM.root";
    TFile* out = new TFile(output_name, "update");
    auto jit_dir = out->mkdir("jitter");
    out->cd("jitter");
    h_Trg_Time->Write(("jit_"+wf_file).c_str());
    out->Close();
  }
}
