#include "../classe.hpp"

// ****************************************************************
// Description
// ****************************************************************

//*** MAIN ************************************
//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Jitter(int st_channel){
//*********************************************
  vector<double> int_wf;
  vector<double> trgs;
  double t;
  int dark_trg_count = 0;

  if(pretrg < prepulse_ticks){
    std::cout << "\n\npretrg<prepulse_ticks !!\n\n" << std::endl;
    return;
  }
  // CompleteWF_Binary_Swap(trg_f, trg_wf, n_wf, memorydepth);
  StructuredWaveformSetReader(wf_file, trg_wf, st_channel, n_wf);

  TH1D* hTrg  = new TH1D("hTrg", Form("%s;%s;%s","hTrg","Ticks","Counts"), memorydepth, 0, memorydepth);

  for (size_t i=0; i<trg_wf.size(); i++){
    trgs = TriggerTime(trg_wf[i]);
    for (auto tr : trgs){
      if (tr<prepulse_ticks) dark_trg_count++;
      if (tr >= 1) hTrg->Fill(tr); //Exclude triggers before the opening window
    }
  }
  auto x_hTrg_max = hTrg->GetMaximumBin();
  auto y_hTrg_max = hTrg->GetBinContent(x_hTrg_max);
  
  TF1* fc = new TF1("fc", "pol0", 2., prepulse_ticks);
  fc->SetNpx(2000); hTrg->Fit(fc, "R");
  double avg_trg_counts = fc->GetParameter(0);
  double thr_counts = avg_trg_counts+sqrt(avg_trg_counts);
 
  TGraph g_hTrg(hTrg);

  std::cout << "xmax" << x_hTrg_max << std::endl;
  for (auto i = x_hTrg_max; i>1; i--){
    pretrg = i;
    if (hTrg->GetBinContent(i-1)<thr_counts) i=1;
  }
  for (auto i = x_hTrg_max; i<memorydepth; i++){
    afttrg = i;
    if (hTrg->GetBinContent(i+1)<thr_counts) i=memorydepth;
  }


  TF1* f1 = new TF1("f1", Form("%f+gaus",fc->GetParameter(0)), pretrg, afttrg);
  if (manual == true ) f1 = new TF1("f1",Form("%f+gaus", fc->GetParameter(0)), fit_low, fit_up);
  t = hTrg->GetBinContent(x_hTrg_max)-fc->GetParameter(0); // Gaussian height
  f1->SetParameters(t, hTrg->GetMaximumBin(), 1.4);
  f1->SetNpx(2000);
  hTrg->Draw(); f1->Draw("SAME");

  // Set t = to the uniform + gauss/2 height
  t = hTrg->GetBinContent(hTrg->GetMaximumBin());
  std::cout << " " << t << std::endl;
  // Compute discrete FWHM
  int fwhm = 0;
  for(int i=0; i<hTrg->GetEntries(); i++){
    if (hTrg->GetBinContent(i)>t/2) {
      while (hTrg->GetBinContent(i)>t/2) {
       fwhm++; i++; 
      }
      continue;
    }
  }

  t = fc->GetParameter(0)+f1->GetParameter(0)*0.5;
  double fhwm_continuous = g_find_x(&g_hTrg, t, f1->GetParameter(1), double(afttrg), 1e-3) -
                           g_find_x(&g_hTrg, t, double(pretrg), f1->GetParameter(1), 1e-3);
  std::cout << "\nFWHM " << fwhm << " " << fhwm_continuous << std::endl;
  if (display == true) DisplayWFs(trg_wf, 1., 10);
  
  TCanvas *c_trg = new TCanvas("c_tr","c_tr",20,20,1000,900);
  c_trg->cd();
  hTrg->SetTitle(This_Directory_Name().c_str());
  hTrg->SetName(This_Directory_Name().c_str());
  hTrg->Draw();
  hTrg->Fit(f1, "R");
  c_trg->Modified();c_trg->Update();
 
  t = dark_trg_count/(tick_len*(prepulse_ticks-1)*trg_wf.size())*1.e6;
  double coincidences = double(dark_trg_count*(afttrg-pretrg) / prepulse_ticks);
  std::cout << "\n\nTrg Abs Time - Jitter - Err - FWHM - Dark trigger rate [Hz] - Coincidences"<< std::endl;
  std::cout << f1->GetParameter(1) << "\t" << f1->GetParameter(2) << "\t" 
    << f1->GetParError(2) << "\t" << fwhm << "\t" << t << "\t" << coincidences << "\n\n" << std::endl;  

if(print== true){
    TString output_name = "../Self_FOM.root";
    TFile* out = new TFile(output_name, "update");
    auto jit_dir = out->mkdir("jitter");
    out->cd("jitter");
    hTrg->Write(("jit_"+wf_file).c_str());
    out->Close();
  }

}
