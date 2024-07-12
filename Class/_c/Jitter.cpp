
//*** MAIN ************************************
void cla::Jitter(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  vector<vector<double>> trg_wf, all_wf, bulk_wf;
  vector<double> int_wf;
  vector<double> trgs;
  double t;
  int dark_trg_count = 0;

  if(pretrg < prepulse_ticks){
    std::cout << "\n\npretrg<prepulse_ticks !!\n\n" << std::endl;
    return;
  }
  CompleteWF_Binary_Swap(trg_f, trg_wf, n_wf, memorydepth);

  TH1D* hTrg  = new TH1D("hTrg" ,"hTrg:Ticks:Counts", memorydepth, 0, memorydepth);
  hTrg->GetXaxis()->SetTitle("Time [ticks]");
  hTrg->GetYaxis()->SetTitle("Counts");

  for (size_t i=0; i<trg_wf.size(); i++){
    trgs = TriggerTime(trg_wf[i]);
    for (auto tr : trgs){
      if (tr<prepulse_ticks) dark_trg_count++;
      if (tr >= 1) hTrg->Fill(tr); //Exclude triggers before the opening window
    }
  }

  TF1* fc = new TF1("fc", "pol0", 2., prepulse_ticks);
  fc->SetNpx(2000); hTrg->Fit(fc, "R");
  
  TF1* f1 = new TF1("f1", Form("%f+gaus",fc->GetParameter(0)), pretrg, afttrg);
  if (manual == true ) f1 = new TF1("f1",Form("%f+gaus", fc->GetParameter(0)), fit_low, fit_up);
  t = hTrg->GetBinContent(hTrg->GetMaximumBin())-fc->GetParameter(0); // Gaussian height
  f1->SetParameters(t, hTrg->GetMaximumBin(), 1.4);
  f1->SetNpx(2000);
  hTrg->Draw(); f1->Draw("SAME");

  // Set t = to the uniform + gauss/2 height
  t = (hTrg->GetBinContent(hTrg->GetMaximumBin())+fc->GetParameter(0))*0.5; 
  int fwhm = 0;
  for(int i=0; i<hTrg->GetEntries(); i++){
    if (hTrg->GetBinContent(i)>t/2) {
      while (hTrg->GetBinContent(i)>t/2) {
       fwhm++; i++; 
      }
      continue;
    }
  }
 
  if (display == true) DisplayWFs(trg_wf, 1., 10);
  
  TCanvas *c_trg = new TCanvas("c_tr","c_tr",20,20,1000,900);
  c_trg->cd();
  hTrg->SetTitle(This_Directory_Name().c_str());
  hTrg->SetName(This_Directory_Name().c_str());
  hTrg->Draw();
  hTrg->Fit(f1, "R");
  c_trg->Modified();
  c_trg->Update();
 
  t = dark_trg_count/(tick_len*(memorydepth-prepulse_ticks-1)*trg_wf.size())*1.e6;
  double coincidences = dark_trg_count*(afttrg-pretrg) / prepulse_ticks;
  std::cout << "\n\nTrg Abs Time - Jitter - Err - FWHM - Dark trigger rate [Hz] - Coincidences"<< std::endl;
  std::cout << f1->GetParameter(1) << "\t" << f1->GetParameter(2) << "\t" 
    << f1->GetParError(2) << "\t" << fwhm << "\t" << t << "\t" << coincidences << "\n\n" << std::endl;  

if(print== true){
    TString output_name = "../Self_FOM.root";
    TFile* out = new TFile(output_name, "update");
    //auto jit_dir = out->mkdir("jitter");
    out->cd("jitter");
    hTrg->Write();
    out->Close();
  }

}
