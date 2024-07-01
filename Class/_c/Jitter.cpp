
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

  CompleteWF_Binary_Swap(trg_f, trg_wf, n_wf, memorydepth);
  //CompleteWF_Binary_Swap(ALLWF_FILE, all_wf, n_wf, memorydepth);
  

  TH1D* hTrg  = new TH1D("hTrg" ,"hTrg:Ticks:Counts", memorydepth, 0, memorydepth);
  hTrg->GetXaxis()->SetTitle("Time [ticks]");
  hTrg->GetYaxis()->SetTitle("Counts");

  for (size_t i=0; i<trg_wf.size(); i++){
    trgs = TriggerTime(trg_wf[i]);
    for (auto tr : trgs){
      if (tr<pretrg || tr>afttrg) dark_trg_count++;
      if (tr >= 1) hTrg->Fill(tr); //Exclude triggers before the opening window
    }
  }

  TF1* f1 = new TF1("f1","gaus", pretrg, afttrg);
  if (manual == true ) f1 = new TF1("f1","gaus", fit_low, fit_up);
  f1->SetParameters(700, (pretrg+afttrg)*0.5, 1.4);
  f1->SetNpx(2000);
  
  if (display == true) DisplayWFs(trg_wf, 1., 10);
  
  TCanvas *c_trg = new TCanvas("c_tr","c_tr",20,20,1000,900);
  c_trg->cd();
  hTrg->Draw();
  hTrg->Fit(f1, "R");
  c_trg->Modified();
  c_trg->Update();
 
  t = dark_trg_count/(tick_len*(memorydepth-afttrg+pretrg-1)*trg_wf.size())*1.e6;

  std::cout << "\n\nTrg Abs Time - Jitter - Err - Dark trigger rate [Hz]"<< std::endl;
  std::cout << f1->GetParameter(1) << "\t" << f1->GetParameter(2) << "\t" 
    << f1->GetParError(2) << "\t" << t << "\n\n" << std::endl;  
}
