  
double jit_low = 1473;
double jit_up = 1482;


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
  int no_trg_count = 0;
  int dark_trg_count = 0;

  CompleteWF_Binary_Swap(trg_f, trg_wf, n_wf, memorydepth);
  //CompleteWF_Binary_Swap(ALLWF_FILE, all_wf, n_wf, memorydepth);
  

  TH1D* hTrg  = new TH1D("hTrg" ,"hTrg", memorydepth, 0, memorydepth);

  for (size_t i=0; i<trg_wf.size(); i++){
    trgs = TriggerTime(trg_wf[i]);
    if (trgs.size() == 0) no_trg_count++;
    else {
      for (auto tr : trgs){
        if (tr<pretrg) dark_trg_count++;
        hTrg->Fill(tr);
      }
    }
    if (t<pretrg && t>-1) dark_trg_count++;
    //if(t>1490 && t<1520) bulk_wf.push_back(all_wf[i]);
  }

  TF1 *f1 = new TF1("f1","gaus", pretrg-3, afttrg+3);
  f1->SetParameters(700, (pretrg+afttrg)*0.5, 1.4);
  f1->SetNpx(2000);
  

  //DisplayWFs(trg_wf, 1., 30);
  


  TCanvas *c_tr = new TCanvas("c_tr","c_tr",20,20,1000,900);
  c_tr->cd();
  hTrg->Draw();
  hTrg->Fit(f1, "R");
  c_tr->Modified();
  c_tr->Update();
 
  t = dark_trg_count/(tick_len*pretrg*trg_wf.size())*1.e6;
  std::cout << "\n\nNo triggers in " << no_trg_count << "/" << n_wf << std::endl;  
  std::cout << "Dark trigger count " << dark_trg_count << std::endl;  
  std::cout << "Dark trigger rate " << t << "\n\n"<< std::endl;

  std::cout << "Jitter - Dark trigger rate - No Trg"<< std::endl;
  std::cout << f1->GetParameter(2) << "\t" << t << "\t" << no_trg_count << std::endl;  
}
