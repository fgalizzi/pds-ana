const int SHIIFT = 20;
std::string ALLWF_FILE = "channel_1.dat";
std::string TRGWF_FILE = "channel_2.dat";
int PRETRG = 1430;
int AFTTRG = 1470;

// With the entire set of WFs (all_wf) it  build the calibration histogram integrating [I_low;I_up]
//*********************************************
void SelfHistos(std::vector<std::vector<double>>& all_wf,
    std::vector<std::vector<double>>& trg_wf, TH1D* h_all, TH1D* h_trg,
    std::vector<double>& int_wf, int I_low, int I_up, double hmin, double hmax, int nbins){
//*********************************************
  int len = all_wf[0].size();
  int count = 0;
  std::vector<int> trg_bool;
  TH1D* hwf = new TH1D("hwf","hwf",len,0,len);
   
  // All wf spectra
  for(auto wf : all_wf){
    for(int i=0; i<len; i++) hwf->SetBinContent(i, wf[i]);
    ComputeIntegral(hwf, int_wf, I_low, I_up);
    hwf->Reset();
  }
  
  
  // Spectra of true positive
  for(auto wf : trg_wf){
   if(wf[PRETRG]>-0.1 && wf[AFTTRG]<-0.5)
   {trg_bool.push_back(1);
     count += 1; 
   }
   else trg_bool.push_back(0);
  }

  for (size_t i=0; i<int_wf.size(); i++){
    h_all->Fill(int_wf[i]);
    if (trg_bool[i] == 1) h_trg->Fill(int_wf[i]);
  }
  std::cout << count << std::endl;
 
}

//*** MAIN ************************************
void cla::SelfTrigger(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  vector<vector<double>> all_wf, trg_wf;
  vector<double> int_wf;

  // Read and subtract the baseline
  CompleteWF_Binary_Swap(ALLWF_FILE, all_wf, n_wf, memorydepth);
  CompleteWF_Binary_Swap(TRGWF_FILE, trg_wf, n_wf, memorydepth);
  
  if(invert == false) SubBaseline(all_wf, prepulse_ticks);
  if(invert == true ) SubBaseline_Invert(all_wf, prepulse_ticks);
  if(invert == false) SubBaseline(trg_wf, prepulse_ticks);
  if(invert == true ) SubBaseline_Invert(trg_wf, prepulse_ticks);
  
  //DisplayWFs(trg_wf, 1., 10);

  TH1D* hAll  = new TH1D("hAll" ,"hAll", nbins, hmin, hmax);
  TH1D* hTrg  = new TH1D("hTrg" ,"hTrg", nbins, hmin, hmax);
  hTrg->SetLineColor(kRed);
  SelfHistos(all_wf, trg_wf, hAll, hTrg, int_wf, int_low, int_up, hmin, hmax, nbins);

   
  TCanvas *c_tr = new TCanvas("c_tr","c_tr",20,20,1000,900);
  c_tr->cd();
  hAll->Draw();
  hTrg->Draw("SAME");
  c_tr->Modified();
  c_tr->Update();
  
}
