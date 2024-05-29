

//*********************************************
double TriggerTimee(std::vector<double>& waveform) {
//*********************************************
 
  double* wf = waveform.data();
  auto* g_wf = new TGraph(waveform.size(), wf);
  double xt, y0;

  y0 = *max_element(std::begin(waveform), std::end(waveform));

  if (y0<0.5) return -5;
  xt = g_find_x(g_wf, 0.5, double(PRETRG), double(AFTTRG), 0.001);

  return xt;
}

double TriggerTime(std::vector<double>& waveform){
  double t = -10;
  for(size_t i=0; i<waveform.size(); i++) if (waveform[i] > 0.5){t = i; break;}
  return t;
}

//*** MAIN ************************************
void cla::Jitter(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  vector<vector<double>> trg_wf, all_wf, bulk_wf;
  vector<double> int_wf;
  double t;
  int no_trg_count = 0;

  CompleteWF_Binary_Swap(TRGWF_FILE, trg_wf, n_wf, memorydepth);
  //CompleteWF_Binary_Swap(ALLWF_FILE, all_wf, n_wf, memorydepth);
  

  TH1D* hTrg  = new TH1D("hTrg" ,"hTrg", memorydepth, 0, memorydepth);

  for (size_t i=0; i<trg_wf.size(); i++){
    t = TriggerTime(trg_wf[i]);
    if(t>PRETRG && t<AFTTRG) hTrg->Fill(t);
    else no_trg_count++;
    //if(t>1490 && t<1520) bulk_wf.push_back(all_wf[i]);
  }

  //DisplayWFs(bulk_wf, 1., bulk_wf.size());
  
  std::cout << "No triggers in " << no_trg_count << "/" << n_wf << std::endl;  

  TCanvas *c_tr = new TCanvas("c_tr","c_tr",20,20,1000,900);
  c_tr->cd();
  hTrg->Draw();
  c_tr->Modified();
  c_tr->Update();
  
}
