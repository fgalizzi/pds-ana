size_t calibration_run = 27563;
int pspe_low = 190;
int pspe_up  = 490;

vector<size_t> read_ch_map(){
  ifstream file_map("/Users/federico/CERN/PDHD/waffles/scripts/cpp_utils/functions/channelmap.txt");
  string line;
  stringstream ssmap;
  Short_t dpch, ch;
  vector<size_t> channels;

  if (file_map.is_open()){
    while (getline(file_map, line)) {
      ssmap.clear();
      ssmap.str(line);

      while (ssmap >> dpch >> ch) {
        channels.push_back(int(dpch));
      }
    }
  }

  return channels;
}
void cla::ProtoDUNE_Calibration(){
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 
  
  //vector<size_t> channels = {11147, 11147, 11147};//, 11145, 11147};
  vector<size_t> channels = read_ch_map();
  vector<vector<double>> sel_wf;
  vector<double> int_wf, spe_avg;

  TFile hf(Form("CalRun2_%zu.root", calibration_run), "recreate");
  hf.mkdir("chargehistos");
  hf.mkdir("spe_wfs");
  
  TTree results("calib_results", "calib_results");
  TH1D h1d; TF1 fun1;
  results.Branch("chg_histo", &h1d);
  results.Branch("fit_fun", &fun1);
   
  vector<TH1D*> h_charge;
  vector<TGraph*> g_spe;
  // Run - Channel - Good/Bad - SNR - gain - spe ampl - 
  vector<tuple<size_t, size_t, bool, double, double, double>> res_tuple;

  for(size_t ch_index=0; ch_index<channels.size(); ch_index++) {
    
    channel = channels[ch_index];
    if (channel < 11100 || channel > 11400) continue;
    // Read the wfs for this channel and subtract the baseline
    read();

    SelCalib_WF(wfs, sel_wf, prepulse_ticks, sat_low, sat_up, bsl);  
    
    TH1D* hI = new TH1D();
    hI = BuildRawChargeHisto(sel_wf, int_wf, int_low, int_up, nbins);

    TH1D* hFind = new TH1D(*hI);
    // Go for peaks: create an instance of TSpectrum
    hFind->Rebin(8);
    TSpectrum *s = new TSpectrum(5);
    s->Background(hFind, 50, "same");
    int npeaks = s->Search(hFind,1,"",0.2);
    Double_t *xpeaks;
    xpeaks = s->GetPositionX();

    std::vector<double> xp;
    for(int i=0; i<npeaks; i++) xp.push_back(xpeaks[i]);
    std::sort(xp.begin(), xp.end());
    if (xp.size()<2){
      cout << "\nWARNING: using spe_low - spe_up\n" << endl; 
      xp = {0, (pspe_low+pspe_up)/2.};
    }
    

    // Parameters for the fit function
    double par[20] = {0}; //Norm_Consts, sigma_0 , sigma_cell, G
    
    // Fit on the first peak to estimate sigma_0
    TF1* f0 = new TF1("f0", "gaus", -xp[1]*0.3, xp[1]*0.3);
    f0->SetParameters(hI->GetBinContent(hI->GetMaximumBin()), xp[0], xp[1]*0.5);
    hI->Fit("f0", "R");
    
    // Initial guess using TSpectrum
    par[0] = xp[0];               //peak 0
    par[1] = xp[1]-xp[0];         //peak 1
    par[2] = f0->GetParameter(2); //sigma_0
    par[3] = par[2]*0.1;          //sigma_cell
   
    if(hI->GetBinContent(hI->GetMaximumBin())<70) hI->Rebin(2);
    if(npeaks<nmaxpeaks) npeaks=nmaxpeaks;
    for(int i = 0 ; i < npeaks ; i++){
      par[i+4] = hI->GetBinContent(hI->GetMaximumBin())*0.8;//peaky[i];
    }

    // Gaussian sumatory
    fgaus = new TF1("fgaus", fRandomName, -par[2]*2., par[1]*(npeaks-1), npeaks+4);
    fRandomName_set(fgaus);
    fgaus->SetParameters(par);
    fgaus->SetParLimits(0, -par[1]*0.2, par[1]*0.2);
    fgaus->SetParLimits(1, par[1]*0.4, par[1]*1.9);
    fgaus->SetParLimits(2, par[2]*0.3, par[2]*2.);
    fgaus->SetParLimits(3, 0., par[2]*2.);
    
    for(int i = 0 ; i < npeaks ; i++) fgaus->SetParLimits(i+4,0,2700);

    hI->Fit("fgaus", "R");
    hI->Draw();//gPad->Update();gPad->WaitPrimitive();
    
    bool goodness = 0;
    // GOODNESS of fit !!
    LoadFitParameters(fgaus);
    Avg_Sel_WF(sel_wf, spe_avg, int_wf, spe_low, spe_up);
    spe_ampl = *max_element(spe_avg.begin()+int_low, spe_avg.begin()+int_up);
   
    h1d = *hI; fun1 = *fgaus;
    results.Fill();

    hI->SetTitle(Form("Channel_%lu", channels[ch_index]));
    h_charge.push_back(hI);
    TGraph* gr = new TGraph(spe_avg.size(), data(spe_avg));
    gr->SetTitle(Form("Channel_%lu", channels[ch_index]));
    g_spe.push_back(gr);
    
    double SNR = fgaus->GetParameter(1)/fgaus->GetParameter(2); 
    double gain = fgaus->GetParameter(1); 
    if (SNR > 5.) goodness = 1;
    res_tuple.push_back(make_tuple(calibration_run, channel, goodness, SNR, gain, spe_ampl));

    int_wf.erase(int_wf.begin(), int_wf.end());
    sel_wf.erase(sel_wf.begin(), sel_wf.end());
  } 

  for(auto tuple : res_tuple) print_tuple(tuple);
 
  // Save in file.root
  // results.Write();
  hf.cd("chargehistos"); for(auto h : h_charge) h->Write();
  hf.cd("spe_wfs");      for(auto g : g_spe)    g->Write();
  
  hf.Close();
}
