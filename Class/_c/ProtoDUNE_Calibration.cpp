//Read a calibration run file, build a charge histogram for each channel
//fit it, estimate the spe amplitude, store the histogram and the average spe
//waveform in a root file, print the results at termina//Read a calibration run file, build a charge histogram for each channel
//fit it, estimate the spe amplitude, store the histogram and the average spe
//waveform in a root file, print the results at terminal

size_t calibration_run = 29760; //Needed create the root file with the results
size_t channel_low = 11100;     //Lower channel to look at (included)
size_t channel_up  = 11400;     //Upper " "

//Output file name, then it adds "calibration_run.root":w
string outfile_name = "CalRun_prep124_bsl18_intlow129_up152_";
bool avoid_auto_peak = true;
int pspe_low = 110; //Lower limit for spe integral (like spe_low)
int pspe_up  = 270;//Upper " " Remember: it depends on the integration window,
                    //the overvoltage and the gain. You can also enable the peak finding
                    //and don't use these


void cla::ProtoDUNE_Calibration(){
  //Choose the channels to read
  //vector<size_t> channels = {11147, 11147, 11147};//, 11145, 11147};
  vector<size_t> channels = read_pdhd_ch_map();
  vector<vector<double>> sel_wf;
  vector<double> int_wf, spe_avg;

  outfile_name = outfile_name+calibration_run+".root";
  TFile hf(outfile_name.c_str(), "recreate");
  hf.mkdir("chargehistos");
  hf.mkdir("spe_wfs");
  
  TTree results("calib_results", "calib_results");
  TH1D h1d; TF1 fun1;
  results.Branch("chg_histo", &h1d);
  results.Branch("fit_fun", &fun1);
   
  vector<TH1D*> h_charge;
  vector<TGraph*> g_spe;
  // Run - Channel - Good/Bad - Auto peak - SNR - gain - spe ampl - spe_peak_peak
  vector<tuple<size_t, size_t, bool, bool, double, double, double, double>> res_tuple;

  for(size_t ch_index=0; ch_index<channels.size(); ch_index++) {
    
    channel = channels[ch_index];
    if (channel < channel_low || channel > channel_up) continue;
    std::cout << "\nReading channel: " << channel << std::endl;
   
   // if (channel >= 10400 && channel < 10500){
   //   prepulse_ticks = 629;
   //   int_low = 631;
   //   int_up  = 661;
   // }
   //
   // if (channel >= 10500 && channel < 10600){
   //   prepulse_ticks = 635;
   //   int_low = 637;
   //   int_up  = 667;
   // }
   //
   // if (channel >= 10700 && channel < 10800){
   //    prepulse_ticks = 612;
   //    int_low = 617;
   //    int_up  = 647;
   // }

    // Read the wfs for this channel and subtract the baseline
    read();
    if (wfs.size() == 0) continue;

    SelCalib_WF(wfs, sel_wf, prepulse_ticks, sat_low, sat_up, bsl);  
    if (sel_wf.size() == 0) continue;
    
    TH1D* hI = new TH1D();
    hI = BuildRawChargeHisto(sel_wf, int_wf, int_low, int_up, nbins);
    if (avoid_auto_peak == true){
      spe_low=pspe_low;
      spe_up=pspe_up;
    }
    TH1D* hFind = new TH1D(*hI); //Delcared to plot it later
    fgaus = set_charge_fit_function(hI, hFind, avoid_auto_peak);
    hI->Fit("fgaus", "R");
    hI->Draw();//gPad->Update();gPad->WaitPrimitive();
    
    bool goodness = 0;
    // GOODNESS of fit !!
    LoadFitParameters(fgaus);
    Avg_Sel_WF(sel_wf, spe_avg, int_wf, spe_low, spe_up);
    spe_ampl = *max_element(spe_avg.begin()+int_low, spe_avg.begin()+int_up);
    double spe_undersh = *min_element(spe_avg.begin()+int_up, spe_avg.end());
   
    h1d = *hI; fun1 = *fgaus;
    results.Fill();

    hI->SetTitle(Form("Channel_%lu", channels[ch_index]));
    h_charge.push_back(hI);
    TGraph* gr = new TGraph(spe_avg.size(), data(spe_avg));
    gr->SetTitle(Form("Channel_%lu", channels[ch_index]));
    g_spe.push_back(gr);
    
    double SNR = fgaus->GetParameter(1)/fgaus->GetParameter(2); 
    double gain = fgaus->GetParameter(1); 
    if (SNR > 4.) goodness = 1;
    res_tuple.push_back(make_tuple(calibration_run, channel, goodness, avoid_auto_peak, SNR, gain, spe_ampl, spe_ampl-spe_undersh));

    int_wf.erase(int_wf.begin(), int_wf.end());
    sel_wf.erase(sel_wf.begin(), sel_wf.end());
    n_wf = 50000;
  } 

  for(auto tuple : res_tuple) print_tuple(tuple);
 
  // Save in file.root
  // results.Write();
  hf.cd("chargehistos"); for(auto h : h_charge) h->Write();
  hf.cd("spe_wfs");      for(auto g : g_spe)    g->Write();
  
  hf.Close();
}
