#include "../classe.hpp"

// ****************************************************************
// Read a calibration run file, build a charge histogram for each channel
// fit it, estimate the spe amplitude, store the histogram and the average spe
// waveform in a root file, print the results at termina//Read a calibration run file, build a charge histogram for each channel
// fit it, estimate the spe amplitude, store the histogram and the average spe
// waveform in a root file, print the results at terminal
// ****************************************************************


//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::ProtoDUNE_Calibration(){
  //////////////////////////////////////////////////////////////////////
  // HARD CODE HERE
  //
  calibration_run = 30003; //Needed create the root file with the results
  channel_low = 10500;     //Lower channel to look at (included)
  channel_up  = 12750;     //Upper " "
  mask = 0;
  int number_of_wf_here = 30000;
  //Output file name, then it adds "calibration_run.root":w
  string outfile_name = "tentative";
  avoid_auto_peak = true;
  pspe_low = 150; //Lower limit for spe integral (like spe_low)
  pspe_up  = 270;//Upper " " Remember: it depends on the integration window,
                    //the overvoltage and the gain. You can also enable the peak finding
                    //and don't use these
  //////////////////////////////////////////////////////////////////////

  //Choose the channels to read
  std::cout << "\n\n READING THE CH LIST \n\n" << std::endl;
  vector<size_t> channels = {10512};//, 11145, 11147};
  // vector<size_t> channels = read_pdhd_ch_map(mask);
  // vector<size_t> channels;
  // for(size_t i=1040; i<1045; i++){
  //   for(size_t j=0; j<8; j++){
  //     channels.push_back(i*10+j);
  //     std::cout << i*10+j << std::endl;
  //   }
  // }

  vector<vector<double>> sel_wf;
  vector<double> int_wf, spe_avg;

  outfile_name = outfile_name+calibration_run+".root";
  TFile hf(outfile_name.c_str(), "recreate");
  hf.mkdir("chargehistos");
  hf.mkdir("spe_wfs");
  
  //TTree results("calib_results", "calib_results");
  TH1D h1d; TF1 fun1;
  //results.Branch("chg_histo", &h1d);
  //results.Branch("fit_fun", &fun1);
   
  vector<TH1D*> h_charge_vec;
  vector<TGraph*> g_spe;
  // Run - Channel - Good/Bad - Auto peak - SNR - gain - spe ampl - spe_peak_peak
  vector<tuple<size_t, size_t, bool, bool, double, double, double, double>> res_tuple;
  vector<pair<string, double>> feature_value; // Store the results of the analysis to be printed 

  for(size_t ch_index=0; ch_index<channels.size(); ch_index++) {
    n_wf = number_of_wf_here;
    channel = channels[ch_index];
    if (channel < channel_low || channel > channel_up) continue;
    std::cout << "\nReading channel: " << channel << std::endl;
  
    /*
   if (channel >= 10400 && channel < 10500){
     prepulse_ticks = 629;
     int_low = 631;
     int_up  = 661;
   }

   if (channel >= 10500 && channel < 10600){
     prepulse_ticks = 635;
     int_low = 637;
     int_up  = 667;
   }

   if (channel >= 10700 && channel < 10800){
      prepulse_ticks = 612;
      int_low = 617;
      int_up  = 647;
   }
    */
    // Read the wfs for this channel and subtract the baseline
    read();
    if (wfs.size() == 0) continue;

    SelCalib_WF(wfs, sel_wf, prepulse_ticks, sat_low, sat_up, bsl);  
    if (sel_wf.size() == 0) continue;
    
    TH1D* hI = new TH1D();
    // hI = BuildRawChargeHisto(sel_wf, int_wf, int_low, int_up, nbins);
    hI = BuildRawChargeHisto(sel_wf, int_wf, int_low, int_up, hmin, hmax, nbins);
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
    //results.Fill();

    hI->SetTitle(Form("Channel_%lu", channels[ch_index]));
    h_charge_vec.push_back(hI);
    TGraph* gr = new TGraph(spe_avg.size(), data(spe_avg));
    gr->SetTitle(Form("Channel_%lu", channels[ch_index]));
    g_spe.push_back(gr);
    
    double SNR = fgaus->GetParameter(1)/fgaus->GetParameter(2); 
    double gain = fgaus->GetParameter(1); 
    if (SNR > 4.) goodness = 1;
    res_tuple.push_back(make_tuple(calibration_run, channel, goodness, avoid_auto_peak, SNR, gain, spe_ampl, spe_ampl-spe_undersh));
    

    feature_value.push_back({"Run",double(calibration_run)});
    feature_value.push_back({"Channel",double(channel)});
    feature_value.push_back({"Goodness",double(goodness)});
    feature_value.push_back({"Avoid_auto_peak",double(avoid_auto_peak)});
    feature_value.push_back({"SNR",SNR});
    feature_value.push_back({Form("Gain low %i up %i", int_low, int_up), gain});
    feature_value.push_back({"Spe ampl",spe_ampl});
    feature_value.push_back({"Spe peak-peak",spe_ampl-spe_undersh});
    int_wf.erase(int_wf.begin(), int_wf.end());
    sel_wf.erase(sel_wf.begin(), sel_wf.end());
  } 

  for(auto tuple : res_tuple) print_tuple(tuple);
  print_vec_pair_csv(Form("CalRun_%zu_pre%i_low%i_up%i_bsl%i.csv",calibration_run,prepulse_ticks,int_low,int_up,int(bsl)), feature_value);
 
  // Save in file.root
  // results.Write();
  hf.cd("chargehistos"); for(auto h : h_charge_vec) h->Write();
  hf.cd("spe_wfs");      for(auto g : g_spe)    g->Write();
  
  hf.Close();
}
