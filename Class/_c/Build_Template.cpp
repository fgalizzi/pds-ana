#include "../classe.hpp"

// ****************************************************************
// Build a template selecting waveforms in two steps: selecting wfs
// with no light in the pre trigger and with no saturation and then
// the ones with integral > spe_low. The (baseline) rms parameter 
// help you in selecting the waveforms with no extra pulses (after
// pulses or scintillation light in coincidence)
// ****************************************************************

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////

string noise_td_file = "/Users/federico/PhD/PDE/Noise_files/Noise_NewEl_20241004_TimeDomain.dat";
string template_files_path = "/Users/federico/PhD/PDE/Template_files/";

///////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Build_Template() {
  vector<double> x, avg_calib, avg_template, int_wf, noise_td;
  vector<vector<double>> calib_wfs, template_wfs;
  size_t nsample = memorydepth;
  TComplex G[nsample];
    
  for (int i = 0; i < memorydepth; i++) x.push_back( (double) i*tick_len);
  
  // Read and subtract the baseline
  read();

  // Subtract the coherent noise of the digitiser (only once, if you re-run the macro)
  CompleteWF_Binary(noise_td_file, noise_td, memorydepth); // t_templ = time domain template
  // if(ite==0){
  //   SubVec_to_WFs(wfs, noise_td);
  //   ite++;
  // }

  // Select calibration waveforms (clean pre-pulse range and no saturation)
  // and compute in integral in [int_low;int_up]
  SelCalib_WF(wfs, calib_wfs, prepulse_ticks, sat_low, sat_up, bsl);
  TH1D* hI = nullptr;
  hI = BuildRawChargeHisto(calib_wfs, int_wf, int_low, int_up, nbins);
  
  // Select wfs to buil the template: integral>spe_low, no after-pulses nor 
  // light signal in the window. Tune rms for a more strict/loose requirement
  SelTemplate_WF(calib_wfs, template_wfs, int_wf, double(spe_low),
                 rms, int_low, int_up);

  // Averages for plots. "avg_template" is our template
  avgWF(calib_wfs, avg_calib);
  avgWF(template_wfs, avg_template);

  double norm; 
  if(print==true){
    string outfile_name;
    cout << "Name of the muon file (without .dat)" << endl;
    cin >> outfile_name;
    outfile_name = template_files_path+outfile_name;
    
    norm = 1./ *max_element(std::begin(avg_template), std::end(avg_template));
    for(auto& e: avg_template) e *= norm*spe_ampl;
    VecDouble_in_Binary(outfile_name, avg_template);
    print = false;
  }
  

  if (plot == true){
    // normalize the averages for the plots
    norm = 1./ *max_element(std::begin(avg_calib), std::end(avg_calib));
    for(auto& e: avg_calib) e*=norm;
    norm = 1./ *max_element(std::begin(avg_template), std::end(avg_template));
    for(auto& e: avg_template) e*=norm;
    
    TGraph *g_avg_calib = new TGraph(avg_calib.size(), &x[0], &avg_calib[0]);
    TGraph *g2 = new TGraph(avg_template.size(), &x[0], &avg_template[0]);
    g_avg_calib->GetXaxis()->SetTitle("Time [#mus]");
    g_avg_calib->GetYaxis()->SetTitle("ADC counts");
    g_avg_calib->SetLineColor(2);
    TCanvas *c_avg = new TCanvas("c_avg","c_avg",20,20,1000,800);
    c_avg->cd();
    g_avg_calib->Draw();
    g2->Draw("same");
    c_avg->Modified();c_avg->Update();

    double y_low, y_up;
    min_max_element(template_wfs, y_low, y_up);
    TH2D* h2_calib = new TH2D("h2", Form("%s;%s;%s", "Calibration Persistence",
                                          "Ticks", "ADC Counts"), memorydepth/2, 0., memorydepth, 
                                          120, y_low, y_up);
  
    for (auto& wf : calib_wfs)
      for (int j=0; j<memorydepth; j=j+2) h2_calib->Fill(j, wf[j]);
  
    TCanvas* c_calib_wfs = new TCanvas("c_calib_wfs","c_calib_wfs",0,0,800,600);
    c_calib_wfs->cd();
    c_calib_wfs->SetLogz();
    h2_calib->Draw("COLZ");
    c_calib_wfs->Modified(); c_calib_wfs->Update();

    TH2D* h2_templ = new TH2D("h2", Form("%s;%s;%s", "Template Persistence",
                                          "Ticks", "ADC Counts"), memorydepth/2, 0., memorydepth, 
                                          120, y_low, y_up);
  
    for (auto& wf : template_wfs)
      for (int j=0; j<memorydepth; j=j+2) h2_templ->Fill(j, wf[j]);
  
    TCanvas* c_templ_wfs = new TCanvas("c_templ_wfs","c_templ_wfs",0,0,800,600);
    c_templ_wfs->cd();
    c_templ_wfs->SetLogz();
    h2_templ->Draw("COLZ");
    c_templ_wfs->Modified(); c_templ_wfs->Update();

  }

}
