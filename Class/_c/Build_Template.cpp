#include "../classe.hpp"
// Define apply_filter, spe_low, spe_up before running this macro!
// You can do it with LED_Analysis.cpp or Filt_Analysis.cpp and updating these
// variables with LoadFitParameters(fgaus)

////////////////////////////////////////////////////////
/////// HARD CODE //////////////////////////////////////

string noise_td_file = "/Users/federico/PhD/PDE/Noise/Noise_NewEl_20241004_TimeDomain.dat";

////////////////////////////////////////////////////////


void cla::Build_Template() {
  double y0, y1, r_time, f_time, undershoot;
  vector<double> x, avg_calib, avg_template, int_wf, noise_td;
  vector<vector<double>> calib_wfs, template_wfs;
  size_t nsample = memorydepth;
  TComplex G[nsample];
    
  for (int i = 0; i < memorydepth; i++) x.push_back( (double) i*tick_len);
  
  // Read and subtract the baseline
  read();

  CompleteWF_Binary(noise_td_file, noise_td, memorydepth); // t_templ = time domain template
  SubVec_to_WFs(wfs, noise_td);

  SelCalib_WF(wfs, calib_wfs, prepulse_ticks, sat_low, sat_up, bsl);
  TH1D* hI = nullptr;
  hI = BuildRawChargeHisto(calib_wfs, int_wf, int_low, int_up, nbins);
  
  SelTemplate_WF(calib_wfs, template_wfs, int_wf, double(spe_low),
                 rms, int_low, int_up);

  avgWF(calib_wfs, avg_calib);
  avgWF(template_wfs, avg_template);

  double norm = *max_element(std::begin(avg_calib), std::end(avg_calib));
  norm = 1/norm;
  for(auto& e: avg_calib) e*=norm;
  norm = *max_element(std::begin(avg_template), std::end(avg_template));
  norm = 1/norm;
  for(auto& e: avg_template) e*=norm;

  if(print==true){
    VecDouble_in_Binary("Template.dat", avg_template);
    print = false;
  }

  if (plot == true){
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
