#include "../classe.hpp"

void cla::AverageWF() {

  double y0, y1, r_time, f_time, undershoot;
  vector<double> x, avg_wf, avg_wf2;
  vector<vector<double>> y2;
  for (size_t i = 0; i < memorydepth; i++) x.push_back( (double) i);
  //for (size_t i = 0; i < memorydepth; i++) x.push_back( (double) i*tick_len);
  
  // Read and subtract the baseline
  read();
  
  SelCalib_WF(wfs, y2, prepulse_ticks, sat_low, sat_up, bsl);
  
  if(mov_win == true) MovingAverageWF<double>(y2, y2, win);
  
  avgWF(y2 , avg_wf);

  std::cout << "WFS size " << wfs.size() << std::endl;  
  std::cout << "WF  size " << avg_wf.size() << std::endl;
  
  RiseFallTimeUndershoot(avg_wf, tick_len, r_time, f_time, undershoot);
  TCanvas* cTime = new TCanvas("pers","pers");
 
  y0 = *min_element(std::begin(avg_wf), std::end(avg_wf));
  y1 = *max_element(std::begin(avg_wf), std::end(avg_wf));
  
  std::cout << "\nUndershoot " << y0/y1*100. << "\n \n";
 
  min_max_element(y2, y0, y1);
  double ymin = y0; // y0, SAT_LOW, whatever
  double ymax = y1;

  TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "Persistence", "Ticks", "ADC Counts"), memorydepth/2, 0., memorydepth, 
      120, ymin, ymax);
  
  for (auto& wf : y2) for (int j=0; j<memorydepth; j=j+2) h2->Fill(j, wf[j]);
  
  h2->Draw("COLZ");
  gPad->SetLogz();

  TH1D* h_bsl = h2->ProjectionY("h_bsl", 0, prepulse_ticks/2);
  std::cout << "\nBaseline RMS = " << h_bsl->GetRMS() << std::endl; 

  //return;
  TGraph *g1 = new TGraph(avg_wf.size(), &x[0], &avg_wf[0]);
  g1->GetXaxis()->SetTitle("Time [#mus]");
  g1->GetYaxis()->SetTitle("ADC counts");
  g1->SetLineColor(2);
  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,800);
  c2->cd();
  g1->Draw();
  c2->Modified();
  c2->Update();
  gPad->Modified();
  gPad->Update();
//  VecDouble_in_Binary("Template50L.dat",avg_wf);
}
