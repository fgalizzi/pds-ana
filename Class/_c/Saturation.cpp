#include "../classe.hpp"
void cla::Saturation() {
  double y0, y1;
  vector<double> x;
  vector<vector<double>> sat_wf;
  for (int i = 0; i < memorydepth; i++) x.push_back( (double) i);
  
  read(); 

  Sat_WF(wfs, sat_wf, sat_up);

  TCanvas* cTime = new TCanvas("wavedec","wavedec");
  
  min_max_element(sat_wf, y0, y1);

  std::cout << "Max " << y1 << "\n";
  std::cout << "Min " << y0 << "\n";
  std::cout << "Max - Min " << y1-y0 << "\n";

  TGraph *g1 = new TGraph(sat_wf[0].size(), &x[0], &sat_wf[0][0]);
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

  return;
}
