#include "../classe.hpp"

void cla::Persistence() {
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);

  double y0, y1;
  vector<double> x;
  vector<vector<double>> y2;
  
  for (size_t i = 0; i < memorydepth; i++) x.push_back( (double) i);
    
  // Read and subtract the baseline
  read();

  std::cout << "wfs size " << wfs.size() << "\nmemorydepth " << memorydepth << std::endl;  
  
  TCanvas* cTime = new TCanvas("wavedec","wavedec");
  
  //MovingAverageWF(WFS, y2, 100);
  DisplayWFs (wfs, tick_len, 10);

  min_max_element(wfs, y0, y1);
  //min_max_element(y2, y0, y1);
  std::cout << "Max - Min " << y1-y0 << "\n";

  // Change here y0->Whatever you want
  double ymin = y0;
  double ymax = y1;

  TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "DecoWF", "Ticks", "ADC Counts"), memorydepth/2, 0., memorydepth, 
      200, ymin, ymax);
  
  for (auto wf : wfs) for (int j=0; j<memorydepth; j=j+2) h2->Fill(j, wf[j]);

  h2->Draw("COLZ");
  gPad->SetLogz();

  return;
}
