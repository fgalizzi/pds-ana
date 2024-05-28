
void cla::Persistence() {
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);

  double y0, y1;
  vector<double> x , y2 , avg_wf, avg_wf2;
  
  for (size_t i = 0; i < MEMORYDEPTH; i++) x.push_back( (double) i);
    
  // Read and subtract the baseline
  read();

  std::cout << "WFS size " << WFS.size() << "\n #WF " << WFS.size()/MEMORYDEPTH << std::endl;  
  std::cout << "\n\navg size " << avg_wf.size() << std::endl;
  
  TCanvas* cTime = new TCanvas("wavedec","wavedec");
  
  //MovingAverageWF(WFS, y2, 100);
  DisplayWFs (WFS, MEMORYDEPTH, TICK_LEN, 10);

  y0 = *min_element(std::begin(WFS), std::end(WFS));
  y1 = *max_element(std::begin(WFS), std::end(WFS));
  //y0 = *min_element(std::begin(y2), std::end(y2));
  //y1 = *max_element(std::begin(y2), std::end(y2));
  std::cout << "Max - Min " << y1-y0 << "\n";

  // Change here y0->Whatever you want
  double ymin = y0;
  double ymax = y1;

  TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "DecoWF", "Ticks", "ADC Counts"), MEMORYDEPTH/2, 0., MEMORYDEPTH, 
      200, ymin, ymax);
  
  for (int i=0; i<N_WF; i++) for (int j=0; j<MEMORYDEPTH; j=j+2) h2->Fill(j, WFS[i*MEMORYDEPTH+j]);
  //for (int i=0; i<N_WF; i++) for (int j=0; j<MEMORYDEPTH; j=j+2) h2->Fill(j, y2[i*MEMORYDEPTH+j]);

  h2->Draw("COLZ");
  gPad->SetLogz();

  return;
}
