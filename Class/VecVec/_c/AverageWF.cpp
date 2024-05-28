
void cla::AverageWF() {
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);

  double y0, y1;
  vector<double> x , y2 , avg_wf, avg_wf2;
  
  for (size_t i = 0; i < MEMORYDEPTH; i++) x.push_back( (double) i);
  //for (size_t i = 0; i < MEMORYDEPTH; i++) x.push_back( (double) i*TICK_LEN);
  
  // Read and subtract the baseline
  read();
  
  SelCalib_WF(WFS, y2, MEMORYDEPTH, PREPULSE_TICKS, SAT_LOW, SAT_UP, BSL);
  avgWF(y2 , avg_wf, MEMORYDEPTH);
  //avgWF(WFS , avg_wf, MEMORYDEPTH);
  //MovingAverageWF(avg_wf, avg_wf, 10);

  std::cout << "WFS size " << WFS.size() << "\n#WF " << WFS.size()/MEMORYDEPTH << std::endl;  
  std::cout << "avg size " << avg_wf.size() << std::endl;
  
  Print_RiseFallTime(avg_wf, TICK_LEN);
  TCanvas* cTime = new TCanvas("pers","pers");
 
  y0 = *min_element(std::begin(avg_wf), std::end(avg_wf));
  y1 = *max_element(std::begin(avg_wf), std::end(avg_wf));
  
  std::cout << "\nUndershoot " << y0/y1*100. << "\n \n";
 
  //y0 = *min_element(std::begin(WFS), std::end(WFS));
  //y1 = *max_element(std::begin(WFS), std::end(WFS));
  y0 = *min_element(std::begin(y2), std::end(y2));
  y1 = *max_element(std::begin(y2), std::end(y2));

  double ymin = y0; // y0, SAT_LOW, whatever
  double ymax = y1;

  TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "Persistence", "Ticks", "ADC Counts"), MEMORYDEPTH/2, 0., MEMORYDEPTH, 
      120, ymin, ymax);
  
  for (int i=0; i*MEMORYDEPTH < y2.size(); i++) for (int j=0; j<MEMORYDEPTH; j=j+2) h2->Fill(j, y2[i*MEMORYDEPTH+j]);
  
  h2->Draw("COLZ");
  gPad->SetLogz();

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
