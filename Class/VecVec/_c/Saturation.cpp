
void cla::Saturation() {
  double y0, y1;
  vector<double> x , sat_wf;
  
  for (size_t i = 0; i < MEMORYDEPTH; i++) x.push_back( (double) i);
  
  if(WFS.size()!=MEMORYDEPTH*N_WF){
    WFS.erase(WFS.begin(),WFS.end());
    if(DATA_FORMAT == "caen")   CAEN_WF_Binary(WF_FILE, WFS, N_WF, MEMORYDEPTH);
    if(DATA_FORMAT == "daphne") CompleteWF_Binary(WF_FILE, WFS, N_WF, MEMORYDEPTH);
    if(DATA_FORMAT == "csv")    CSV_WF_Binary(WF_FILE, WFS, N_WF, MEMORYDEPTH);
    
    if(INVERT == false) SubBaseline(WFS, MEMORYDEPTH, PREPULSE_TICKS);
    if(INVERT == true ) SubBaseline_Invert(WFS, MEMORYDEPTH, PREPULSE_TICKS);
  } 

  Sat_WF(WFS, sat_wf, MEMORYDEPTH, SAT_UP);

  TCanvas* cTime = new TCanvas("wavedec","wavedec");
  
  y0 = *min_element(std::begin(sat_wf), std::end(sat_wf));
  y1 = *max_element(std::begin(sat_wf), std::end(sat_wf));
  
  std::cout << "Max " << y1 << "\n";
  std::cout << "Min " << y0 << "\n";
  std::cout << "Max - Min " << y1-y0 << "\n";

  TGraph *g1 = new TGraph(sat_wf.size(), &x[0], &sat_wf[0]);
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
