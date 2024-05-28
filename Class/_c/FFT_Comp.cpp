#include "flc.hpp"

const int F = 2;

void FFT_Comp() {
  double y0, y1;
  vector<double> x , y[F];
  
  for (size_t i = 0; i < MEMORYDEPTH; i++) {
    x.push_back( (double) i/(TICK_LEN*MEMORYDEPTH));
  }
 
  y0 = MEMORYDEPTH*0.5+1; //MEM = 25000
                          //
  CompleteWF_Binary("FFTwave0_rms_bdeon_notrun.dat", y[2], 1, (int)y0);
  CompleteWF_Binary("FFTwave0_rms_bdeon_run_ledon.dat", y[5], 1, (int)y0);
  //CompleteWF_Binary("FFTwave0_rms_bdeon_run.dat", y[3], 1, (int)y0);
  //CompleteWF_Binary("FFTwave0_rms_bdeon_runoff_ledon.dat", y[4], 1, (int)y0);
  //CompleteWF_Binary("FFTwave0_rms_bdewibon.dat", y[1], 1, (int)y0);
  //CompleteWF_Binary("FFTwave0_rms_dark_bde0ff.dat", y[0], 1, (int)y0);


  TCanvas* cTime = new TCanvas("wavedec","wavedec");
 
 
  //y0 = *min_element(std::begin(y), std::end(y));
  //y1 = *max_element(std::begin(y), std::end(y));
 
  TGraph* g_fft[F];

  for(int i=0; i<F; i++) {
    for (auto &e : y[i]) e = 20*TMath::Log10(e/(pow(2,14)));
    g_fft[i] = new TGraph((int)y0, &x[0], &y[i][0]);
  }
  auto mg = new TMultiGraph;
  
  for(int i=0; i<F; i++) mg->Add(g_fft[i]); 
  
  g_fft[0]->SetTitle("Dark BDE-off");
  g_fft[1]->SetTitle("Dark BDE WIB-on");
  g_fft[2]->SetTitle("Dark BDE-on Not Run");
  g_fft[3]->SetTitle("Dark BDE-on Run");
  g_fft[4]->SetTitle("LED BDE-on Not Run");
  g_fft[5]->SetTitle("LED BDE-on Run");


  TCanvas* c = new TCanvas("c", "Noise spectral density", 100, 100, 800, 600);
//  c->SetLogy(1);
  c->SetLogx(1);
  mg->Draw("A pmc plc");
  mg->GetXaxis()->SetTitle("Frequency [MHz]");
  mg->GetYaxis()->SetTitle("Power Spectral Density");
  c->BuildLegend();
  c->Modified(); c->Update();
//  VecDouble_in_Binary("Template50L.dat",avg_wf);
}
