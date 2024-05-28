#include "flc.hpp"


auto INT_WIN = INT_UP-INT_LOW;
auto MOV_WIN = 200;
double X_UP = 4.;
double Y_UP = 1654;

void Sel_WF(const std::vector<double>& all_wf, std::vector<double>& x, std::vector<double>& y, int wfs, int len, int low, int up){
  double x0, x1, dev;
  double sum, sumsq, mean, variance;
  double pre = up-low;
  for(int i=0; i<N_WF; i++){
    auto in = all_wf.begin()+i*len+low;
    auto fin = all_wf.begin()+ i*len + up;
    sum = std::reduce(in, fin, 0.0);
    sumsq = std::accumulate(in, fin, double(),  [](double accumulator, double value) {return accumulator + value*value; });
    mean = sum/pre;
    x0 = *min_element(in, fin);
    x1 = *max_element(in, fin);
    //x.push_back(sumsq/pre - mean*mean);
    x.push_back(x1-x0);
    y.push_back(mean);
  }
  std::cout << sumsq << " " << len << " " << mean << std::endl;
}

TH1D* Good_WF(std::vector<double>& all_wf,  std::vector<double>& x, std::vector<double>& y, int wfs, int len){
  std::vector<double> sel_wf, int_wf;

  for (int i=0; i<N_WF; i++){
    if ( x[i]<X_UP && y[i]<Y_UP) for (int j=0; j<len; j++) sel_wf.push_back(all_wf[i*len+j]);
  }

  SubBaseline_Range(sel_wf, len, INT_LOW-INT_WIN-100, INT_UP-INT_WIN-100);
  //TH1D* h = BuildRawChargeHisto(sel_wf, int_wf, MEMORYDEPTH, INT_LOW, INT_UP);
  TH1D* h = BuildRawChargeHisto(sel_wf, int_wf, MEMORYDEPTH, INT_LOW, INT_UP, HMIN, HMAX);
  return h;
}

//*** MAIN ************************************
void MAX_Analysis(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");

  double t, x0, x1, y0, y1;
  std::vector<double> all_wf, int_wf, x, y, sm_wf;
  std::vector<double>::iterator aa;
  std::vector<double>::iterator bb;

  //CompleteWF_Binary(WF_FILE, all_wf, N_WF, MEMORYDEPTH);
  CAEN_WF_Binary(WF_FILE, all_wf, N_WF, MEMORYDEPTH);
  //SubBaseline(all_wf, MEMORYDEPTH, PREPULSE_TICKS);
 

  // A moving avegare just to select the whs without light before the LED PULSE
  MovingAverageWF(all_wf, sm_wf, MOV_WIN);
  Sel_WF(sm_wf, x, y, N_WF, MEMORYDEPTH, INT_LOW-INT_WIN+MOV_WIN-100, INT_UP-INT_WIN+MOV_WIN-100);
  //Sel_WF(sm_wf, x, y, N_WF, MEMORYDEPTH, INT_LOW+MOV_WIN, INT_UP+MOV_WIN);

  x0 = *min_element(std::begin(x), std::end(x));
  x1 = *max_element(std::begin(x), std::end(x));
  y0 = *min_element(std::begin(y), std::end(y));
  y1 = *max_element(std::begin(y), std::end(y));
  


/*
  //TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "DecoWF", "Ticks", "Amplitude [A.U.]"), 400, x0, x1, 400, y0, y1);
  TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "DecoWF", "Ticks", "Amplitude [A.U.]"), 400, 0., 50., 400, 1630., 1690.);

  for (int i=0; i<N_WF; i++)  h2->Fill(x[i], y[i]);
  


  TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
  c1->cd();
 
  h2->Draw("COLZ");
  c1->Modified();
  c1->Update();

  return;
*/



  //SubBaseline_Range(all_wf, MEMORYDEPTH, INT_LOW-800, INT_UP-800);
  //TH1D *h1 = h2->ProjectionY("h1");
  TH1D *h1 = Good_WF(all_wf, x, y, N_WF, MEMORYDEPTH);


  double par[NMAXPEAKS+4] = {0}; //Norm_Consts, sigma_0 , sigma_cell, G

    //go for peaks: create an instance of TSpectrum
  TSpectrum *s = new TSpectrum(NMAXPEAKS);
  int npeaks = s->Search(h1,1,"goff",0.05);
  par[0] = 0;     //peakx[0];
  par[1] = (SPE_LOW+SPE_UP)/2;   //peakx[1]-peakx[0];
  par[2] = (S0_LOW+S0_UP)/2;   //sigma_0
  par[3] = (SC_LOW+SC_UP)/2;    //sigma_cell
  for(int i = 0 ; i < npeaks ; i++){
    par[i+4] = 140;//peaky[i];
    //par[i*3+0] = peaky[i];
    //par[i*3+1] = peakx[i];
    //par[i*3+2] = 0.000000001;
  }

  //Gaussian sumatory
  TF1* fgaus = new TF1("fgaus", fRandomName, FIT_LOW, FIT_UP, npeaks+3);
  fRandomName_set(fgaus);
  fgaus->SetParameters(par);
  fgaus->SetParLimits(0, MU0_LOW , MU0_UP);
  //fgaus->FixParameter(0,0);
  fgaus->SetParLimits(1,SPE_LOW,SPE_UP);
  fgaus->SetParLimits(2,S0_LOW,S0_UP);
  fgaus->SetParLimits(3,SC_LOW,SC_UP);
  for(int i = 0 ; i < npeaks ; i++) fgaus->SetParLimits(i+4,0,700);

  TCanvas *c2 = new TCanvas("c2","c2",0,0,1000,900);
  c2->cd();
  h1->Draw();

  //fit histogram
  cout << "Fit ... " << endl;
  h1->Fit("fgaus", "R");
  cout << "... end fit. " << endl;

  c2->Modified();
  c2->Update();

}
