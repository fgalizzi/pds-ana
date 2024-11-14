#include "../classe.hpp"
#include <cmath>
#include <ctime>
#include <vector>

const unsigned NSAMPLE = 5000;

double func_fit(double *x, double *p) {
  double val =(p[0]*exp(-(x[0]-p[5])/p[1])*exp(p[4]*p[4]/(2*p[1]*p[1])))*TMath::Erfc(((p[5]-x[0])/p[4]+p[4]/p[1])/TMath::Power(2,0.5))/2. + (p[2]*exp(-(x[0]-p[5])/p[3])*exp(p[4]*p[4]/(2*p[3]*p[3])))*TMath::Erfc(((p[5]-x[0])/p[4]+p[4]/p[3])/TMath::Power(2,0.5))/2. + p[6];
  return val;
}

void cla::Muon_PDHD(){

  double t=0;
 
  templ_f = "/Users/federico/PhD/PDE/Templates/Template.dat";
  muon_f = "/Users/federico/PhD/PDE/Muon.dat";
  std::vector<std::vector<double>> templ_v, avg_muon_v;
  std::vector<double> templ, avg_muon, y, noise;

  // Create the template vector
  CompleteWF_Binary(templ_f, templ_v, 1, NSAMPLE);
  CompleteWF_Binary(muon_f, avg_muon_v, 1, NSAMPLE);
  
  for(size_t i=0; i<NSAMPLE; i++){
    templ.push_back(templ_v[0][i]);
    avg_muon.push_back(avg_muon_v[0][i]);
  }

  double max_avg = *max_element(std::begin(templ), std::end(templ));
  for (int i=0; i<NSAMPLE; i++) templ[i] = templ[i]*(spe_ampl/max_avg);
  
  // prepare auxiliary arrays
  double t0 = 0.;      // in us
  double t1 = tick_len*NSAMPLE;    // in us
  const unsigned nsample = NSAMPLE;
  const int tmpl_prepulse_tick = prepulse_ticks; // template pre-pulse ticks:w
  t = t0;
  double xt[NSAMPLE] = {0}; // time array

  std::vector<Double_t> time;// time vector
  double xv[nsample] = {0}; // waveform array
  double xh[nsample] = {0}; // impulse response function array (spe template)
  double xs[nsample] = {0}; // original signal (δ-like function)
  double xm[nsample] = {0}; // template array
  double* xy;               // deconvoluted signal
  
  // fill the above arrays
  for (int i=0; i<nsample; i++) {
    xv[i] = avg_muon[i];     //Deconvolution of the avg_muon wf
    xm[i] = 0.;
    xs[i] = NSAMPLE*TMath::Gaus(i, 6, 0.09, false); // signal = delta function
    xt[i] = t;
    time.push_back(t);
    t+=tick_len;
   }
  
  for (int ip=0; ip<NSAMPLE; ip++) xh[ip] = templ[ip]; // template=spe
  
  //******************************
  //  Perform FFT
  //******************************
  int nsample_ = nsample;
  // xV: FFT of waveform
  TComplex xV[nsample]; double xV_re[nsample]; double xV_im[nsample];  //waveform
  // xM: FFT of waveform
  TComplex xM[nsample]; double xM_re[nsample]; double xM_im[nsample];  //template
  // xH: FFT of spe response
  TComplex xH[nsample]; double xH_re[nsample]; double xH_im[nsample];  //spe
  // xS: FFT of original signal
  TComplex xS[nsample]; double xS_re[nsample]; double xS_im[nsample];  //delta
  // xY: FFT of the filtered signal
  TComplex xY[nsample]; double xY_re[nsample]; double xY_im[nsample];  //filtered signal
  
  // Instance the FFT engine
  TVirtualFFT* fft = TVirtualFFT::FFT(1, &nsample_, "M R2C");
  
  // *************************FFT Waveform*******************
  fft->SetPoints(xv);
  fft->Transform();
  fft->GetPointsComplex(xV_re, xV_im);
  
  // ************************FFT Template********************
  fft->SetPoints(xm);
  fft->Transform();
  fft->GetPointsComplex(xM_re, xM_im);
  
  // ************************FFT SPE********************
  fft->SetPoints(xh);
  fft->Transform();
  fft->GetPointsComplex(xH_re, xH_im);
  
  //  *********************FFT DELTA *******************
  fft->SetPoints(xs);
  fft->Transform();
  fft->GetPointsComplex(xS_re, xS_im);

  // Fill FFT arrays and perform Wiener deconvolution
  double c_scale = 1./nsample;
  double H2[nsample] = {0}; // spe response spectral density
  double M2[nsample] = {0}; // template spectral density
  double N2[nsample] = {0}; // noise spectral density
  double S2[nsample] = {0}; // original signal spectral density
  double G2[nsample] = {0}; // Wiener filter profile
  double xf[nsample] = {0}; // frequency array
  TComplex G[nsample];
  
  for (int i=0; i<nsample*0.5+1; i++) {
    // fill FFT arrays
    xH[i] = TComplex(xH_re[i], xH_im[i])* c_scale;
    xV[i] = TComplex(xV_re[i], xV_im[i])* c_scale;
    xS[i] = TComplex(xS_re[i], xS_im[i])* c_scale;
    xM[i] = TComplex(xM_re[i], xM_im[i])* c_scale;
    // cout << "x2" << xH[i] <<endl ;
    
    // Compute spectral density
    H2[i] = xH[i].Rho2();
    M2[i] = xM[i].Rho2();
    N2[i] = n2_;
    S2[i] = xS[i].Rho2();
    
  //******************************
  // Compute Wiener filter
  //******************************
    xf[i] = i/t1;
    G[i]  = TComplex::Conjugate(xH[i])*S2[i] / (H2[i]*S2[i] + N2[i]); 
    //G[i]  = TComplex::Conjugate(xH[i])/H2[i]; // If you want to switch to 1/H
    
    // Compute filtered signal
    xY[i] = G[i]*xV[i];
    xY_re[i] = xY[i].Re(); xY_im[i] = xY[i].Im();
  }

  //*****************************************
  // Backward transform of the filtered signal
  //*****************************************
  
  fft = TVirtualFFT::FFT(1, &nsample_, "M C2R");
  fft->SetPointsComplex(xY_re, xY_im);
  fft->Transform();
  xy = fft->GetPointsReal();

  vector<double> deco_wf(memorydepth, 0.0);
  vector<double> e_x(memorydepth, 0.0);
  vector<double> e_y(memorydepth, 0.0);
  for (int i=0; i<nsample; i++){
    deco_wf[i] = xy[i]*0.01;
    e_x[i] = 0.004;
    e_y[i] = sqrt(abs(deco_wf[i]));
    std::cout << deco_wf[i] << " " << e_y[i] << std::endl;
  }

  rotate(deco_wf.begin(), deco_wf.begin()+deco_wf.size()-prepulse_ticks, deco_wf.end());
  //Consider to subtract the baseline

  //TF1* f1 = new TF1("f1", func_fit , FIT_L , FIT_U , 7);
  TF1 *f1 = new TF1("f1","([0]*exp(-(x-[5])/[1])*exp([4]*[4]/(2*[1]*[1])))*TMath::Erfc((([5]-x)/[4]+[4]/[1])/TMath::Power(2,0.5))/2. + ([2]*exp(-(x-[5])/[3])*exp([4]*[4]/(2*[3]*[3])))*TMath::Erfc((([5]-x)/[4]+[4]/[3])/TMath::Power(2,0.5))/2. + [6]",
                    fit_l,fit_u);
  f1->SetParameters(a_fast, tau_fast, a_slow, tau_slow, sigma, t_0);
  f1->SetParNames("A_{s}", "#tau_{s}", "A_{t}", "#tau_{t}", "#sigma", "t_{0}", "c");
  f1->SetNpx(2000);
  if(FIX_CONST == true) f1->FixParameter( 6 , 0. );

  TGraphErrors* g_er = new TGraphErrors(time.size(), &time[0], &deco_wf[0], &e_x[0], &e_y[0]);
  TGraph* gy = new TGraphErrors(time.size(), &time[0], &deco_wf[0]);
  
  TGraphSmooth* sk1 = new TGraphSmooth("normal");
  gy = sk1->SmoothKern(g_er, "normal", deco_sm);

  if (FFUNC == true){
    auto fitResult = g_er->Fit(f1 ,"RS");
    f1->Draw("SAME");
  }
  if (FFUNC == false){ f1->Draw("SAME");}


  double A_s   = f1->GetParameter(0);
  double tau_s = f1->GetParameter(1);
  double A_t   = f1->GetParameter(2);
  double tau_t = f1->GetParameter(3);
  std::cout <<"\n \nFast/Slow intensities" << std::endl;
  std::cout <<"Is = " << A_s*tau_s/(A_s*tau_s+A_t*tau_t) << std::endl;  
  std::cout <<"It = " << A_t*tau_t/(A_s*tau_s+A_t*tau_t) << std::endl;  

  std::cout << "\n \nDeco_sm - As - tau_s - At - tau_t - Sigma - [ns]" << std::endl;  
  std::cout << deco_sm*1000 << "\t" << A_s << "\t" << tau_s*1000 << "\t" 
    << A_t << "\t" << tau_t*1000 << "\t" << f1->GetParameter(4)*1000 << std::endl;  





    





    //---------------P L O T S----------------------
      
    // display waveform and noise
    TGraph* gv = new TGraph(nsample, xt, xv);  //waveform
    TGraph* gs = new TGraph(nsample, xt, xs);  //delta funtion
    TGraph* gh = new TGraph(nsample, xt, xh);  //spe
    TGraph* gm = new TGraph(nsample, xt, xm);  //spe template
    
      
    gv->SetLineColor(kRed+1);
    gs->SetLineColor(kBlue+1);
    gh->SetLineColor(kMagenta+1);

    gStyle->SetOptTitle(0);
    TCanvas* cTime = new TCanvas("wavedec","wavedec");
    cTime->Divide(1, 2);
    cTime->cd(1);
    gv->Draw("awlx+");
    gh->Draw("same");
    gv->SetNameTitle("gv", "Syntetic waveform");
    gv->GetXaxis()->SetTitle("Time [#mus]");
    gv->GetYaxis()->SetTitle("Amplitude (ADC counts)");
    gv->GetXaxis()->CenterTitle();
    gv->GetYaxis()->CenterTitle();
    gv->GetYaxis()->SetTitleSize(0.06);
    gv->GetYaxis()->SetLabelSize(0.06);
    gv->GetXaxis()->SetTitleSize(0.06);
    gv->GetXaxis()->SetLabelSize(0.06);
    //gs->Draw("l");
    gPad->SetTopMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.01);
    gPad->SetTicks(1, 1); gPad->SetGrid(1, 1);



     
  // display delta and spe template in the time domain
  TCanvas* cDelta = new TCanvas("cDelta","cDelta");
  gh->Draw("aDwpl");
  gs->Draw("l");
  
  
  // Display (normalized) spectral densities
  TGraph* gN2 = new TGraph(nsample*0.5+1, xf, N2);
  TGraph* gH2 = new TGraph(nsample*0.5+1, xf, H2);
  TGraph* gM2 = new TGraph(nsample*0.5+1, xf, M2);
  TGraph* gS2 = new TGraph(nsample*0.5+1, xf, S2);

  TCanvas* cPower = new TCanvas();
  cPower->SetLogy(1);
  cPower->SetLogx(1);
  cPower->cd();
  gN2->SetLineColor(kGray+1);
  gH2->SetLineColor(kRed+1);
  gS2->SetLineColor(kBlue+1);
  gM2->SetLineColor(kOrange+1);
  gH2->GetXaxis()->SetTitle("Frequency [MHz]");
  gH2->GetYaxis()->SetTitle("Power Spectral Density");
  gH2->Draw("awl");
  gN2->Draw("l");
  gM2->Draw("l");
  gS2->Draw("l");

  cTime->cd(2);
  
  TCanvas *f_canv = new TCanvas("FitCanv","FitCanv",20,20,1000,900);
  f_canv->cd();
  
  //f_canv->SetLogy(1);
  g_er->GetXaxis()->SetTitle("Time [#mus]");
  g_er->GetYaxis()->SetTitle("Amplitude [A.U.]");
  g_er->SetTitle("Deconvolved waveform");
  //gy->GetXaxis()->CenterTitle();
  //gy->GetYaxis()->CenterTitle();
  //gy->GetYaxis()->SetTitleSize(0.06);
  //gy->GetYaxis()->SetLabelSize(0.06);
  //gy->GetXaxis()->SetTitleSize(0.06);
  //gy->GetXaxis()->SetLabelSize(0.06);
  g_er->SetLineColor(kBlue);
  g_er->SetLineWidth(2);
  g_er->Draw("al");
 

  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.14);
  gPad->SetTicks(1, 1); gPad->SetGrid(1, 1);
  gPad->BuildLegend(0.5, 0.88, 0.88, 0.8, "", "l");


}
