//
//  G_Func.hpp
//
//  Created by Federico Galizzi on 28/09/23
//

#ifndef G_Func_hpp
#define G_Func_hpp

#include <stdio.h>

// Charge Histogram - Calibration Function (multi-Gaussians)
//**********************************************************
double fRandomName(Double_t *x, Double_t *par){
//**********************************************************

  double mu0  = par[0];
  double gain = par[1];
  double sg0  = par[2];
  double sg1  = par[3];

  Double_t result = 0;
  for (int i = 0 ; i < 10; i++){
    Double_t norm  = par[i+4];
    Double_t mean  = mu0+i*gain;
    Double_t sigma = sqrt(pow(sg0,2)+i*pow(sg1,2));
    result += norm*TMath::Gaus(x[0],mean,sigma);
  }
  return result;
}

void fRandomName_set(TF1* &f){
  f->SetNpx(2000);
  f->SetParNames("#mu_{0}" , "G" , "#sigma_{0}" , "#sigma_{cel}");
}


// Function "B_ik" for the Cross-talk probability
//**********************************************************
double fB(int i, int k){
//**********************************************************
  if (i==0 && k==0) return 1;
  if (i==0 && k>0)  return 0;
  return TMath::Factorial(k-1)/( TMath::Factorial(i)*TMath::Factorial(i-1)*TMath::Factorial(k-i) );
}

// Function "P_k(l,p)" for the Cross-talk probability
//**********************************************************
double fPK(int i, int k, double L, double p){
//**********************************************************
  return exp(-L)*fB(i,k)*pow(L*(1-p),i)*pow(p,k-i);
}

// Function for the Cross-talk estimation (x == k-th peak)
//**********************************************************
double fCX(double *x, double *par){
//**********************************************************
  double result = 0;
  for (int i=0; i<=x[0]; i++) result += fPK(i, (int)x[0], par[0], par[1]);
  return result;
}

void fCX_set(TF1 *f){
  f->SetNpx(2000);
  f->SetParNames("L", "p");
  f->SetParLimits(1, 0., 1.);
}

// Function for deconvolved LAr scint. light waveform +c
//**********************************************************
double fScintLight(double *x, double *p) {
//**********************************************************
  double val =(p[0]*exp(-(x[0]-p[5])/p[1])*exp(p[4]*p[4]/(2*p[1]*p[1]))) *TMath::Erfc(((p[5]-x[0])/p[4]+p[4]/p[1])/TMath::Power(2,0.5))/2. + (p[2]*exp(-(x[0]-p[5])/p[3]) *exp(p[4]*p[4]/(2*p[3]*p[3])))*TMath::Erfc(((p[5]-x[0])/p[4]+p[4]/p[3])/TMath::Power(2,0.5))/2. + p[6];
  return val;
}

void fScintLight_set(TF1 *f){
  f->SetNpx(2000);
  f->SetParNames("A_{s}", "#tau_{s}", "A_{t}", "#tau_{t}", "#sigma", "t_{0}", "c");
}



#endif /* G_Func_hpp */
