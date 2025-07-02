#include "../../../plotter/DUNEStyle.h"

#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"

using namespace std;

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////
vector<int> modules = {3}; // M1, M2, M3, M4

// --- INPUT -----------------------------------------------------
// --- OUTPUT ----------------------------------------------------

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void finddacvoltconversion(){
  dunestyle::SetDuneStyle();
  auto style = gROOT->GetStyle("duneStyle");
  style->SetOptFit(1111); style->SetOptTitle(0);
  style->SetStatX(0.9);   style->SetStatY(0.9);
  style->SetStatBorderSize(0);

  vector<double> dacs, volts;  

  int module = 3;
  if(module==1 || module==2 || module==11){
    dacs = {1148, 1161, 1174, 1187, 1200};
    volts = {45.06, 45.54, 46.09, 46.55, 47.03};
  } else {
//    dacs = {754, 767, 780, 793, 806};
//    volts = {30.54, 31.11, 31.63, 32.01, 32.52};
  dacs = {300, 400, 500};
  volts = {12.45, 16.38, 20.28};
  
  }
  vector<double> err_volts(volts.size(), 0.07);
  vector<double> err_dacs(dacs.size(), 0.);
  
  // --- DAC TO VOLT CONVERSION ---------------------------------
  TF1* f1 = new TF1 ( "f1", "[0]+[1]*x"); // y = q + m*x
  f1->SetParName(0, "Offset");
  f1->SetParName(1, "Slope");
  f1->SetParameter(0, 0.);
  f1->SetParameter(1, (volts[volts.size()-1]-volts[0])/(dacs[dacs.size()-1]-dacs[0]));
  
  TGraphErrors* g_Volt_DAC = new TGraphErrors(dacs.size(), &dacs[0], &volts[0],
                                              &err_dacs[0], &err_volts[0]);
  
  cout << "\n\n----------- FIT VOLT-DAC -------------------------\n" << endl;
  TFitResultPtr r1 = g_Volt_DAC->Fit(f1 ,"S");
  cout << "\n---------------------------------------------------  \n" << endl;
  
  // --- PLOT ---------------------------------------------------
  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);
  c1->cd();
  g_Volt_DAC->Draw("AP");
  c1->Modified(); c1->Update();
  
  return;
}
