#include "../../../Class/_c/class_include.hpp"
#include "../../../plotter/DUNEStyle.h"
using namespace std;

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////
int module = 3; 
double channel = 0; // M1 (20,27), M2 (21,26), M3 (0,2), M4 (1,3)

TString input_ana_folder = "/eos/home-f/fegalizz/ColdBox_VD/December24/Daphne_DAQ/FineBiasScan/";
string input_file((input_ana_folder+Form("VBias_Scan_Module_%i.csv", module)).Data());

size_t channel_colunm = 0;
size_t gains_colunm = 7;
size_t biases_colunm = 1;
size_t err_gains_colunm = 8;
// size_t err_biases_colunm = 0;

// M1 and M2
vector<double> dacs = {1148, 1161, 1174, 1187, 1200};
vector<double> volts = {45.06, 45.54, 46.09, 46.55, 47.03};
// M3 and M4
// vector<double> dacs = {754, 767, 780, 793, 806};
// vector<double> volts = {30.54, 31.11, 31.64, 32.01, 32.52};
vector<double> err_volts(volts.size(), 0.03);
vector<double> err_dacs(dacs.size(), 0.5);
///////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void Gain_vs_VBias_fit(){
  dunestyle::SetDuneStyle();
  auto style = gROOT->GetStyle("duneStyle");
  style->SetOptFit(1111); style->SetOptTitle(0);
  style->SetStatX(0.9);   style->SetStatY(0.9);
  style->SetStatBorderSize(0);

  // --- VARIABLES -----------------------------------------------
  double q, m; // y = q + m*x
  vector<double> gains, biases, err_gains, err_biases;

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


  // --- READ ----------------------------------------------------
  vector<pair<string, vector<double>>> input_data = read_vec_pair_CSV(input_file.c_str());
  for(size_t j=0; j<input_data[0].second.size(); j++){
    if (input_data[channel_colunm].second[j] != channel) continue;
    gains.push_back(input_data[gains_colunm].second[j]);
    biases.push_back(f1->Eval(input_data[biases_colunm].second[j]));
    err_gains.push_back(input_data[err_gains_colunm].second[j]);
    // err_biases.push_back(input_data[err_biases_colunm].second[j]);
    err_biases.push_back(0.015);
  }


  // --- FIT -----------------------------------------------------
  TGraphErrors* g_Gain_VBias = new TGraphErrors(biases.size(), &biases[0], &gains[0],
                                                &err_biases[0], &err_gains[0]);

  TF1* f2 = new TF1 ( "f2", "[0]+[1]*x"); // y = q + m*x
  f2->SetNpx(3000);
  f2->SetParName(0, "Offset" );
  f2->SetParName(1, "Slope" );
  f2->SetParameter(0, 0.);
  f2->SetParameter(1, (gains[gains.size()-1]-gains[0])/(biases[biases.size()-1]-biases[0]));

  cout << "\n\n----------- FIT GAIN-BIAS ------------------------\n" << endl;
  TFitResultPtr r = g_Gain_VBias->Fit(f2 ,"S");
  r->Print("V");
  TMatrixDSym cov = r->GetCovarianceMatrix();
  cout << "\n---------------------------------------------------  \n" << endl;

  q = f2->GetParameter(0); m = f2->GetParameter(1);
  double v_br = -q/m;
  double err_v_br = error_propagation(r, f2, 0, 1, "div");

  cout << "\n V_br " << v_br << " +/- " << err_v_br << endl;
 
  TH1D *h_Confidence = new TH1D("h_Confidence", "h_Confidence", 1000, -q/m, biases[biases.size()-1]);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(h_Confidence);
  h_Confidence->SetStats(false);  h_Confidence->SetFillColor(2);
  h_Confidence->SetMarkerSize(0); h_Confidence->SetLineWidth(0);
  h_Confidence->SetFillStyle(3001);
  h_Confidence->SetFillColorAlpha(kRed, 0.01);
   


  // --- PLOT ----------------------------------------------------
  TCanvas* c_Volt_DAC = new TCanvas("c_Volt_DAC","c_Volt_DAC",0,0,800,600);
  c_Volt_DAC->cd();
  g_Volt_DAC->GetXaxis()->SetTitle("DAC"); 
  g_Volt_DAC->GetYaxis()->SetTitle("V_{Bias} [V]");
  g_Volt_DAC->SetTitle("");
  g_Volt_DAC->Draw("AL");
  c_Volt_DAC->Modified(); c_Volt_DAC->Update();

  TCanvas* c_Gain_VBias = new TCanvas("c_Gain_VBias","c_Gain_VBias",0,0,800,600);
  c_Gain_VBias->cd();
  dunestyle::ApplyDuneStyle(g_Gain_VBias->GetHistogram());
  dunestyle::CenterTitles(g_Gain_VBias->GetHistogram());
  g_Gain_VBias->GetXaxis()->SetTitle("V_{Bias} [V]");
  g_Gain_VBias->GetYaxis()->SetTitle("Gain [ADC#times ticks]");
  g_Gain_VBias->GetYaxis()->SetTitleOffset(2);
  g_Gain_VBias->SetTitle("");
  g_Gain_VBias->SetLineWidth(1);
  g_Gain_VBias->Draw("AP");
  h_Confidence->Draw("e3 same");
  f2->Draw("same");
  g_Gain_VBias->Fit(f2, "QN");
  g_Gain_VBias->Draw("P same");

  auto pave = new TPaveText(0.6, 0.45, 0.85, 0.65, "NDC NB");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextSizePixels(35);
  pave->SetTextFont(6); pave->SetTextSize(18);
  pave->AddText(Form("V Breakdown = %.2f #pm %.2f", v_br, err_v_br));
  pave->Draw();

  c_Gain_VBias->Modified(); c_Gain_VBias->Update();
}
