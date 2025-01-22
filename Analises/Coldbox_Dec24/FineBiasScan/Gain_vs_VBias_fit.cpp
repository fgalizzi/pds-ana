#include "../../../Class/_c/class_include.hpp"
#include "../../../plotter/DUNEStyle.h"
using namespace std;

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////
vector<int> modules = {11,2,3,4}; // M1, M2, M3, M4
vector<pair<double,double>> channels = {{20,27},{21,26},{0,2},{1,3}}; // M1 (20,27), M2 (21,26), M3 (0,2), M4 (1,3)
double err_volt_m1_m2 = 0.03; // to have chi2 ~1
double err_volt_m3_m4 = 0.07; // to have chi2 ~1
// --- INPUT -----------------------------------------------------
TString input_ana_folder = "/eos/home-f/fegalizz/ColdBox_VD/December24/Daphne_DAQ/FineBiasScan/";

size_t channel_colunm = 0;
size_t biases_colunm = 1;
size_t gains_colunm = 7;
size_t err_gains_colunm = 8;
size_t spe_ampls_colunm = 9;
size_t drs_colunm = 10;
size_t snrs_colunm = 11;
size_t err_snrs_colunm = 12;
size_t cxs_colunm = 13;
size_t err_cxs_colunm = 14;
size_t navg_cx_phs_colunm = 15;
size_t err_navg_cx_phs_colunm = 16;
size_t navg_phs_colunm = 17;
size_t navg_pes_colunm = 18;

// --- OUTPUT ----------------------------------------------------
TString output_folder = "/eos/home-f/fegalizz/ColdBox_VD/December24/Daphne_DAQ/";

///////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void Gain_vs_VBias_fit(){
  dunestyle::SetDuneStyle();
  auto style = gROOT->GetStyle("duneStyle");
  style->SetOptFit(1111); style->SetOptTitle(0);
  style->SetStatX(0.9);   style->SetStatY(0.9);
  style->SetStatBorderSize(0);

  // Create and open a file.root where to store the canvas
  TFile* file = new TFile(output_folder+"Gain_vs_VBias_fit.root", "RECREATE");
  file->cd();
  // Create a folder to store the canvas
  TDirectory* dir_gains = file->mkdir("Gains");
  TDirectory* dir_spe_ampl = file->mkdir("SPE_amplitude");
  TDirectory* dir_cx = file->mkdir("CX");
  TDirectory* dir_snr = file->mkdir("SNR");
  TDirectory* dir_dr = file->mkdir("DR");
  TDirectory* dir_navg_cx_ph = file->mkdir("Navg_cx_ph");
  TDirectory* dir_navg_ph = file->mkdir("Navg_ph");
  TDirectory* dir_navg_pe = file->mkdir("Navg_pe");
  TDirectory* dir_Vbr_canv= file->mkdir("VBr_canvas");

  vector<double> dacs, volts;  

  // --- LOOP OVER MODULES ---------------------------------------
  for(int idx_module=0; idx_module<modules.size(); idx_module++){
    int module = modules[idx_module];
    string input_file((input_ana_folder+Form("VBias_Scan_Module_%i.csv", module)).Data());
    vector<double> channel_this_module = {channels[idx_module].first, channels[idx_module].second};
    
    double err_volt;
    if(module==1 || module==2 || module==11){
      dacs = {1148, 1161, 1174, 1187, 1200};
      volts = {45.06, 45.54, 46.09, 46.55, 47.03};
      err_volt = err_volt_m1_m2;
    } else {
      dacs = {754, 767, 780, 793, 806};
      volts = {30.54, 31.11, 31.64, 32.01, 32.52};
      err_volt = err_volt_m3_m4;
    }
    vector<double> err_volts(volts.size(), err_volt);
    vector<double> err_dacs(dacs.size(), 0.);
    
    for(int idx_channel=0; idx_channel<channel_this_module.size(); idx_channel++){
       // --- VARIABLES -----------------------------------------------
      double q, m; // y = q + m*x
      vector<double> gains, biases, err_gains, err_biases, err_zeros;
      vector<double> spe_ampls, drs, snrs, cxs, navg_cx_phs, navg_phs, navg_pes, err_snrs, err_cxs, err_navg_cx_phs;
      double channel = channel_this_module[idx_channel];

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
        err_gains.push_back(input_data[err_gains_colunm].second[j]);
        // Convert Bias DAC in Volts and propagate the error
        double this_bias_dac[1] = {input_data[biases_colunm].second[j]};
        biases.push_back(f1->Eval(this_bias_dac[0]));
        double err_this_bias[1];
        r1->GetConfidenceIntervals(1, 1, 1, this_bias_dac, err_this_bias, 0.683, false);
        err_biases.push_back(err_this_bias[0]);
        
        err_zeros.push_back(0.); // for variable with no error

        spe_ampls.push_back(input_data[spe_ampls_colunm].second[j]);
        drs.push_back(input_data[drs_colunm].second[j]);
        snrs.push_back(input_data[snrs_colunm].second[j]);
        cxs.push_back(input_data[cxs_colunm].second[j]);
        navg_cx_phs.push_back(input_data[navg_cx_phs_colunm].second[j]);
        navg_phs.push_back(input_data[navg_phs_colunm].second[j]);
        navg_pes.push_back(input_data[navg_pes_colunm].second[j]);
        err_snrs.push_back(input_data[err_snrs_colunm].second[j]);
        err_cxs.push_back(input_data[err_cxs_colunm].second[j]);
        err_navg_cx_phs.push_back(input_data[err_navg_cx_phs_colunm].second[j]);

      }


      // --- FIT -----------------------------------------------------
      TGraphErrors* g_Gain_VBias = new TGraphErrors(biases.size(), &biases[0], &gains[0],
                                                    &err_biases[0], &err_gains[0]);

      TF1* f2 = new TF1 ( "f2", "[0]+[1]*x"); // y = q + m*x
      f2->SetNpx(5000);
      f2->SetParName(0, "Offset");
      f2->SetParName(1, "Slope");
      f2->SetParameter(0, 0.);
      f2->SetParameter(1, (gains[gains.size()-1]-gains[0])/(biases[biases.size()-1]-biases[0]));

      cout << "\n\n----------- FIT GAIN-BIAS ------------------------\n" << endl;
      TFitResultPtr r = g_Gain_VBias->Fit(f2 ,"S");
      r->Print("V");
      TMatrixDSym cov = r->GetCovarianceMatrix();
      cout << "\n---------------------------------------------------  \n" << endl;

      q = f2->GetParameter(0); m = f2->GetParameter(1);
      double v_br = -q/m;
      double v_br_low, v_br_up, err_v_br;

      cout << "\nV_br " << v_br << " +/- " << err_v_br << endl;
     
      TH1D *h_Confidence = new TH1D("h_Confidence", "h_Confidence", 10000, v_br-1., v_br+1.);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(h_Confidence);
      // loop over the h_Confidence bins and find when content+error is 0
      for (int i = 0; i < h_Confidence->GetNbinsX(); i++) {
        if (h_Confidence->GetBinContent(i)+h_Confidence->GetBinError(i) > 0) {
          v_br_low = h_Confidence->GetBinCenter(i-1);
          break;
        }
      }
      for (int i = h_Confidence->GetNbinsX(); i > 0; i--) {
        if (h_Confidence->GetBinContent(i)-h_Confidence->GetBinError(i) < 0) {
          v_br_up = h_Confidence->GetBinCenter(i+1);
          break;
        }
      }
      cout << "V_br_low " << v_br_low << " V_br_up " << v_br_up << endl;
      err_v_br = (v_br_up-v_br_low)/2.;

      // --- MANY TGRAPH ---------------------------------------------
      vector<double> overvoltages, err_overvoltages;
      for(size_t i=0; i<biases.size(); i++){
        overvoltages.push_back(biases[i]-v_br);
        err_overvoltages.push_back(error_propagation(biases[i], err_biases[i], v_br, err_v_br, "sub"));
      }

      file->cd("Gains");
      g_Gain_VBias->SetName(Form("M%i_Ch%i", module, int(channel)));
      g_Gain_VBias->SetTitle(Form("M%i Ch.%i", module, int(channel)));
      g_Gain_VBias->GetXaxis()->SetTitle("Bias [V]");
      g_Gain_VBias->GetYaxis()->SetTitle("Gain [ADC#times ticks]");
      g_Gain_VBias->Write();

      file->cd("SPE_amplitude");
      TGraphErrors* g_SPEampl_OV = new TGraphErrors(biases.size(), &overvoltages[0], &spe_ampls[0],
                                                    &err_overvoltages[0], &err_zeros[0]);
      g_SPEampl_OV->SetName(Form("M%i_Ch%i", module, int(channel)));
      g_SPEampl_OV->SetTitle(Form("M%i Ch.%i", module, int(channel)));
      g_SPEampl_OV->GetXaxis()->SetTitle("Overvoltage [V]");
      g_SPEampl_OV->GetYaxis()->SetTitle("SPE amplitude [ADC]");
      g_SPEampl_OV->Write();

      file->cd("CX");
      TGraphErrors* g_CX_OV = new TGraphErrors(biases.size(), &overvoltages[0], &cxs[0],
                                                    &err_overvoltages[0], &err_cxs[0]);
      g_CX_OV->SetName(Form("M%i_Ch%i", module, int(channel)));
      g_CX_OV->SetTitle(Form("M%i Ch.%i", module, int(channel)));
      g_CX_OV->GetXaxis()->SetTitle("Overvoltage [V]");
      g_CX_OV->GetYaxis()->SetTitle("CX");
      g_CX_OV->Write();

      file->cd("SNR");
      TGraphErrors* g_SNR_OV = new TGraphErrors(biases.size(), &overvoltages[0], &snrs[0],
                                                    &err_overvoltages[0], &err_snrs[0]);
      g_SNR_OV->SetName(Form("M%i_Ch%i", module, int(channel)));
      g_SNR_OV->SetTitle(Form("M%i Ch.%i", module, int(channel)));
      g_SNR_OV->GetXaxis()->SetTitle("Overvoltage [V]");
      g_SNR_OV->GetYaxis()->SetTitle("SNR");
      g_SNR_OV->Write();

      file->cd("DR");
      TGraphErrors* g_DR_OV = new TGraphErrors(biases.size(), &overvoltages[0], &drs[0],
                                                    &err_overvoltages[0], &err_zeros[0]);
      g_DR_OV->SetName(Form("M%i_Ch%i", module, int(channel)));
      g_DR_OV->SetTitle(Form("M%i Ch.%i", module, int(channel)));
      g_DR_OV->GetXaxis()->SetTitle("Overvoltage [V]");
      g_DR_OV->GetYaxis()->SetTitle("DR");
      g_DR_OV->Write();

      file->cd("Navg_cx_ph");
      TGraphErrors* g_Navg_CX_PH_OV = new TGraphErrors(biases.size(), &overvoltages[0], &navg_cx_phs[0],
                                                    &err_overvoltages[0], &err_navg_cx_phs[0]);
      g_Navg_CX_PH_OV->SetName(Form("M%i_Ch%i", module, int(channel)));
      g_Navg_CX_PH_OV->SetTitle(Form("M%i Ch.%i", module, int(channel)));
      g_Navg_CX_PH_OV->GetXaxis()->SetTitle("Overvoltage [V]");
      g_Navg_CX_PH_OV->GetYaxis()->SetTitle("Navg_CX_PH");
      g_Navg_CX_PH_OV->Write();

      file->cd("Navg_ph");
      TGraphErrors* g_NavgPH_OV = new TGraphErrors(biases.size(), &overvoltages[0], &navg_phs[0],
                                                    &err_overvoltages[0], &err_zeros[0]);
      g_NavgPH_OV->SetName(Form("M%i_Ch%i", module, int(channel)));
      g_NavgPH_OV->SetTitle(Form("M%i Ch.%i", module, int(channel)));
      g_NavgPH_OV->GetXaxis()->SetTitle("Overvoltage [V]");
      g_NavgPH_OV->GetYaxis()->SetTitle("Navg_PH");
      g_NavgPH_OV->Write();

      file->cd("Navg_pe");
      TGraphErrors* g_NavgPE_OV = new TGraphErrors(biases.size(), &overvoltages[0], &navg_pes[0],
                                                    &err_overvoltages[0], &err_zeros[0]);
      g_NavgPE_OV->SetName(Form("M%i_Ch%i", module, int(channel)));
      g_NavgPE_OV->SetTitle(Form("M%i Ch.%i", module, int(channel)));
      g_NavgPE_OV->GetXaxis()->SetTitle("Overvoltage [V]");
      g_NavgPE_OV->GetYaxis()->SetTitle("Navg_PE");
      g_NavgPE_OV->Write();

      // --- CONFIDENCE INTERVAL -------------------------------------


      h_Confidence = new TH1D("h_Confidence", "h_Confidence", 1000, v_br, biases[biases.size()-1]);
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

      TString canvas_name(Form("VBreadown_Module_%i_Channel_%i", module, int(channel)));
      TCanvas* c_Gain_VBias = new TCanvas(canvas_name,canvas_name,0,0,800,600);
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
      pave->AddText(Form("M%i Ch.%i - V Breakdown = %.2f #pm %.2f", module, int(channel), v_br, err_v_br));
      pave->Draw();

      // h_Confidence->Draw("e3");
      c_Gain_VBias->Modified(); c_Gain_VBias->Update();
      
      file->cd("VBr_canvas");
      c_Gain_VBias->Write();
       
    }
  }

  file->Close();

}
