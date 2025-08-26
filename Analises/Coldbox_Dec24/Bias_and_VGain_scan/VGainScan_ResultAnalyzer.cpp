#include "../../../Class/_c/class_include.hpp"
#include "../../../plotter/DUNEStyle.h"
#include "TGraphErrors.h"
using namespace std;

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////
//vector<int> modules = {1,2}; // M1, M2, M3, M4
vector<int> modules = {1,2,3,4}; // M1, M2, M3, M4
//vector<vector<int>> module_channels = {{20,27},{21,26}}; // M1 (20,27), M2 (21,26), M3 (0,2), M4 (1,3)
vector<vector<int>> module_channels = {{20,27},{21,26},{0,2},{1,3}}; // M1 (20,27), M2 (21,26), M3 (0,2), M4 (1,3)

vector<int> m1_m2_biases = {1148, 1161, 1174, 1187, 1200};
vector<int> m3_m4_biases = {754, 767, 780, 793, 806};

bool verbose = false;
// --- INPUT -----------------------------------------------------
TString ana_folder = "/eos/home-g/gpiemont/ColdBox_VD/December24/Daphne_DAQ/VGain_Scans/";

size_t channel_colunm = 2;
size_t bias_volt_colunm = 4;
size_t overvoltage_colunm = 6;
size_t vgain_colunm = 8;
size_t gain_colunm = 14;
size_t err_gain_colunm = 15;
size_t spe_ampl_colunm = 16;
size_t dr_colunm = 17;
size_t snr_colunm = 18;
size_t err_snr_colunm = 19;
size_t rms_colunm = 20;
size_t cx_colunm = 21;
size_t err_cx_colunm = 22;
size_t navg_cx_phs_colunm = 23;
size_t err_navg_cx_phs_colunm = 24;
size_t navg_phs_colunm = 25;
size_t navg_pes_colunm = 26;

///////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void VGainScan_ResultAnalyzer(){
  dunestyle::SetDuneStyle();
  auto style = gROOT->GetStyle("duneStyle");
  style->SetOptFit(1111); style->SetOptTitle(0);
  style->SetStatX(0.9);   style->SetStatY(0.9);
  style->SetStatBorderSize(0);

  // Create and open a file.root where to store the canvas
  TFile* file = new TFile(ana_folder+"VGainScans_ResultSummary_3RMS_noverbose.root", "RECREATE");
  file->cd();

  // --- LOOP OVER MODULES ---------------------------------------
  for (int idx_module=0; idx_module<modules.size(); idx_module++){
    int module = modules[idx_module];
    vector<int> channels = module_channels[idx_module];
    vector<int> biases;
    TString sipm;
    if(module==1 || module==2){
      biases = m1_m2_biases;
      sipm = "HPK";
    } else {
      biases = m3_m4_biases;
      sipm = "FBK";
    }

    // --- LOOP OVER CHANNELS ------------------------------------
    for (const auto& channel : channels){
      TDirectory* dir_ch = file->mkdir(Form("M%i_%s_Ch%i", module, sipm.Data(), channel));

      // --- LOOP OVER BIASES ------------------------------------
      for (const auto& bias : biases){
        dir_ch->cd();
        
        string input_file((ana_folder+Form("Module_%i_Bias_%i.csv", module, bias)).Data());
        vector<pair<string, vector<double>>> input_data = read_vec_pair_CSV(input_file.c_str());
        float bias_volt, overvoltage;
        for (size_t j=0; j<input_data[0].second.size(); j++){
          if (input_data[channel_colunm].second[j] == channel){
            bias_volt = input_data[bias_volt_colunm].second[j];
            overvoltage = input_data[overvoltage_colunm].second[j];
            break;
          }
        }

        TDirectory* dir_bias;
        if (verbose) dir_bias = dir_ch->mkdir(Form("BiasDAC_%i_BiasVolt_%.2f_OV_%.2f", bias, bias_volt, overvoltage));
        else         dir_bias = dir_ch->mkdir(Form("BiasDAC_%i", bias));
        dir_bias->cd();

        // --- VARIABLES -----------------------------------------
        vector<double> vgains, gains, spe_ampls, drs, snrs, rmses, cxs, navg_cx_phs, navg_phs, navg_pes;
        vector<double> err_gains, err_snrs, err_cxs, err_navg_cx_phs;
        for (size_t j=0; j<input_data[0].second.size(); j++){
          if (input_data[channel_colunm].second[j] != channel) continue;
          vgains.push_back(input_data[vgain_colunm].second[j]);
          gains.push_back(input_data[gain_colunm].second[j]);
          err_gains.push_back(input_data[err_gain_colunm].second[j]);
          spe_ampls.push_back(input_data[spe_ampl_colunm].second[j]);
          drs.push_back(input_data[dr_colunm].second[j]);
          snrs.push_back(input_data[snr_colunm].second[j]);
          err_snrs.push_back(input_data[err_snr_colunm].second[j]);
          rmses.push_back(input_data[rms_colunm].second[j]);
          cxs.push_back(input_data[cx_colunm].second[j]*100);
          err_cxs.push_back(input_data[err_cx_colunm].second[j]*100);
          navg_cx_phs.push_back(input_data[navg_cx_phs_colunm].second[j]);
          err_navg_cx_phs.push_back(input_data[err_navg_cx_phs_colunm].second[j]);
          navg_phs.push_back(input_data[navg_phs_colunm].second[j]);
          navg_pes.push_back(input_data[navg_pes_colunm].second[j]);
        }
        vector<double> err_zeros(vgains.size(), 0.);

        // --- TGRAPHS -------------------------------------------
        TGraphErrors* g_Gain_VGain = new TGraphErrors(vgains.size(), &vgains[0], &gains[0],
                                                      &err_zeros[0], &err_gains[0]);
        if (verbose) g_Gain_VGain->SetTitle(Form("Gain vs VGain - M%i-%s Ch.%i Bias= %.2f OV= %.2f ", module, sipm.Data(), channel, bias_volt, overvoltage));
        else         g_Gain_VGain->SetTitle(Form("Bias= %.2f OV= %.2f ", bias_volt, overvoltage));
        g_Gain_VGain->GetXaxis()->SetTitle("VGain");
        g_Gain_VGain->GetYaxis()->SetTitle("Gain [ADC#times ticks]");
        if (verbose) g_Gain_VGain->Write(Form("Gain_VGain_M%i_Ch%i_Bias%i", module, channel, bias));
        else         g_Gain_VGain->Write("Gain_VGain");

        TGraphErrors* g_SPEampl_VGain = new TGraphErrors(vgains.size(), &vgains[0], &spe_ampls[0],
                                                      &err_zeros[0], &err_zeros[0]);
        if (verbose) g_SPEampl_VGain->SetTitle(Form("SPE Amplitude vs VGain - M%i-%s Ch.%i Bias= %.2f OV= %.2f ", module, sipm.Data(), channel, bias_volt, overvoltage));
        else         g_SPEampl_VGain->SetTitle(Form("Bias= %.2f OV= %.2f ", bias_volt, overvoltage));
        g_SPEampl_VGain->GetXaxis()->SetTitle("VGain");
        g_SPEampl_VGain->GetYaxis()->SetTitle("SPE Amplitude [ADC]");
        if (verbose) g_SPEampl_VGain->Write(Form("SPEampl_VGain_M%i_Ch%i_Bias%i", module, channel, bias));
        else         g_SPEampl_VGain->Write("SPEampl_VGain");

        TGraphErrors* g_DR_VGain = new TGraphErrors(vgains.size(), &vgains[0], &drs[0],
                                                    &err_zeros[0], &err_zeros[0]);
        if (verbose) g_DR_VGain->SetTitle(Form("Dynamic Range vs VGain - M%i-%s Ch.%i Bias= %.2f OV= %.2f ", module, sipm.Data(), channel, bias_volt, overvoltage));
        else         g_DR_VGain->SetTitle(Form("Bias= %.2f OV= %.2f ", bias_volt, overvoltage));
        g_DR_VGain->GetXaxis()->SetTitle("VGain");
        g_DR_VGain->GetYaxis()->SetTitle("Dynamic Range [pe]");
        if (verbose) g_DR_VGain->Write(Form("DR_VGain_M%i_Ch%i_Bias%i", module, channel, bias));
        else         g_DR_VGain->Write("DR_VGain");

        TGraphErrors* g_SNR_VGain = new TGraphErrors(vgains.size(), &vgains[0], &snrs[0],
                                                     &err_zeros[0], &err_snrs[0]);
        if (verbose) g_SNR_VGain->SetTitle(Form("SNR vs VGain - M%i-%s Ch.%i Bias= %.2f OV= %.2f ", module, sipm.Data(), channel, bias_volt, overvoltage));
        else         g_SNR_VGain->SetTitle(Form("Bias= %.2f OV= %.2f ", bias_volt, overvoltage));
        g_SNR_VGain->GetXaxis()->SetTitle("VGain");
        g_SNR_VGain->GetYaxis()->SetTitle("SNR");
        if (verbose) g_SNR_VGain->Write(Form("SNR_VGain_M%i_Ch%i_Bias%i", module, channel, bias));
        else         g_SNR_VGain->Write("SNR_VGain");

        TGraphErrors* g_RMS_VGain = new TGraphErrors(vgains.size(), &vgains[0], &rmses[0],
                                                     &err_zeros[0], &err_zeros[0]);
        if (verbose) g_RMS_VGain->SetTitle(Form("RMS vs VGain - M%i-%s Ch.%i Bias= %.2f OV= %.2f ", module, sipm.Data(), channel, bias_volt, overvoltage));
        g_RMS_VGain->GetXaxis()->SetTitle("VGain");
        g_RMS_VGain->GetYaxis()->SetTitle("RMS [ADC]");
        if (verbose) g_RMS_VGain->Write(Form("RMS_VGain_M%i_Ch%i_Bias%i", module, channel, bias));
        else         g_RMS_VGain->Write("RMS_VGain");

        TGraphErrors* g_CX_VGain = new TGraphErrors(vgains.size(), &vgains[0], &cxs[0],
                                                    &err_zeros[0], &err_cxs[0]);
        if (verbose) g_CX_VGain->SetTitle(Form("Cross Talk vs VGain - M%i-%s Ch.%i Bias= %.2f OV= %.2f ", module, sipm.Data(), channel, bias_volt, overvoltage));
        else         g_CX_VGain->SetTitle(Form("Bias= %.2f OV= %.2f ", bias_volt, overvoltage));
        g_CX_VGain->GetXaxis()->SetTitle("VGain");
        g_CX_VGain->GetYaxis()->SetTitle("Cross Talk [%]");
        if (verbose) g_CX_VGain->Write(Form("CX_VGain_M%i_Ch%i_Bias%i", module, channel, bias));
        else         g_CX_VGain->Write("CX_VGain");

        TGraphErrors* g_Navg_CX_PH_VGain = new TGraphErrors(vgains.size(), &vgains[0], &navg_cx_phs[0],
                                                    &err_zeros[0], &err_navg_cx_phs[0]);
        if (verbose) g_Navg_CX_PH_VGain->SetTitle(Form("Navg_CX_PH vs VGain - M%i-%s Ch.%i Bias= %.2f OV= %.2f ", module, sipm.Data(), channel, bias_volt, overvoltage));
        else         g_Navg_CX_PH_VGain->SetTitle(Form("Bias= %.2f OV= %.2f ", bias_volt, overvoltage));
        g_Navg_CX_PH_VGain->GetXaxis()->SetTitle("VGain");
        g_Navg_CX_PH_VGain->GetYaxis()->SetTitle("Navg CX PH");
        if (verbose) g_Navg_CX_PH_VGain->Write(Form("Navg CX PH vs VGain_M%i_Ch%i_Bias%i", module, channel, bias));
        else         g_Navg_CX_PH_VGain->Write("Navg_CX_PH_VGain");

        TGraphErrors* g_Navg_PH_VGain = new TGraphErrors(vgains.size(), &vgains[0], &navg_phs[0],
                                                    &err_zeros[0], &err_zeros[0]);
        if (verbose) g_Navg_PH_VGain->SetTitle(Form("Navg PH vs VGain - M%i-%s Ch.%i Bias= %.2f OV= %.2f ", module, sipm.Data(), channel, bias_volt, overvoltage));
        else         g_Navg_PH_VGain->SetTitle(Form("Bias= %.2f OV= %.2f ", bias_volt, overvoltage));
        g_Navg_PH_VGain->GetXaxis()->SetTitle("VGain");
        g_Navg_PH_VGain->GetYaxis()->SetTitle("Navg PH");
        if (verbose) g_Navg_PH_VGain->Write(Form("Navg_PH_VGain_M%i_Ch%i_Bias%i", module, channel, bias));
        else         g_Navg_PH_VGain->Write("Navg_PH_VGain");

        TGraphErrors* g_Navg_PE_VGain = new TGraphErrors(vgains.size(), &vgains[0], &navg_pes[0],
                                                    &err_zeros[0], &err_zeros[0]);
        if (verbose) g_Navg_PE_VGain->SetTitle(Form("Navg PE vs VGain - M%i-%s Ch.%i Bias= %.2f OV= %.2f ", module, sipm.Data(), channel, bias_volt, overvoltage));
        g_Navg_PE_VGain->GetXaxis()->SetTitle("VGain");
        g_Navg_PE_VGain->GetYaxis()->SetTitle("Navg PE");
        if (verbose) g_Navg_PE_VGain->Write(Form("Navg_PE_VGain_M%i_Ch%i_Bias%i", module, channel, bias));
        else         g_Navg_PE_VGain->Write("Navg_PE_VGain");

        TGraphErrors* g_SNR_DR = new TGraphErrors(drs.size(), &drs[0], &snrs[0],
                                                  &err_zeros[0], &err_snrs[0]);
        if (verbose) g_SNR_DR->SetTitle(Form("SNR vs DR - M%i-%s Ch.%i Bias= %.2f OV= %.2f ", module, sipm.Data(), channel, bias_volt, overvoltage));
        else         g_SNR_DR->SetTitle(Form("Bias= %.2f OV= %.2f ", bias_volt, overvoltage));
        g_SNR_DR->GetXaxis()->SetTitle("Dynamic Range [pe]");
        g_SNR_DR->GetYaxis()->SetTitle("SNR");
        if (verbose) g_SNR_DR->Write(Form("SNR_DR_M%i_Ch%i_Bias%i", module, channel, bias));
        else         g_SNR_DR->Write("SNR_DR");

        TGraphErrors* g_SNR_SPEampl = new TGraphErrors(spe_ampls.size(), &spe_ampls[0], &snrs[0],
                                                      &err_zeros[0], &err_snrs[0]);
        if (verbose) g_SNR_SPEampl->SetTitle(Form("SNR vs SPE Amplitude - M%i-%s Ch.%i Bias= %.2f OV= %.2f ", module, sipm.Data(), channel, bias_volt, overvoltage));
        else         g_SNR_SPEampl->SetTitle(Form("Bias= %.2f OV= %.2f ", bias_volt, overvoltage));
        g_SNR_SPEampl->GetXaxis()->SetTitle("SPE Amplitude [ADC]");
        g_SNR_SPEampl->GetYaxis()->SetTitle("SNR");
        if (verbose) g_SNR_SPEampl->Write(Form("SNR_SPEampl_M%i_Ch%i_Bias%i", module, channel, bias));
        else         g_SNR_SPEampl->Write("SNR_SPEampl");
      
      } // end loop over biases
    } // end loop over channels
  } // end loop over modules

  std::cout << "Close file" << std::endl;
  file->Close();
}
