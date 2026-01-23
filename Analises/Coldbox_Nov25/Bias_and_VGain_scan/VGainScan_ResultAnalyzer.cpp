#include "../../../Class/_c/class_include.hpp"
#include "../../../plotter/DUNEStyle.h"
#include "TGraphErrors.h"
#include <cstddef>
#include <string>
#include <vector>
using namespace std;

// --- HARDCODED DECLARATIONS ---------------------------------------
double ce_snr_edge = 25.; // From SNR vs SPE ampl. plots we observed
                          // that for SPE > 25 ADC the CE noise dominates


// --- FUNCTIONS ----------------------------------------------------
vector<double> unique_values_in_column(const vector<pair<string, vector<double>>>& data, size_t column_index){
  set<double> unique_values_set;
  for (const auto& value : data[column_index].second){
    unique_values_set.insert(value);
  }
  return vector<double>(unique_values_set.begin(), unique_values_set.end());
}

vector<double> get_channel_data(const vector<pair<string, vector<double>>>& data, size_t column_index, double channel){
  vector<double> channel_data;
  for (size_t i=0; i<data[0].second.size(); i++){
    if (data[1].second[i] == channel){
      channel_data.push_back(data[column_index].second[i]);
    }
  }
  return channel_data;
}

vector<double> get_channel_data(const vector<pair<string, vector<double>>>& data, size_t column_index, double channel,
                                size_t filter_columnn, double filter_value, string condition = ">"){
  vector<double> channel_data;
  for (size_t i=0; i<data[0].second.size(); i++){
    if (data[1].second[i] == channel){
      if (condition == ">" && data[filter_columnn].second[i] <= filter_value) continue;
      channel_data.push_back(data[column_index].second[i]);
    }
  }
  return channel_data;
}

TGraphErrors* form_tgraph(string title, string name, string x_title, string y_title,
                          const vector<pair<string, vector<double>>>& data,
                          double channel, size_t x_colunm, size_t y_colunm,
                          size_t x_err_colunm = SIZE_MAX, size_t y_err_colunm = SIZE_MAX){

  vector<double> x_values = get_channel_data(data, x_colunm, channel);
  vector<double> y_values = get_channel_data(data, y_colunm, channel);
  vector<double> x_errors(x_values.size(), 0.);
  vector<double> y_errors(y_values.size(), 0.);

  if (x_err_colunm != SIZE_MAX){
    x_errors = get_channel_data(data, x_err_colunm, channel);
  }
  if (y_err_colunm != SIZE_MAX){
    y_errors = get_channel_data(data, y_err_colunm, channel);
  }

  TGraphErrors* graph = new TGraphErrors(x_values.size(),
                                         &x_values[0], &y_values[0],
                                         &x_errors[0], &y_errors[0]);
  graph->SetTitle(title.c_str());
  graph->SetName(name.c_str());
  graph->GetXaxis()->SetTitle(x_title.c_str());
  graph->GetYaxis()->SetTitle(y_title.c_str());

  return graph;
}

// Function to compute the weighted average and its error given two vectors
// of values and their errors. Returns a pair (average, error).
pair<double, double> weighted_average(const vector<double>& values, const vector<double>& errors) {
  if (values.size() != errors.size() || values.empty()) {
    return make_pair(0.0, 0.0);
  }
  double sum_weighted_values = 0.0;
  double sum_weights = 0.0;

  for (size_t i = 0; i < values.size(); ++i) {
    if (errors[i] == 0) continue; // Avoid division by zero
    double weight = 1.0 / (errors[i] * errors[i]);
    sum_weighted_values += values[i] * weight;
    sum_weights += weight;
  }

  double average = sum_weighted_values / sum_weights;
  double error = sqrt(1.0 / sum_weights);

  return make_pair(average, error);
}

pair<double, double> weighted_average_from_data(const vector<pair<string, vector<double>>>& data,
                                               size_t value_colunm, size_t error_colunm,
                                               double channel) {
  vector<double> values = get_channel_data(data, value_colunm, channel);
  vector<double> errors = get_channel_data(data, error_colunm, channel);
  return weighted_average(values, errors);
}

pair<double, double> weighted_average_from_data(const vector<pair<string, vector<double>>>& data,
                                               size_t value_colunm, size_t error_colunm,
                                               double channel, size_t filter_columnn, double filter_value, string condition = ">") {
  vector<double> values = get_channel_data(data, value_colunm, channel, filter_columnn, filter_value, condition);
  vector<double> errors = get_channel_data(data, error_colunm, channel, filter_columnn, filter_value, condition);
  return weighted_average(values, errors);
}

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////
vector<int> modules = {1, 2, 3, 4, 5, 6};
// --- INPUT -----------------------------------------------------
string ana_folder = "/eos/home-f/fegalizz/ColdBox_VD/November25/SpyBuffer/VGain_Scans/allowed_bsl_rms_3/";

size_t channel_colunm         = 1;
size_t bias_volt_colunm       = 3;
size_t overvoltage_colunm     = 4;
size_t vgain_colunm           = 5;
size_t gain_colunm            = 11;
size_t err_gain_colunm        = 12;
size_t spe_ampl_colunm        = 13;
size_t dr_colunm              = 14;
size_t snr_colunm             = 15;
size_t err_snr_colunm         = 16;
size_t cx_colunm              = 18;
size_t err_cx_colunm          = 19;
size_t navg_cx_phs_colunm     = 20;
size_t err_navg_cx_phs_colunm = 21;

///////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void VGainScan_ResultAnalyzer(){
  dunestyle::SetDuneStyle();
  auto style = gROOT->GetStyle("duneStyle");
  style->SetOptFit(1111); style->SetOptTitle(0);
  style->SetStatX(0.9);   style->SetStatY(0.9);
  style->SetStatBorderSize(0);

  // --- LOOP OVER MODULES ----------------------------------------
  for (const auto& module : modules){

    // Create and open a file.root where to store the canvas
    TFile* out_file = new TFile((ana_folder+"M"+to_string(module)+"_VGainScans_ResultSummary.root").c_str(), "RECREATE");
    out_file->cd();

    string module_folder = ana_folder + Form("M%i/", module);
    // Store all the csv files in the module folder in a vector
    vector<string> csv_files;
    for (const auto& entry : filesystem::directory_iterator(module_folder)){
      string file_path = entry.path().string();
      if (file_path.find(".csv") != string::npos){
        csv_files.push_back(file_path);
      }
    }

    vector<TGraphErrors*> g_CX_OV_vec = {};
    vector<TGraphErrors*> g_CESNR_OV_vec = {};
    for (size_t i = 1; i <= 2; i++) {
      g_CX_OV_vec.push_back(new TGraphErrors());
      g_CX_OV_vec.back()->SetTitle(Form("M%i_CX_OV_Ch%i", module, int(i)));
      g_CX_OV_vec.back()->SetName(Form("M%i_CX_OV_Ch%i", module, int(i)));
      g_CX_OV_vec.back()->GetXaxis()->SetTitle("Overvoltage [V]");
      g_CX_OV_vec.back()->GetYaxis()->SetTitle("Cross Talk [%]");

      g_CESNR_OV_vec.push_back(new TGraphErrors());
      g_CESNR_OV_vec.back()->SetTitle(Form("M%i_CESNR_OV_Ch%i", module, int(i)));
      g_CESNR_OV_vec.back()->SetName(Form("M%i_CESNR_OV_Ch%i", module, int(i)));
      g_CESNR_OV_vec.back()->GetXaxis()->SetTitle("Overvoltage [V]");
      g_CESNR_OV_vec.back()->GetYaxis()->SetTitle("CE SNR");
    }

    // --- LOOP OVER CSV FILES ------------------------------------
    for (const auto& file : csv_files){
      // remove ".csv" from the file name
      std::cout << file << std::endl;
      string dir_name = file.substr(file.find_last_of("/")+1);
      dir_name = dir_name.substr(0, dir_name.find(".csv"));
      std::cout << dir_name << std::endl;

      TDirectory* dir_setting = out_file->mkdir(dir_name.c_str());
      dir_setting->cd();

      vector<pair<string, vector<double>>> input_data = read_vec_pair_CSV(file.c_str());
      vector<double> channels = unique_values_in_column(input_data, channel_colunm);

      // --- LOOP OVER CHANNELS -----------------------------------
      for (const auto& channel : channels){
        TDirectory* dir_ch = dir_setting->mkdir(Form("ch_%i", int(channel)));
        dir_ch->cd();

        double bias_volt, overvoltage;
        bias_volt = input_data[bias_volt_colunm].second[0];
        overvoltage = input_data[overvoltage_colunm].second[0];

        int ch_index = distance(channels.begin(),
                                find(channels.begin(), channels.end(), channel));

        // Fill CX vs OV graph
        pair<double, double> cx_avg = weighted_average_from_data(input_data,
                                                                 cx_colunm, err_cx_colunm,
                                                                 channel);
        g_CX_OV_vec[ch_index]->SetPoint(g_CX_OV_vec[ch_index]->GetN(),
                                        overvoltage, cx_avg.first);
        g_CX_OV_vec[ch_index]->SetPointError(g_CX_OV_vec[ch_index]->GetN()-1,
                                             0., cx_avg.second);

        pair<double, double> snr_avg = weighted_average_from_data(input_data,
                                                       snr_colunm, err_snr_colunm,
                                                       channel,
                                                       spe_ampl_colunm, ce_snr_edge, ">");
        if (snr_avg.first != 0.) {
          g_CESNR_OV_vec[ch_index]->SetPoint(g_CESNR_OV_vec[ch_index]->GetN(),
                                             overvoltage, snr_avg.first);
          g_CESNR_OV_vec[ch_index]->SetPointError(g_CESNR_OV_vec[ch_index]->GetN()-1,
                                                  0., snr_avg.second);
        }

        // --- TGRAPHS --------------------------------------------
        string graph_title = Form("M%i_Ch%i_BiasVolt_%.2f_OV_%.2f", module, int(channel), bias_volt, overvoltage);
        TGraphErrors* g_Gain_VGain = form_tgraph(graph_title, "Gain_VGain", "VGain", "Gain [ADC#times ticks]",
                                                 input_data, channel,
                                                 vgain_colunm, gain_colunm,
                                                 SIZE_MAX, err_gain_colunm);

        TGraphErrors* g_SPEampl_VGain = form_tgraph(graph_title, "SPEampl_VGain", "VGain", "SPE Amplitude [ADC]",
                                                    input_data, channel,
                                                    vgain_colunm, spe_ampl_colunm);

        TGraphErrors* g_DR_VGain = form_tgraph(graph_title, "DR_VGain", "VGain", "Dynamic Range [pe]",
                                               input_data, channel,
                                               vgain_colunm, dr_colunm);

        TGraphErrors* g_SNR_VGain = form_tgraph(graph_title, "SNR_VGain", "VGain", "SNR",
                                                input_data, channel,
                                                vgain_colunm, snr_colunm,
                                                SIZE_MAX, err_snr_colunm);

        TGraphErrors* g_CX_VGain = form_tgraph(graph_title, "CX_VGain", "VGain", "Cross Talk [%]",
                                               input_data, channel,
                                               vgain_colunm, cx_colunm,
                                               SIZE_MAX, err_cx_colunm);

        TGraphErrors* g_Navg_CX_PH_VGain = form_tgraph(graph_title, "Navg_CX_PH_VGain", "VGain", "Navg CX PH",
                                                       input_data, channel,
                                                       vgain_colunm, navg_cx_phs_colunm,
                                                       SIZE_MAX, err_navg_cx_phs_colunm);

        TGraphErrors* g_SNR_DR = form_tgraph(graph_title, "SNR_DR", "Dynamic Range [pe]", "SNR",
                                             input_data, channel,
                                             dr_colunm, snr_colunm,
                                             SIZE_MAX, err_snr_colunm);

        TGraphErrors* g_SNR_SPEampl = form_tgraph(graph_title, "SNR_SPEampl", "SPE Amplitude [ADC]", "SNR",
                                                  input_data, channel,
                                                  spe_ampl_colunm, snr_colunm,
                                                  SIZE_MAX, err_snr_colunm);

        // --- WRITE TGRAPHS --------------------------------------
        g_Gain_VGain->Write();
        g_SPEampl_VGain->Write();
        g_DR_VGain->Write();
        g_SNR_VGain->Write();
        g_CX_VGain->Write();
        g_Navg_CX_PH_VGain->Write();
        g_SNR_DR->Write();
        g_SNR_SPEampl->Write();
      } // end loop over biases
    } // end loop over channels
    
    out_file->cd();
    for (auto g_CX_OV : g_CX_OV_vec) g_CX_OV->Write();
    for (auto g_CESNR_OV : g_CESNR_OV_vec) g_CESNR_OV->Write();

    out_file->Close();
  }
}
