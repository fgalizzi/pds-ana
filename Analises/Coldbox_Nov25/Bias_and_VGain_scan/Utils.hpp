#include "../../../Class/_c/class_include.hpp"
#include <nlohmann/json.hpp>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>


using namespace std;
using json = nlohmann::json;

// --- ANA CONFIG ---------------------------------------------------
struct AnaConfig {
  bool display;
  bool plot;
  bool print;
  double tick_len;
  int nbins;
  string runs_folder;
  string data_format;
  string input_ana_folder;
  string rms_result_file;
  string output_ana_folder;
  double allowed_bsl_rms;
  bool   print_results;
  string bias_fit_csv;
  // int memorydepth;
};

inline AnaConfig load_ana_config(const std::string& filename) {
  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error("Could not open config file: " + filename);
  }

  json j;
  file >> j;
  AnaConfig config;
  config.runs_folder = j.at("runs_folder").get<std::string>();
  config.input_ana_folder = j.at("input_ana_folder").get<std::string>();
  config.rms_result_file = j.at("rms_result_file").get<std::string>();
  config.output_ana_folder = j.at("output_ana_folder").get<std::string>();
  config.allowed_bsl_rms = j.at("allowed_bsl_rms").get<double>();
  config.print_results = j.at("print_results").get<bool>();
  config.display = j.at("display").get<int>();
  config.print = j.at("print").get<int>();
  config.plot = j.at("plot").get<int>();
  config.tick_len = j.at("tick_len").get<double>();
  config.nbins = j.at("nbins").get<int>();
  config.data_format = j.at("data_format").get<std::string>();
  config.bias_fit_csv = j.at("bias_fit_csv").get<std::string>();

  return config;
}

// --- MODULE CONFIG ------------------------------------------------
struct ModuleConfig {
  int memorydepth;
  bool invert;
  int n_wf;
  int prepulse_ticks;
  int int_low;
  int int_up;
  int nmaxpeaks;
  int module;
  vector<int> module_channels;
  vector<double> daphne_biases;
  vector<int> vgains;
  string modules_in_foldername;
  string led_afe_extension;
  string custom_vgain_folder = "";
};

inline ModuleConfig load_module_config(const std::string& filename) {
  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error("Could not open config file: " + filename);
  }
  
  json j;
  file >> j;
  
  ModuleConfig config;
  config.memorydepth = j.at("memorydepth").get<int>();
  config.invert = j.at("invert").get<bool>();
  config.n_wf = j.at("n_wf").get<int>();
  config.prepulse_ticks = j.at("prepulse_ticks").get<int>();
  config.int_low = j.at("int_low").get<int>();
  config.int_up = j.at("int_up").get<int>();
  config.nmaxpeaks = j.at("nmaxpeaks").get<int>();
  config.module = j.at("module").get<int>();
  config.module_channels = j.at("module_channels").get<vector<int>>();
  config.daphne_biases = j.at("daphne_biases").get<vector<double>>();
  config.vgains = j.at("vgains").get<vector<int>>();
  config.modules_in_foldername = j.at("modules_in_foldername").get<string>();
  config.led_afe_extension = j.at("led_afe_extension").get<string>();
  if (j.contains("custom_vgain_folder")) {
    config.custom_vgain_folder = j.at("custom_vgain_folder").get<string>();
  }
  
  return config;
}

// --- BIAS-CONVERSION ------------------------------------------------
struct FitParams {
  double slope;
  double offset;
};


using FitKey = std::pair<int, std::string>;

struct FitKeyHash {
  std::size_t operator()(const FitKey& k) const {
    return std::hash<int>()(k.first) ^ std::hash<std::string>()(k.second);
  }
};

inline std::unordered_map<FitKey, FitParams, FitKeyHash>
read_bias_fit_csv(const std::string& filename){
  std::unordered_map<FitKey, FitParams, FitKeyHash> table;

  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error("Cannot open CSV file: " + filename);
  }

  std::string line;
  std::getline(file, line);  // skip header

  while (std::getline(file, line)) {
    std::stringstream ss(line);

    std::string afe_str, fit, slope_str, offset_str;

    std::getline(ss, afe_str, ',');
    std::getline(ss, fit, ',');
    std::getline(ss, slope_str, ',');
    std::getline(ss, offset_str, ',');

    int afe = std::stoi(afe_str);
    double slope = std::stod(slope_str);
    double offset = std::stod(offset_str);

    table[{afe, fit}] = {slope, offset};
  }

  return table;
}

inline double real_bias(int AFE,
                        const std::string& fit_type,
                        double bias,
                        const std::unordered_map<FitKey, FitParams, FitKeyHash>& table){

  auto it = table.find({AFE, fit_type});
  if (it == table.end()) {
    throw std::runtime_error(
      "No fit parameters for AFE=" + std::to_string(AFE) +
      ", Fit=" + fit_type
    );
  }

  const FitParams& p = it->second;
  return p.slope * bias + p.offset;
}
