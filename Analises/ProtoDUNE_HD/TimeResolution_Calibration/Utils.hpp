#include "../../../Class/_c/class_include.hpp"
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;

// --- ANA CONFIG ---------------------------------------------------
struct AnaConfig {
  string runs_folder;
  string input_ana_folder;
  string rms_result_file;
  string channel_map_file;
  string output_ana_folder;
  double allowed_bsl_rms;
  bool   print_results;
  bool display;
  bool print;
  bool plot;
  int memorydepth;
  double tick_len;
  int nbins;
  int nmaxpeaks;
  string data_format;
};

AnaConfig load_ana_config(const std::string& filename) {
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
  config.channel_map_file = j.at("channel_map_file").get<std::string>();
  config.output_ana_folder = j.at("output_ana_folder").get<std::string>();
  config.allowed_bsl_rms = j.at("allowed_bsl_rms").get<double>();
  config.print_results = j.at("print_results").get<bool>();
  config.display = j.at("display").get<bool>();
  config.print = j.at("print").get<bool>();
  config.plot = j.at("plot").get<bool>();
  config.memorydepth = j.at("memorydepth").get<int>();
  config.tick_len = j.at("tick_len").get<double>();
  config.nbins = j.at("nbins").get<int>();
  config.nmaxpeaks = j.at("nmaxpeaks").get<int>();
  config.data_format = j.at("data_format").get<std::string>();
  return config;
}

// --- MODULE CONFIG ------------------------------------------------
struct ModuleConfig {
  bool invert;
  int n_wf;
  int prepulse_ticks;
  int int_low;
  int int_up;
  std::vector<int> channels;
  std::vector<int> runs;
};

ModuleConfig load_module_config(const std::string& filename) {
  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error("Could not open config file: " + filename);
  }
  
  json j;
  file >> j;
  
  ModuleConfig config;
  config.invert = j.at("invert").get<bool>();
  config.n_wf = j.at("n_wf").get<int>();
  config.prepulse_ticks = j.at("prepulse_ticks").get<int>();
  config.int_low = j.at("int_low").get<int>();
  config.int_up = j.at("int_up").get<int>();

  config.channels = j.at("channels").get<std::vector<int>>();
  config.runs = j.at("runs").get<std::vector<int>>();
  
  return config;
}
