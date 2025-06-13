#include <iostream>
#include <fstream>
#include <TString.h>
#include <TError.h>
#include <TBufferJSON.h>

#include "BestAnalysisConfiguration.hpp"

int main(int argc, char* argv[])
{
  if (argc > 2 or argc < 2) {
    printf("FATAL:Wrong number of arguments. Use TestBestAnalysisConfiguration read/write\n");
    return 0;
  }
  if (std::string(argv[1]) == "write") {
    std::ofstream outf("configuration.json");
    BestAnalysisConfiguration* conf = new BestAnalysisConfiguration();

    /* let's produce the configuration string */
    TString first_json = TBufferJSON::ToJSON(conf);
    /* write it */
    outf << first_json << std::endl;
    outf.close();
  } else if (std::string(argv[1]) == "read") {
    /* read the json configuration string */
    std::ifstream inf("configuration.json");
    if (inf.is_open()) {
      TString reco_json;
      std::string line;
      while (getline(inf, line)) {
        reco_json += line;
      }
      inf.close();
      /* let's produce the configuration object out of the string */
      BestAnalysisConfiguration* confreco = nullptr;
      TBufferJSON::FromJSON(confreco, reco_json);
      /* show it to chek it */
      TString test_json = TBufferJSON::ToJSON(confreco);
      std::cout << test_json << std::endl;
    } else {
      printf("ERROR:File not found\n");
    }
  } else {
    printf("FATAL:Wrong arguments. Use TestBestAnalysisConfiguration read/write\n");
  }
  return 1;
}
