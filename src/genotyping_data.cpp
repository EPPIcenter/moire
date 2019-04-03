
#include "genotyping_data.h"
#include "mcmc_utils.h"

//------------------------------------------------
// declare static member variables for class Data_multiallelic
std::vector<std::vector<std::vector<int> > > GenotypingData::observed_alleles;
std::vector<int> GenotypingData::observed_coi;
std::vector<int> GenotypingData::num_alleles;
size_t GenotypingData::num_samples;
size_t GenotypingData::num_loci;
int GenotypingData::max_alleles;

//------------------------------------------------
// constructor for Data_multiallelic class
GenotypingData::GenotypingData(const Rcpp::List &args) {

  // read in data
  observed_alleles = UtilFunctions::r_to_array_int(args["data"]);
  
  num_loci = observed_alleles.size();
  num_samples = observed_alleles[0].size();
  
  num_alleles = std::vector<int>(num_loci);
  observed_coi = std::vector<int>(num_samples, 0);

  for(size_t i = 0; i < num_loci; i++) {
    num_alleles[i] = observed_alleles[i][0].size();
    if (num_alleles[i] > max_alleles) {
      max_alleles = num_alleles[i];
    }
    for(size_t j = 0; j < num_samples; j++) {
      int total_alleles = 0;
      for(size_t k = 0; k < observed_alleles[i][j].size(); k++) {
        total_alleles += observed_alleles[i][j][k];
      }

      if(total_alleles > observed_coi[j]) {
        observed_coi[j] = total_alleles;
      }
    }
    
  }

}
