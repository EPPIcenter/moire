
#include "Data.h"
// #include "misc_v1.h"

using namespace std;

//------------------------------------------------
// declare static member variables for class Data_biallelic

vector<vector<int> > Data_biallelic::data;
int Data_biallelic::n;
int Data_biallelic::L;

//------------------------------------------------
// constructor for Data_biallelic class
Data_biallelic::Data_biallelic(const Rcpp::List &args) {

  // data = rcpp_to_mat_int(args["data"]);
  // n = rcpp_to_int(args["n"]);
  // L = rcpp_to_int(args["L"]);
}

//------------------------------------------------
// declare static member variables for class Data_multiallelic

vector<vector<vector<int> > > Data_multiallelic::data;
vector<int> Data_multiallelic::observed_COI;
vector<int> Data_multiallelic::alleles;
int Data_multiallelic::n;
int Data_multiallelic::L;

//------------------------------------------------
// constructor for Data_multiallelic class
Data_multiallelic::Data_multiallelic(const Rcpp::List &args) {

  // read in data
  // data = rcpp_to_array_int(args["data"]);
  // alleles = rcpp_to_vector_int(args["alleles"]);
  // n = rcpp_to_int(args["n"]);
  // L = rcpp_to_int(args["L"]);

  // get observed number of alleles in each sample
  // observed_COI = vector<int>(n);
  // for (int i=0; i<n; ++i) {
  //   for (int l=0; l<L; ++l) {
  //     if (int(data[i][l].size()) > observed_COI[i]) {
  //       observed_COI[i] = data[i][l].size();
  //     }
  //   }
  // }

}
