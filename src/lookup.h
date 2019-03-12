#pragma once

#ifndef LOOKUP_H_
#define LOOKUP_H_

#include <Rcpp.h>
#include "parameters.h"

class Lookup {
    private:
    static bool lgamma_initialized;
    
    public:
    Lookup() {};
    Lookup(int max_coi);
    static std::vector<double> lookup_lgamma;
    void init_lgamma(int max_coi);
};


#endif // LOOKUP_H_