//
//  SimulationResult.hpp
//  fits
//
//  Created by Tal Zinger on 15/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#ifndef SimulationResult_hpp
#define SimulationResult_hpp

#include <vector>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/accumulators/numeric/functional.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/framework/features.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/format.hpp>

#include "CMulator.h"

struct SimulationResult {
    std::string sim_id;
    FLOAT_TYPE distance_from_actual;
    std::vector<FLOAT_TYPE> individual_distance;
    std::vector<FLOAT_TYPE> fitness_values;
    
    MATRIX_TYPE mutation_rates;
    int N;
    
    int wt_index;
    
    // in order to uniquely identify the sample from the prior used to simulate this result
    std::size_t prior_sample_index;
    
    // need to save memory, reaching over 13GB 2016-04-20
    // reduced number of simulation so checking this again
    std::string raw_data;
    
    MATRIX_TYPE sim_data_matrix;
    
    
    // constructors
    SimulationResult();
    SimulationResult(std::string id, FLOAT_TYPE distance);
    SimulationResult(const SimulationResult &original); // copy constructor
    SimulationResult(SimulationResult&& other) noexcept; // move constructor
    SimulationResult(const CMulator& sim_object);
    
    // swap for sorting
    void swap(SimulationResult& other);
    void swap(SimulationResult& res1, SimulationResult& res2);
    void swap(SimulationResult *res1, SimulationResult *res2);
    
    // for sorting so the rest are redundant
    bool operator<(const SimulationResult& result) const;
    
    SimulationResult& operator=(SimulationResult other);
    
    // this is used when the simulation result is taken as pseudo-data
    std::vector<FLOAT_TYPE> GetSDForEachAllele();
    void DivideEachAllele( std::vector<FLOAT_TYPE> value_vector );
};

#endif /* SimulationResult_hpp */
