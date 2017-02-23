//
//  SimulationResult.cpp
//  fits
//
//  Created by Tal Zinger on 15/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "SimulationResult.hpp"


/*
 std::string sim_id;
 FLOAT_TYPE distance_from_actual;
 std::vector<FLOAT_TYPE> individual_distance;
 std::vector<FLOAT_TYPE> fitness_values;
 
 MATRIX_TYPE mutation_rates;
 int N;
 
 int wt_index;
 
 // need to save memory, reaching over 13GB 2016-04-20
 // reduced number of simulation so checking this again
 std::string raw_data;
 */


SimulationResult::SimulationResult()
: sim_id(""),
distance_from_actual(-1.0),
fitness_values(),
raw_data(),
wt_index(-1),
individual_distance(),
N(0),
mutation_rates(),
prior_sample_index(0),
sim_data_matrix()
{}


SimulationResult::SimulationResult(std::string id, FLOAT_TYPE distance)
: sim_id(id),
distance_from_actual(distance),
fitness_values(),
raw_data(),
wt_index(-1),
individual_distance(),
N(0),
mutation_rates(),
prior_sample_index(0),
sim_data_matrix()
{}


SimulationResult::SimulationResult(const SimulationResult &original)
: sim_id(original.sim_id),
distance_from_actual(original.distance_from_actual),
fitness_values(original.fitness_values),
raw_data(original.raw_data),
wt_index(original.wt_index),
individual_distance(original.individual_distance),
N(original.N),
mutation_rates(original.mutation_rates),
prior_sample_index(original.prior_sample_index),
sim_data_matrix(original.sim_data_matrix)
{}


SimulationResult::SimulationResult(const CMulator& sim_object)
: sim_id(sim_object.GetSimUID()),
distance_from_actual(-1.0f),
fitness_values(sim_object.GetAlleleFitnessValues()),
raw_data(sim_object.GetAllOutputAsText()),
wt_index(sim_object.GetWTAllele()),
individual_distance(sim_object.GetAlleleNumber(), -1.0f),
N(sim_object.GetPopulationSize()),
mutation_rates(sim_object.GetMutationRateMatrix()),
prior_sample_index(0),
sim_data_matrix(sim_object.GetAllOutputAsMatrix())
{
}


// move constructor
// creates a minimal object that will be destructed with minimal effort
SimulationResult::SimulationResult( SimulationResult&& other ) noexcept
{
    // simple data types
    std::swap(N, other.N);
    std::swap(wt_index, other.wt_index);
    std::swap(distance_from_actual, other.distance_from_actual);
    std::swap(prior_sample_index, other.prior_sample_index);
    
    // complex, use swap
    sim_id.swap(other.sim_id);
    raw_data.swap(other.raw_data);
    fitness_values.swap(other.fitness_values);
    individual_distance.swap(other.individual_distance);
    mutation_rates.swap(other.mutation_rates);
    sim_data_matrix.swap(other.sim_data_matrix);
}


// Operators
// the less-than is the only one used for comparisons
bool SimulationResult::operator<(const SimulationResult& result) const
{
    
    // 2016-10-25 try to take individual allleles into account
    
    // old version
    /*
     if ( individual_distance.empty() || result.individual_distance.empty()) {
     std::cout << "OLD Comparison" << ( individual_distance.empty() ? "1" : " " )
     << ( result.individual_distance.empty() ? "2" : " " ) << std::endl;
     return distance_from_actual < result.distance_from_actual;
     }
     
     if ( individual_distance.size() != result.individual_distance.size() ) {
     return distance_from_actual < result.distance_from_actual;
     }
     
     // skip calculations if possible
     if ( distance_from_actual >= result.distance_from_actual ) {
     return false;
     }
     
     // more strict that simply sum of distances
     auto all_alleles_shorter_distance = true;
     for ( auto i=0; i<individual_distance.size(); ++i ) {
     all_alleles_shorter_distance = all_alleles_shorter_distance &&
     ( individual_distance[i] <= result.individual_distance[i] );
     }
     
     return all_alleles_shorter_distance;
     */
    
    return distance_from_actual < result.distance_from_actual;
}

/*
 bool CMulatorToolbox::SimulationResult::operator>(const SimulationResult& result) const
 {
	return distance_from_actual > result.distance_from_actual;
 }
 
 bool CMulatorToolbox::SimulationResult::operator==(const SimulationResult& result) const
 {
	return distance_from_actual == result.distance_from_actual;
 }
 */

// Swap functions
void SimulationResult::swap( SimulationResult& other )
{
    //auto tmp_distance = distance_from_actual;
    //distance_from_actual = other.distance_from_actual;
    //other.distance_from_actual = tmp_distance;
    
    
    // simple data
    std::swap(N, other.N);
    std::swap(wt_index, other.wt_index);
    std::swap(distance_from_actual, other.distance_from_actual);
    std::swap(prior_sample_index, other.prior_sample_index);
    
    // complex data types
    sim_id.swap(other.sim_id);
    fitness_values.swap(other.fitness_values);
    raw_data.swap(other.raw_data);
    individual_distance.swap(other.individual_distance);
    mutation_rates.swap(other.mutation_rates);
    sim_data_matrix.swap(other.sim_data_matrix);
}


void SimulationResult::swap(SimulationResult& res1, SimulationResult& res2)
{
    std::swap(res1.N, res2.N);
    std::swap(res1.wt_index, res2.wt_index);
    std::swap(res1.distance_from_actual, res2.distance_from_actual);
    std::swap(res1.prior_sample_index, res2.prior_sample_index );
    
    res1.sim_id.swap(res2.sim_id);
    res1.raw_data.swap(res2.raw_data);
    res1.fitness_values.swap(res2.fitness_values);
    res1.individual_distance.swap(res2.individual_distance);
    res1.sim_data_matrix.swap(res2.sim_data_matrix);
}


void SimulationResult::swap(SimulationResult *res1, SimulationResult *res2)
{
    std::swap(res1->N, res2->N);
    std::swap(res1->wt_index, res2->wt_index);
    std::swap(res1->distance_from_actual, res2->distance_from_actual);
    std::swap(res1->prior_sample_index, res2->prior_sample_index );
    
    res1->sim_id.swap(res2->sim_id);
    res1->raw_data.swap(res2->raw_data);
    res1->fitness_values.swap(res2->fitness_values);
    res1->individual_distance.swap(res2->individual_distance);
    res1->sim_data_matrix.swap(res2->sim_data_matrix);
}


// hack I learned that unifies assignment for copy and move
// other is passed by value, so it's a copy
// swapping will essentialy move the data, leaving other with garbage
// either way, "other" will be destructed, leaving us with the data
// and the original rvalue is untouched
SimulationResult& SimulationResult::operator=(SimulationResult other)
{
    swap(other);
    return *this;
}

// TODO: move this to result stats - it needs to go through all 100 best results or something
std::vector<FLOAT_TYPE> SimulationResult::GetSDForEachAllele()
{
    
    std::vector<FLOAT_TYPE> allele_sd_vec(sim_data_matrix.size2(), 0.0f);
    
    for ( auto current_allele=0; current_allele<sim_data_matrix.size2(); ++current_allele ) {
        
        boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col( sim_data_matrix, current_allele );
        
        boost::accumulators::accumulator_set<
        FLOAT_TYPE,
        boost::accumulators::stats<
        boost::accumulators::tag::median,
        boost::accumulators::tag::variance,
        boost::accumulators::tag::mean,
        boost::accumulators::tag::min,
        boost::accumulators::tag::max> > allele_accumulator;
        
        //std::cout << "Accumulating allele " << current_allele << std::endl;
        for ( auto current_freq : current_col ) {
            allele_accumulator( current_freq );
            std::cout << current_freq << std::endl;
        }
        
        //current_col = current_col * 2.0f;
        
        allele_sd_vec[current_allele] = std::sqrt( boost::accumulators::variance(allele_accumulator) );
        
        //std::cout << " sd is " << allele_sd_vec[current_allele] << std::endl;
    }
    
    
    return allele_sd_vec;
}

void SimulationResult::DivideEachAllele( std::vector<FLOAT_TYPE> value_vector )
{
    for ( auto current_allele=0; current_allele<sim_data_matrix.size2(); ++current_allele ) {
        
        boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col( sim_data_matrix, current_allele );
        current_col = current_col / value_vector[current_allele];
    }
    
}
