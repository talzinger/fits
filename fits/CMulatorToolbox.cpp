//
//  CMulatorToolbox.cpp
//  CMulatorLib
//
//  Created by Tal Zinger on 03/01/2016.
//  Copyright © 2016 Tal Zinger. All rights reserved.
//

// #define DEBUG_VERBOSE

#include "CMulatorToolbox.hpp"


// distance function will calculate total distance,
// one far allele and a few close ones could be similar to
// few medicore-distance alleles.
// this kernel would allow to shrink distance for the latter
// TODO: REPLACE WITH A REAL KERNEL FUNCTION
// https://en.wikipedia.org/wiki/Kernel_(statistics)
//
FLOAT_TYPE CalculateDistKernel( std::vector<FLOAT_TYPE> actual_data, std::vector<FLOAT_TYPE> sim_data )
{
    // divide them between each other
    
    std::vector<FLOAT_TYPE> ratios(actual_data.size(), 0.0f);
    
    std::transform( sim_data.cbegin(), sim_data.cend(), actual_data.cbegin(), ratios.begin(),
                   []( FLOAT_TYPE num1, FLOAT_TYPE num2 ) {
                       return ( num1 / num2 );
                   } );
    
    // mean
    FLOAT_TYPE tmp_val = 0.0f;
    for ( auto val : ratios ) {
        tmp_val += val;
    }
    tmp_val = tmp_val / static_cast<FLOAT_TYPE>( ratios.size() );
    
    // closer to 1 is better, closer to 0 is worse
    // we want the better distance to be smaller than the worse
    // thus, we need to return the reciprocal
    tmp_val = 1.0f / tmp_val;
    
    return tmp_val;
}

FLOAT_TYPE CMulatorToolbox::FreqsSameMagnitude( FLOAT_TYPE original_val )
{
    while (original_val<0.1 && original_val>0.0) {
        original_val *= 10.0f;
    }
    
    return original_val;
}


std::vector<FLOAT_TYPE> CMulatorToolbox::GetDistanceVectorSimActual( std::vector<FLOAT_TYPE> actual_data, std::vector<FLOAT_TYPE> sim_data, std::size_t num_alleles )
{
    
    // sized don't match
    if ( actual_data.size() != sim_data.size() ) {
        throw "CMulatorBatch-GetDstanceSimActual: size mismatch. actual "
        + std::to_string(actual_data.size()) + ", sim " + std::to_string(sim_data.size());
    }
    
    // at least one of them is empty
    if ( actual_data.size() < 1 || sim_data.size() < 1 ) {
        throw "CMulatorBatch-GetDstanceSimActual: illegal size. actual "
        + std::to_string(actual_data.size()) + ", sim " + std::to_string(sim_data.size());
    }
    
    
    // sum for each allele. init with 0.
    std::vector<FLOAT_TYPE> allele_distance(num_alleles, 0);
    std::vector<FLOAT_TYPE> distances(actual_data.size(), 0.0f);
    
    std::transform( sim_data.cbegin(), sim_data.cend(), actual_data.cbegin(), distances.begin(),
                   []( FLOAT_TYPE num1, FLOAT_TYPE num2 ) {
                       return ( std::fabs(num1-num2) );
                   } );
    
    for ( auto i=0; i<distances.size(); i++ ) {
        
        auto current_allele = i % num_alleles;
        
        allele_distance[current_allele] += distances[i];
    }
    
    
    /*
     FLOAT_TYPE dist = 0.0f;
    for ( auto num : allele_distance ) {
    
        std::cout << "Allele " << num << std::endl;
        dist += num;
        
    }
    std::cout << "Total=" << dist << std::endl;
    */
    return allele_distance;
}


// This assumes generations are exactly the same in both
FLOAT_TYPE CMulatorToolbox::GetDistanceSimActual( MATRIX_TYPE actual_data, MATRIX_TYPE sim_data )
{
    MATRIX_TYPE diff_matrix = actual_data - sim_data;
    
    diff_matrix = boost::numeric::ublas::abs(diff_matrix);
    
    // trick - sum can be for vectors, not matrices. multiply matrix by unit vector to get a vector
    auto sum_diff = sum(prod(boost::numeric::ublas::scalar_vector<float>(diff_matrix.size1()), diff_matrix));
    
    return sum_diff;
}

FLOAT_TYPE CMulatorToolbox::GetDistanceSimActual( std::vector<FLOAT_TYPE> actual_data, std::vector<FLOAT_TYPE> sim_data )
{

	// sized don't match
	if ( actual_data.size() != sim_data.size() ) {
		throw "CMulatorBatch-GetDstanceSimActual: size mismatch. actual " 
		+ std::to_string(actual_data.size()) + ", sim " + std::to_string(sim_data.size());
	}
	
	// at least one of them is empty
	if ( actual_data.size() < 1 || sim_data.size() < 1 ) {
		throw "CMulatorBatch-GetDstanceSimActual: illegal size. actual " 
		+ std::to_string(actual_data.size()) + ", sim " + std::to_string(sim_data.size());
	}
	
	
    // TODO: combine distance and summation into a single command
    std::vector<FLOAT_TYPE> tmp_distance(0);
    tmp_distance.resize( actual_data.size() );
    
    
    std::transform( sim_data.cbegin(), sim_data.cend(), actual_data.cbegin(), tmp_distance.begin(),
                   []( FLOAT_TYPE num1, FLOAT_TYPE num2 ) {
                       // return( std::fabs( FreqsSameMagnitude(num1) - FreqsSameMagnitude(num2) ) ));
                       return ( std::fabs(num1-num2) );
                   } );
    
    FLOAT_TYPE tmpzero = 0.0;
    auto sum_distance = std::accumulate( tmp_distance.cbegin(), tmp_distance.cend(), tmpzero );
    
    return sum_distance;
}


/*
 FLOAT_TYPE CMulatorToolbox::GetDistanceSimActual( ActualDataEntry actual_data, SimulationResult sim_data )
{
	return 0.0;
}
*/


/*
FLOAT_TYPE CMulatorToolbox::OrderResultsAndUpdateThreshold(std::vector<SimulationResult> &result_vector, int max_idx_to_keep, bool erase_over_threshold)
{
	if (max_idx_to_keep > result_vector.size()) {
		std::cerr << "Error: result vector too small for sorting, using vector size as nth-element" << std::endl;
		throw "Error: result vector too small for sorting, using vector size as nth-element";
	}


	std::nth_element(result_vector.begin(), result_vector.begin() + max_idx_to_keep, result_vector.end());

	auto new_threshold_value = (result_vector.begin() + max_idx_to_keep)->distance_from_actual;

	// at least one item to erase
	if ( erase_over_threshold && (result_vector.size() > max_idx_to_keep)) {
		result_vector.erase(result_vector.begin() + max_idx_to_keep, result_vector.end());
	}
	
	return new_threshold_value;
}
*/


// TODO: use accumulators and save this
std::vector<FLOAT_TYPE> CMulatorToolbox::GetMinFitness(const std::vector<SimulationResult> &result_vector)
{
    std::vector<FLOAT_TYPE> tmp_min_vector;

	tmp_min_vector = result_vector[0].fitness_values;

	for ( auto res : result_vector ) {

		for ( auto i = 0; i < res.fitness_values.size(); ++i ) {

			if (tmp_min_vector[i] > res.fitness_values[i]) {

				tmp_min_vector[i] = res.fitness_values[i];
			}
		}
	}

	return tmp_min_vector;
}


// TODO: use accumulators and save this
std::vector<FLOAT_TYPE> CMulatorToolbox::GetMaxFitness(const std::vector<SimulationResult> &result_vector)
{
    std::vector<FLOAT_TYPE> tmp_max_vector;

	tmp_max_vector = result_vector[0].fitness_values;

	for ( auto res : result_vector ) {

		for ( auto i = 0; i < res.fitness_values.size(); ++i ) {

			if (tmp_max_vector[i] > res.fitness_values[i]) {

				tmp_max_vector[i] = res.fitness_values[i];
			}
		}
	}

	return tmp_max_vector;
}


std::vector<std::vector<FLOAT_TYPE>> CMulatorToolbox::BuildPosteriorFitness(const std::vector<SimulationResult>& result_vector)
{
    
    std::vector<std::vector<FLOAT_TYPE>> posterior_vector;
    
    for ( auto tmp_result : result_vector ) {
        
        posterior_vector.push_back( tmp_result.fitness_values );
    }
    
    return posterior_vector;
}

/*
ZParams CMulatorToolbox::ReadParametersFromFile( std::string param_filename )
{
    ZParams my_zparams;
    
    try {
        // parameters should be read only
        std::cout << "Reading parameters file: " << param_filename << std::endl;
        my_zparams.ReadParameters(param_filename, true);
    }
    catch( const char* txt ) {
        std::cerr << "Error in zparams: " << txt << std::endl;
    }
    catch (...) {
        std::cerr << "Error while reading parameter file: " << param_filename << std::endl;
    }
    
    if (!my_zparams.IsEmpty()) {
        std::cout << "Zparams:" << std::endl;
        std::cout << my_zparams.GetAllParameters() << std::endl;
    }
    
    //_initialized_with_parameters = true;
    
    //InitMemberVariables(my_zparams);
    
    return my_zparams;
}
*/
