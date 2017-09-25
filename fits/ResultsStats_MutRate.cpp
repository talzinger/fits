//
//  ResultsStats_MutRate.cpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "ResultsStats.hpp"

void ResultsStats::CalculateStatsMutation(const std::vector<SimulationResult>& result_vector)
{
    _num_alleles = result_vector[0].fitness_values.size();
    
    _num_results = result_vector.size();
    
    boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_distance;
    
    boost::numeric::ublas::matrix< boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > acc_matrix(_num_alleles,_num_alleles);
    
    boost::numeric::ublas::matrix< std::vector<FLOAT_TYPE> > median_matrix(_num_alleles,_num_alleles);
    
    min_mutation_rates.resize(_num_alleles, _num_alleles);
    max_mutation_rates.resize(_num_alleles, _num_alleles);
    mean_mutation_rates.resize(_num_alleles, _num_alleles);
    median_mutation_rates.resize(_num_alleles, _num_alleles);
    normalized_median_mutation_rates.resize(_num_alleles, _num_alleles);
    
    for ( auto sim_result : result_vector ) {
        
        acc_distance(sim_result.distance_from_actual);
        
        for ( auto row=0; row<_num_alleles; ++row ) {
            
            for ( auto col=0; col<_num_alleles; ++col ) {
                
                acc_matrix(row,col)( sim_result.mutation_rates(row,col) );
                median_matrix(row,col).push_back(sim_result.mutation_rates(row,col));
            }
        }
    }
    
    for ( auto row=0; row<_num_alleles; ++row ) {
        
        for ( auto col=0; col<_num_alleles; ++col ) {
            
            min_mutation_rates(row,col) = boost::accumulators::min( acc_matrix(row,col) );
            max_mutation_rates(row,col) = boost::accumulators::max( acc_matrix(row,col) );
            
            mean_mutation_rates(row,col) = boost::accumulators::mean( acc_matrix(row,col) );
            
            median_mutation_rates(row,col) = GetMedian(median_matrix(row,col));
        }
    }
    
    _distance_min = boost::accumulators::min(acc_distance);
    _distance_max = boost::accumulators::max(acc_distance);
    _distance_mean = boost::accumulators::mean(acc_distance);
    _distance_sd = std::sqrt(boost::accumulators::variance(acc_distance));
    
    for ( auto row=0; row<_num_alleles; ++row ) {
        
        float tmp_sum = 0.0f;
        
        for ( auto col=0; col<_num_alleles; ++col ) {
            
            tmp_sum += median_mutation_rates(row,col);
        }
        
        for ( auto col=0; col<_num_alleles; ++col ) {
            
            normalized_median_mutation_rates(row,col) = median_mutation_rates(row,col) / tmp_sum;
        }
    }
}


std::string ResultsStats::GetSummaryMutRate()
{
    std::stringstream ss;
    
    ss.imbue(std::locale(fits_constants::used_locale));
    
    ss << "Mutation Rate Report" << std::endl;
    ss << GetSummaryHeader();
    
    //auto current_time_raw = std::chrono::system_clock::now();
    //auto current_time = std::chrono::system_clock::to_time_t(current_time_raw);
    //auto current_time_final = *std::localtime(&current_time);
    //ss << std::put_time(&current_time_final, "%F %T") << std::endl;
    
    //std::cout << "Simulation results used for calculations: " << _num_results << std::endl;
    
    //ss << "=====================" << std::endl;
    
    //if ( _rejection_threshold > 0.0f ) {
     //   ss << "Rejection threshold set to " << boost::format("%-10d") % _rejection_threshold << std::endl;
    //}
    
    
    ss << "Population size (N) is " << _zparams.GetInt(fits_constants::PARAM_POPULATION_SIZE, -1);
    if ( _zparams.GetInt(fits_constants::PARAM_SAMPLE_SIZE, 0) > 0 ) {
        ss << " (sampled " << _zparams.GetInt(fits_constants::PARAM_SAMPLE_SIZE, 0) << ")";
    }
    ss << std::endl;
    
    if ( _single_mutrate_inferred ) {
        ss << "Inferred a single mutation rate." << std::endl;
    }
    
    ss << "====================" << std::endl;
    
    // first column - no header
    ss << boost::format("%-10s") % "";
    
    for (auto col=0; col<_num_alleles; ++col ) {
        ss << boost::format("%-10s") % ("allele" + std::to_string(col));
    }
    ss << boost::format("%-10s") % "minldist";
    ss << boost::format("%-10s") % "maxldist";
    ss << std::endl;
    
    for ( auto row=0; row<_num_alleles; ++row ) {
        
        ss << boost::format("%-10s") % ("allele" + std::to_string(row));
        
        
        for (auto col=0; col<_num_alleles; ++col ) {
            
            auto tmp_power = median_mutation_rates(row,col);
            
            if ( row == col ) {
                ss << boost::format("%-10s") % "---";
            }
            else {
                ss << boost::format("%-10.3d") % tmp_power;
            }
            
            
        }
        
        ss << boost::format("%-10.3d") % _distance_min;
        ss << boost::format("%-10.3d") % _distance_max;
        ss << std::endl;
    }
    ss << std::endl;
    
    if ( _zparams.GetInt( fits_constants::PARAM_DUMP_PARAMETERS, 0) > 1 ) {
        ss << _zparams.GetAllParameters() << std::endl;
    }
    
    return ss.str();

}
