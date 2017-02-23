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
    
    boost::accumulators::accumulator_set<
    float,
    boost::accumulators::stats<
    boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_distance;
    
    boost::numeric::ublas::matrix< boost::accumulators::accumulator_set<
    float,
    boost::accumulators::stats<
    boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > acc_matrix(_num_alleles,_num_alleles);
    
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
            }
        }
    }
    
    for ( auto row=0; row<_num_alleles; ++row ) {
        
        for ( auto col=0; col<_num_alleles; ++col ) {
            
            min_mutation_rates(row,col) = boost::accumulators::min( acc_matrix(row,col) );
            max_mutation_rates(row,col) = boost::accumulators::max( acc_matrix(row,col) );
            
            mean_mutation_rates(row,col) = boost::accumulators::mean( acc_matrix(row,col) );
            median_mutation_rates(row,col) = boost::accumulators::median( acc_matrix(row,col) );
        }
    }
    
    _distance_min = boost::accumulators::min(acc_distance);
    _distance_max = boost::accumulators::max(acc_distance);
    _distance_mean = boost::accumulators::mean(acc_distance);
    _distance_sd = std::sqrt(boost::accumulators::variance(acc_distance));
    _distance_median = boost::accumulators::median(acc_distance);
    
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
    
    ss << "Mutation Rate Report" << std::endl;
    ss << "=====================" << std::endl;
    
    ss << "Rejection threshold set to " << boost::format("%-10f") % _rejection_threshold << std::endl;
    
    // first column - no header
    ss << boost::format("%-10s") % "";
    
    for (auto col=0; col<_num_alleles; ++col ) {
        ss << boost::format("%-15s") % ("allele" + std::to_string(col));
    }
    ss << boost::format("%-15s") % "min dist";
    ss << boost::format("%-15s") % "max dist";
    ss << std::endl;
    
    for ( auto row=0; row<_num_alleles; ++row ) {
        
        ss << boost::format("%-10s") % ("allele" + std::to_string(row));
        
        for (auto col=0; col<_num_alleles; ++col ) {
            
            //ss << "1.00e" << ( median_mutation_rates(row,col) > 0 ? "+" : "-" ) << std::round(median_mutation_rates(row,col));
            auto tmp_power = median_mutation_rates(row,col);
            auto complement = 1.0f - tmp_power;
            
            // I want to prevent from the mutation to appear as rate of 1.00 dur to rounding
            //if ( col==row ) {
            //    ss << "1-" << boost::format("%-20.10s") % complement;
            //}
            //else {
                ss << boost::format("%-15.10f") % tmp_power;
            //}
            
        }
        
        ss << boost::format("%-15.10f") % _distance_min;
        ss << boost::format("%-15.10f") % _distance_max;
        ss << std::endl;
    }
    ss << std::endl;
    
    /* No need for that, I've already normalized right after prior sampled
     
    ss << "normalized inferred rates:" << std::endl;
    ss << boost::format("%-10s") % "";
    for (auto col=0; col<_num_alleles; ++col ) {
        ss << boost::format("%-10s") % ("allele" + std::to_string(col));
    }
    ss << std::endl;
    
    for ( auto row=0; row<_num_alleles; ++row ) {
        
        ss << boost::format("%-10s") % ("allele" + std::to_string(row));
        
        for (auto col=0; col<_num_alleles; ++col ) {
            
            //ss << "1.00e" << ( median_mutation_rates(row,col) > 0 ? "+" : "-" ) << median_mutation_rates(row,col);
            //ss << boost::format("%-10.2e") % normalized_median_mutation_rates(row,col);
            auto tmp_power =  normalized_median_mutation_rates(row,col);
            ss << boost::format("%-10.2e") % tmp_power;
        }
        
        ss << std::endl;
    }
    
    ss << std::endl;
    */
    return ss.str();

}
