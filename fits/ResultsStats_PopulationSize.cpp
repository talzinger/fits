//
//  ResultsStats_PopulationSize.cpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "ResultsStats.hpp"


void ResultsStats::CalculateStatsPopulationSize(const std::vector<SimulationResult>& result_vector)
{
    boost::accumulators::accumulator_set<
    float,
    boost::accumulators::stats<
    boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_distance;
    
    boost::accumulators::accumulator_set<
    float,
    boost::accumulators::stats<
    boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_population_size;
    
    for ( auto sim_result : result_vector ) {
        acc_distance(sim_result.distance_from_actual);
        
        acc_population_size(sim_result.N);
    }
    
    _pop_min = boost::accumulators::min(acc_population_size);
    _pop_max = boost::accumulators::max(acc_population_size);
    _pop_mean = boost::accumulators::mean(acc_population_size);
    _pop_sd = std::sqrt(boost::accumulators::variance(acc_population_size));
    _pop_median = boost::accumulators::median(acc_population_size);
    
    _distance_min = boost::accumulators::min(acc_distance);
    _distance_max = boost::accumulators::max(acc_distance);
    _distance_mean = boost::accumulators::mean(acc_distance);
    _distance_sd = std::sqrt(boost::accumulators::variance(acc_distance));
    _distance_median = boost::accumulators::median(acc_distance);
}


std::string ResultsStats::GetSummaryPopSize()
{
    std::stringstream ss;
    
    ss << "Population Size Report" << std::endl;
    ss << "=======================" << std::endl;
    
    ss << "Rejection threshold set to " << boost::format("%-10d") % _rejection_threshold << std::endl;
    
    ss << boost::format("%-12s") % "median";
    ss << boost::format("%-12s") % "mean";
    ss << boost::format("%-12s") % "low";
    ss << boost::format("%-12s") % "high";
    ss << boost::format("%-10s") % "min dist";
    ss << boost::format("%-10s") % "max dist";
    ss << std::endl;
    
    ss << boost::format("%-12.2e") % _pop_median;
    ss << boost::format("%-12.2e") % _pop_mean;
    ss << boost::format("%-12.2e") % _pop_min;
    ss << boost::format("%-12.2e") % _pop_max;
    ss << boost::format("%-10s") % _distance_min;
    ss << boost::format("%-10s") % _distance_max;
    ss << std::endl;
    ss << std::endl;
    
    return ss.str();
}
