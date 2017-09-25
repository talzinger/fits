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
    _num_results = result_vector.size();
    
    boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_distance;
    
    boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_population_size;
    

    std::vector<int> popsize_storage;
    for ( auto sim_result : result_vector ) {
        acc_distance(sim_result.distance_from_actual);
        
        acc_population_size(sim_result.N);
        
        std::cout << "N=" << sim_result.N << std::endl;
        popsize_storage.push_back(sim_result.N);
    }
    
    _pop_min = boost::accumulators::min(acc_population_size);
    _pop_max = boost::accumulators::max(acc_population_size);
    _pop_mean = boost::accumulators::mean(acc_population_size);
    _pop_sd = std::sqrt(boost::accumulators::variance(acc_population_size));
    
    _pop_median = GetMedian(popsize_storage);
    
    _distance_min = boost::accumulators::min(acc_distance);
    _distance_max = boost::accumulators::max(acc_distance);
    _distance_mean = boost::accumulators::mean(acc_distance);
    _distance_sd = std::sqrt(boost::accumulators::variance(acc_distance));
}


std::string ResultsStats::GetSummaryPopSize()
{
    std::stringstream ss;
    
    ss.imbue(std::locale(fits_constants::used_locale));
    
    
    ss << "Population Size Report" << std::endl;
    ss << GetSummaryHeader();
    
    //ss << "FITS v"<< fits_constants::current_version_str << std::endl;
    
    //auto current_time_raw = std::chrono::system_clock::now();
    //auto current_time = std::chrono::system_clock::to_time_t(current_time_raw);
    //auto current_time_final = *std::localtime(&current_time);
    //ss << std::put_time(&current_time_final, "%F %T") << std::endl;
    
    //std::cout << "Simulation results used for calculations: " << _num_results << std::endl;
    
    //ss << "=======================" << std::endl;
    
    //if ( _rejection_threshold > 0.0f ) {
      //  ss << "Rejection threshold set to " << boost::format("%-10d") % _rejection_threshold << std::endl;
    //}
    
    
    if ( _single_mutrate_used ) {
        ss << "Used a single mutation rate." << std::endl;
    }
    
    ss << "====================" << std::endl;
    
    ss << boost::format("%-12s") % "median";
    ss << boost::format("%-12s") % "mean";
    ss << boost::format("%-12s") % "low";
    ss << boost::format("%-12s") % "high";
    ss << boost::format("%-10s") % "minldist";
    ss << boost::format("%-10s") % "maxldist";
    ss << std::endl;
    
    ss << boost::format("%-12.2e") % _pop_median;
    ss << boost::format("%-12.2e") % _pop_mean;
    ss << boost::format("%-12.2e") % _pop_min;
    ss << boost::format("%-12.2e") % _pop_max;
    ss << boost::format("%-10s") % _distance_min;
    ss << boost::format("%-10s") % _distance_max;
    ss << std::endl;
    ss << std::endl;
    
    if ( _zparams.GetInt( fits_constants::PARAM_DUMP_PARAMETERS, 0) > 1 ) {
        ss << _zparams.GetAllParameters() << std::endl;
    }
    
    return ss.str();
}
