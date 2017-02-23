//
//  ResultsStats_Fitness.cpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "ResultsStats.hpp"


void ResultsStats::CalculateStatsFitness(const std::vector<SimulationResult>& result_vector)
{
    auto wt_allele = result_vector[0].wt_index;
    _num_alleles = result_vector[0].fitness_values.size();
    
    // This is meant to be used to measure deviations of the median from the actual result
    // Median of abs(Xi-median(X)) is called Median Absolute Deviation (MAD) and is
    // a measure of distribution and thus relaiability of the result
    std::vector< boost::accumulators::accumulator_set<
    float,
    boost::accumulators::stats<
    boost::accumulators::tag::median,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > acc_vec_MAD(_num_alleles);

    
    boost::accumulators::accumulator_set<
    float,
    boost::accumulators::stats<
    boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_distance;
    
    std::vector< boost::accumulators::accumulator_set<
    float,
    boost::accumulators::stats<
    boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > acc_vec_fitness(_num_alleles);

    allele_MAD.resize(_num_alleles, 0.0f);
    allele_max_fitness.resize(_num_alleles, 0.0f);
    allele_mean_fitness.resize(_num_alleles, 0.0f);
    allele_sd_fitness.resize(_num_alleles, 0.0f);
    allele_median_fitness.resize(_num_alleles, 0.0f);
    allele_min_fitness.resize(_num_alleles, 0.0f);
    
    allele_min_95percentile_fitness.resize(_num_alleles, 0.0f);
    allele_max_95percentile_fitness.resize(_num_alleles, 0.0f);
    
    allele_pval.resize(_num_alleles, 0.0f);
    lethal_counter.resize(_num_alleles, 0.0f);
    deleterious_counter.resize(_num_alleles, 0.0f);
    neutral_counter.resize(_num_alleles, 0.0f);
    advantageous_counter.resize(_num_alleles, 0.0f);
    
    lethal_percent.resize(_num_alleles, 0.0f);
    deleterious_percent.resize(_num_alleles, 0.0f);
    neutral_percent.resize(_num_alleles, 0.0f);
    advantageous_percent.resize(_num_alleles, 0.0f);
    
    allele_category.resize(_num_alleles, AlleleCategory::Undefined);
    
    // sort according to distance - used to bound 95 percentile
    // std::sort(result_vector.begin(), result_vector.end());
    
    // get the index of the upper bound of the 95th percentile
    //auto tmp_ubound_idx_95th_percentile = std::ceil(static_cast<double>(result_vector.size()) * 0.95);
    //auto ubound_idx_95th_percentile = static_cast<int>(tmp_ubound_idx_95th_percentile);
    
    
    for (auto current_sim_data : result_vector) {
        
        acc_distance(current_sim_data.distance_from_actual);
    
        for (auto current_allele = 0; current_allele<current_sim_data.fitness_values.size(); current_allele++) {
         
            auto tmpval = current_sim_data.fitness_values[current_allele];
           
            acc_vec_fitness[current_allele](tmpval);
            
            // the ranges are written explicitlly for readibility, by be redundant
            
           
            // 2016-10-04 - removed categories for lethal and neutral
            if ( tmpval > FITNESS_NEUTRAL ) {
                ++advantageous_counter[current_allele];
            }
            else if ( tmpval < FITNESS_NEUTRAL ) {
                ++deleterious_counter[current_allele];
            }
            else {
                ++neutral_counter[current_allele];
            }
        } // for with _counters
        
    }
    
    _distance_min = boost::accumulators::min(acc_distance);
    _distance_max = boost::accumulators::max(acc_distance);
    _distance_mean = boost::accumulators::mean(acc_distance);
    _distance_sd = std::sqrt(boost::accumulators::variance(acc_distance));
    _distance_median = boost::accumulators::median(acc_distance);
    
    for (auto current_allele = 0; current_allele<_num_alleles; ++current_allele) {
        
        // added boost::accumulators to prevent any ambiguity with std library
        allele_median_fitness[current_allele] = boost::accumulators::median(acc_vec_fitness[current_allele]);
        allele_mean_fitness[current_allele] = boost::accumulators::mean(acc_vec_fitness[current_allele]);
        allele_sd_fitness[current_allele] = sqrt(boost::accumulators::variance(acc_vec_fitness[current_allele]));
        allele_max_fitness[current_allele] = boost::accumulators::max(acc_vec_fitness[current_allele]);
        allele_min_fitness[current_allele] = boost::accumulators::min(acc_vec_fitness[current_allele]);
        
        
        // calculate Bayesian pval - P(w<1|data)
        allele_pval[current_allele] =
        static_cast<double>(lethal_counter[current_allele] +
                            deleterious_counter[current_allele]) * 100.0 /
        static_cast<double>(_num_results);
        
        // 95th percentile
        //allele_min_95percentile_fitness[current_allele] = result_vector[0].fitness_values[current_allele];
        //allele_max_95percentile_fitness[current_allele] = result_vector[ubound_idx_95th_percentile].fitness_values[current_allele];
    }
    
    
    // Calculate MAD for the alleles
    for (auto current_sim_data : result_vector) {
        for (auto current_allele = 0; current_allele<current_sim_data.fitness_values.size(); current_allele++) {
            
            auto current_allele_value = current_sim_data.fitness_values[current_allele];
            auto current_allele_median = allele_median_fitness[current_allele];
            
            auto tmpval = std::fabs(current_allele_value-current_allele_median);
            
            acc_vec_MAD[current_allele](tmpval);
        }
    }
    
    for (auto current_allele = 0; current_allele < _num_alleles; ++current_allele) {
        allele_MAD[current_allele] = boost::accumulators::median(acc_vec_MAD[current_allele]);
    }
    
    for (auto current_allele = 0; current_allele < _num_alleles; ++current_allele) {
        
        // first so it will be the default
        allele_category[current_allele] = AlleleCategory::Undefined;
        /*
         lethal_percent[current_allele] = lethal_counter[current_allele] * 100 / result_vector.size();
         if (lethal_percent[current_allele] >= CATEGORY_INCLUSION_STRICT_THRESHOLD)
         allele_category[current_allele] = AlleleCategory::Lethal;
         else if (lethal_percent[current_allele] >= CATEGORY_INCLUSION_RELAXED_THRESHOLD)
         allele_category[current_allele] = AlleleCategory::Possible_lethal;
         */
        deleterious_percent[current_allele] = deleterious_counter[current_allele] * 100 / result_vector.size();
        if (deleterious_percent[current_allele] >= CATEGORY_INCLUSION_STRICT_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Deleterious;
        else if (deleterious_percent[current_allele] >= CATEGORY_INCLUSION_RELAXED_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Possible_deleterious;
        
        neutral_percent[current_allele] = neutral_counter[current_allele] * 100 / result_vector.size();
        if (neutral_percent[current_allele] >= CATEGORY_INCLUSION_STRICT_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Neutral;
        else if (neutral_percent[current_allele] >= CATEGORY_INCLUSION_RELAXED_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Possible_neutral;
        
        advantageous_percent[current_allele] = advantageous_counter[current_allele] * 100 / result_vector.size();
        if (advantageous_percent[current_allele] >= CATEGORY_INCLUSION_STRICT_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Adventageous;
        else if (advantageous_percent[current_allele] >= CATEGORY_INCLUSION_RELAXED_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Possible_advantageous;
        
        // last in order to override any fitness value assignments
        if (current_allele == wt_allele) {
            allele_category[current_allele] = AlleleCategory::WT;
        }
    }
}


std::string ResultsStats::GetSummaryFitness()
{
    
    std::stringstream ss;
    
    ss << "Fitness Report" << std::endl;
    ss << "===============" << std::endl;
    
    ss << "Rejection threshold set to " << boost::format("%-10d") % _rejection_threshold << std::endl;
    
    ss << boost::format("%-10s") % "allele";
    ss << boost::format("%-10s") % "median";
    ss << boost::format("%-10s") % "mean";
    ss << boost::format("%-10s") % "low";
    ss << boost::format("%-10s") % "high";
    ss << boost::format("%-10s") % "DEL(%)";
    ss << boost::format("%-10s") % "NEU(%)";
    ss << boost::format("%-10s") % "ADV(%)";
    ss << boost::format("%-10s") % "category";
    ss << boost::format("%-10s") % "min dist";
    ss << boost::format("%-10s") % "max dist";
    ss << boost::format("%-10s") % "MAD";
    ss << std::endl;
    
    for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
        ss << boost::format("%-10d") % current_allele;
        ss << boost::format("%-10.3d") % allele_median_fitness[current_allele];
        ss << boost::format("%-10.3d") % allele_mean_fitness[current_allele];
        ss << boost::format("%-10.3d") % allele_min_fitness[current_allele];
        ss << boost::format("%-10.3d") % allele_max_fitness[current_allele];
        ss << boost::format("%-10s") % deleterious_percent[current_allele];
        ss << boost::format("%-10s") % neutral_percent[current_allele];
        ss << boost::format("%-10s") % advantageous_percent[current_allele];
        ss << boost::format("%-10s") % AlleleCategory2String(allele_category[current_allele]);
        ss << boost::format("%-10s") % _distance_min;
        ss << boost::format("%-10s") % _distance_max;
        ss << boost::format("%-10s") % allele_MAD[current_allele];
        ss << std::endl;
    }
    
    ss << std::endl;
    
    return ss.str();
}


std::string ResultsStats::AlleleCategory2String( AlleleCategory category )
{
    std::string tmp_str = "";
    
    switch ( category ) {
            
        case AlleleCategory::Undefined:
            tmp_str = "???";
            break;
            
        case AlleleCategory::Adventageous:
            tmp_str = "ADV";
            break;
            
        case AlleleCategory::Possible_advantageous:
            tmp_str = "?ADV";
            break;
            
        case AlleleCategory::Deleterious:
            tmp_str = "DEL";
            break;
            
        case AlleleCategory::Possible_deleterious:
            tmp_str = "?DEL";
            break;
            
        case AlleleCategory::Neutral:
            tmp_str = "NEU";
            break;
        case AlleleCategory::Possible_neutral:
            tmp_str = "?NEU";
            break;
            
        case AlleleCategory::Lethal:
            tmp_str = "LTH";
            break;
            
        case AlleleCategory::Possible_lethal:
            tmp_str = "?LTH";
            break;
            
        case AlleleCategory::WT:
            tmp_str = "WT";
            break;
    }
    
    return tmp_str;
}

