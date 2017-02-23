//
//  ReportStats.hpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#ifndef ReportStats_hpp
#define ReportStats_hpp

#include <vector>
#include <string>
#include <sstream>

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
#include "SimulationResult.hpp"

enum AlleleCategory {
    // cannot be assigned with fitness value
    WT,
    Undefined,
    
    // pass relaxed threshold
    Possible_lethal,
    Possible_deleterious,
    Possible_neutral,
    Possible_advantageous,
    
    // pass strict threshold
    Lethal,
    Deleterious,
    Neutral,
    Adventageous
};


class ResultsStats {
    
    const FLOAT_TYPE FITNESS_LETHAL_MIN = 0.0;
    const FLOAT_TYPE FITNESS_LETHAL_MAX = 0.5;
    const FLOAT_TYPE FITNESS_DELETERIOUS_MIN = 0.5;
    const FLOAT_TYPE FITNESS_DELETERIOUS_MAX = 0.95;
    const FLOAT_TYPE FITNESS_NEUTRAL_MIN = 0.95;
    const FLOAT_TYPE FITNESS_NEUTRAL_MAX = 1.05;
    const FLOAT_TYPE FITNESS_ADVANTAGEOUS_MIN = 1.05;
    const FLOAT_TYPE FITNESS_ADVANTAGEOUS_MAX = 2.0;
    const FLOAT_TYPE FITNESS_NEUTRAL = 1.0f;
    
    const int CATEGORY_INCLUSION_RELAXED_THRESHOLD = 50;
    const int CATEGORY_INCLUSION_STRICT_THRESHOLD = 95;

public:    
    ResultsStats();

    void SetResultsCount( std::size_t count );
    
    std::size_t _results_count;
    
    // These will calculate and store in internal variables
    // will later be private

    void CalculateStatsFitness(const std::vector<SimulationResult>& result_vector);
    void CalculateStatsPopulationSize(const std::vector<SimulationResult>& result_vector);
    void CalculateStatsMutation(const std::vector<SimulationResult>& result_vector);
    
     
    std::string GetSummaryFitness();
    std::string GetSummaryPopSize();
    std::string GetSummaryMutRate();
    
    void SetRejectionThreshold(FLOAT_TYPE new_val);
    FLOAT_TYPE GetRejectionThreshold();
    
    // Accumulators
    double _pop_min;
    double _pop_max;
    double _pop_mean;
    double _pop_sd;
    double _pop_median;
    

    FLOAT_TYPE _distance_min;
    FLOAT_TYPE _distance_max;
    FLOAT_TYPE _distance_mean;
    FLOAT_TYPE _distance_sd;
    FLOAT_TYPE _distance_median;
    
    
    std::string AlleleCategory2String( AlleleCategory category );
    
    /* general data */
    std::size_t _num_results;
    std::size_t _num_alleles;
    FLOAT_TYPE _rejection_threshold;

    /* fitness inference */
    std::vector<double> allele_mean_fitness;
    std::vector<double> allele_sd_fitness;
    std::vector<double> allele_median_fitness;
    
    std::vector<double> allele_min_fitness;
    std::vector<double> allele_max_fitness;
    
    std::vector<double> allele_min_95percentile_fitness;
    std::vector<double> allele_max_95percentile_fitness;
    
    // measure of disperssion
    std::vector<double> allele_MAD;
    
    
    // P(w<1|data)
    std::vector<double> allele_pval;
    
    std::vector<unsigned int> lethal_counter;
    std::vector<unsigned int> deleterious_counter;
    std::vector<unsigned int> neutral_counter;
    std::vector<unsigned int> advantageous_counter;
    
    std::vector<unsigned int> lethal_percent;
    std::vector<unsigned int> deleterious_percent;
    std::vector<unsigned int> neutral_percent;
    std::vector<unsigned int> advantageous_percent;
    
    std::vector<AlleleCategory> allele_category;
    
    /* mutation rate inference */
    boost::numeric::ublas::matrix<float> min_mutation_rates;
    boost::numeric::ublas::matrix<float> max_mutation_rates;
    boost::numeric::ublas::matrix<float> mean_mutation_rates;
    boost::numeric::ublas::matrix<float> median_mutation_rates;
    boost::numeric::ublas::matrix<float> normalized_median_mutation_rates;
    
    /* population size inference */
    int max_population_size;
    int min_population_size;
    int mean_population_size;
    int median_population_size;
    int normalized_median_population_size;
};


#endif /* ReportStats_hpp */
