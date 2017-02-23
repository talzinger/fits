//
//  ReportStats.cpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "ResultsStats.hpp"

ResultsStats::ResultsStats()
: _num_results(0),
_num_alleles(0),
allele_mean_fitness(_num_alleles),
allele_sd_fitness(_num_alleles),
allele_median_fitness(_num_alleles),
allele_min_fitness(_num_alleles, -1.0),
allele_max_fitness(_num_alleles, -1.0),
allele_min_95percentile_fitness(_num_alleles, -1.0),
allele_max_95percentile_fitness(_num_alleles, -1.0),
allele_pval(_num_alleles, 0),
lethal_counter(_num_alleles, 0),
deleterious_counter(_num_alleles, 0),
neutral_counter(_num_alleles, 0),
advantageous_counter(_num_alleles, 0),
lethal_percent(_num_alleles, -1),
deleterious_percent(_num_alleles, -1),
neutral_percent(_num_alleles, -1),
advantageous_percent(_num_alleles, -1),
allele_category(_num_alleles, AlleleCategory::Undefined),
_rejection_threshold(0.0f),
_results_count(0)
{
}

void ResultsStats::SetRejectionThreshold(FLOAT_TYPE new_val)
{
    _rejection_threshold = new_val;
}

FLOAT_TYPE ResultsStats::GetRejectionThreshold()
{
    return _rejection_threshold;
}

void ResultsStats::SetResultsCount( std::size_t count )
{
    _results_count = count;
}
/*
std::string ReportStatsFitness(const std::vector<SimulationResult>& result_vector, std::string filename = "")
{
    std::string output_str;
    
    if (result_vector.empty()) {
        std::cerr << "ReportStats: empty result vector." << std::endl;
        throw "ReportStats: empty result vector.";
    }
    
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    
    // open the file now, so if it'll fail we won't even start generating data
    if (!outfile.is_open() && !filename.empty() ) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    output_str += "Fitness Report\n";
    output_str += "===============\n";
    
    //cout << "Processing " << result_vector.size() << " results." << endl;
    output_str += "Processing " + std::to_string(result_vector.size()) + " results.\n";
    
    // int _num_alleles = result_vector[0].fitness_values.size();
    
    ResultsStats results_stats;
    results_stats.CalculateStatsFitness(result_vector);
    
    for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
        output_str += "\tallele" + std::to_string(i);
        //cout << "\t" << "allele" << i;
        
    }
    output_str += "\n";
    //cout << endl;
    
    //cout << "median";
    
    // replacing boost's format so that we could put this inside a string
    //std::setprecision(3);
    
    output_str += "median";
    for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
        std::ostringstream tmp_out;
        tmp_out << std::fixed << std::setprecision(REPORT_FLOAT_PRECISION) << results_stats.allele_median_fitness[i];
        output_str += "\t" + tmp_out.str();
        
        // cout << "\t" << boost::format("%1.2f") % results_stats.allele_median_fitness[i];
        // output_str += "\t" + to_string(results_stats.allele_median_fitness[i]);
    }
    output_str += "\n";
    // cout << endl;
    
    output_str += "mean";
    // cout << "mean";
    for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
        std::ostringstream tmp_out;
        tmp_out << std::fixed << std::setprecision(REPORT_FLOAT_PRECISION) << results_stats.allele_mean_fitness[i];
        output_str += "\t" + tmp_out.str();
        
        // cout << "\t" << boost::format("%1.2f") % results_stats.allele_mean_fitness[i];
        // output_str += "\t" + to_string(results_stats.allele_mean_fitness[i]);
    }
    output_str += "\n";
    //cout << endl;
    
    //cout << "lbound";
    output_str += "lbound";
    for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
        std::ostringstream tmp_out;
        tmp_out << std::fixed << std::setprecision(REPORT_FLOAT_PRECISION) << results_stats.allele_min_fitness[i];
        output_str += "\t" + tmp_out.str();
        
        
        // output_str += "\t" + to_string(results_stats.allele_min_fitness[i]);
        // cout << "\t" << boost::format("%1.2f") % results_stats.allele_min_fitness[i];
    }
    //cout << endl;
    output_str += "\n";
    
    output_str += "ubound";
    //cout << "ubound";
    for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
        std::ostringstream tmp_out;
        tmp_out << std::fixed << std::setprecision(REPORT_FLOAT_PRECISION) << results_stats.allele_max_fitness[i];
        output_str += "\t" + tmp_out.str();
        
        // output_str += "\t" + to_string(results_stats.allele_max_fitness[i]);
        // cout << "\t" << boost::format("%1.2f") % results_stats.allele_max_fitness[i];
    }
    //cout << endl;
    output_str += "\n";
    
    // distribution
    // todo: do this more generally with boost? only 4 bins so here it may be more work than benefit
 
    output_str += "\nDistribution:\n";
    // output_str += "\tlet.\tdel.\tneu.\tadv.\n";
    output_str += "\tdel.\tneu.\tadv.\n";
    
    for (int current_allele = 0; current_allele < results_stats._num_alleles; current_allele++) {
        
        output_str += "allele" + std::to_string(current_allele);
        //cout << "allele" << current_allele;
        
        // output_str += "\t" + to_string(results_stats.lethal_percent[current_allele]) + "%";
        output_str += "\t" + std::to_string(results_stats.deleterious_percent[current_allele]) + "%";
        output_str += "\t" + std::to_string(results_stats.neutral_percent[current_allele]) + "%";
        output_str += "\t" + std::to_string(results_stats.advantageous_percent[current_allele]) + "%";
        
        
        output_str += "\t=>";
        //cout << "\t=> ";
        
        switch ( results_stats.allele_category[current_allele] ) {
                
            case AlleleCategory::Undefined:
                output_str += "???";
                break;
                
            case AlleleCategory::Adventageous:
                output_str += "ADV";
                break;
                
            case AlleleCategory::Possible_advantageous:
                output_str += "?ADV";
                break;
                
            case AlleleCategory::Deleterious:
                output_str += "DEL";
                break;
                
            case AlleleCategory::Possible_deleterious:
                output_str += "?DEL";
                break;
                
            case AlleleCategory::Neutral:
                output_str += "NEU";
                break;
            case AlleleCategory::Possible_neutral:
                output_str += "?NEU";
                break;
                
            case AlleleCategory::Lethal:
                output_str += "LTH";
                break;
                
            case AlleleCategory::Possible_lethal:
                output_str += "?LTH";
                break;
                
            case AlleleCategory::WT:
                output_str += "WT";
                break;
        }
        
        output_str += "\n";
    }
    
    
    if ( outfile.is_open() ) {
        outfile << output_str;
    }
    
    std::cout << output_str;
    
    return output_str;
}


*/
