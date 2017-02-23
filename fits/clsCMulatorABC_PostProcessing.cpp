/****************************
  Post Processing of ABC
  Assert results of the ABC
 ****************************/

#include "clsCMulatorABC.h"

std::vector<int> clsCMulatorABC::GetUniqueIndexSet( int num_items )
{
    std::cout << "Generating " << num_items << " unique items..." << std::endl;
    
    std::vector<bool> chosen_array(num_items, false);
    
    auto chosen_counter = 0;;
    boost::random::uniform_int_distribution<> local_int_distrib(0, num_items);
    
    while (chosen_counter<num_items) {
        
        auto chosen_index = local_int_distrib(_boost_gen);
        
        if (!chosen_array[chosen_index]) {
            
            chosen_array[chosen_index] = true;
            ++chosen_counter;
            
            //std::cout << chosen_index << std::endl;
        }
        
        
    }
    
    std::vector<int> return_vec(num_items, 0);
    
    for ( auto i=0; i<chosen_array.size(); ++i ) {
        
        if ( chosen_array[i] ) {
            return_vec.push_back(i);
        }
    }
    
    return return_vec;
    
}

std::vector<FLOAT_TYPE> clsCMulatorABC::GetSDPerAllele( std::size_t start_idx, std::size_t end_idx )
{
    
    std::vector<FLOAT_TYPE> allele_sd_vec( _num_alleles, 0.0f);
    
    std::vector< boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > allele_accumulator_vec;
    
    allele_accumulator_vec.resize(_num_alleles);
    
    // for each of the unchosen-ones [ ones that were not pseudo data
    for ( auto current_idx=start_idx; current_idx<end_idx; ++current_idx ) {
        
        for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
            
            boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col( _simulation_result_vector[current_idx].sim_data_matrix, current_allele );
            
            for ( auto current_freq : current_col ) {
                allele_accumulator_vec[current_allele]( current_freq );
            }
        }
    }
    
    for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
        allele_sd_vec[current_allele] = std::sqrt( boost::accumulators::variance(allele_accumulator_vec[current_allele]) );
    }
    
    return allele_sd_vec;
}

// Calculate Median Absolute Deviation (MAD)
std::vector<FLOAT_TYPE> clsCMulatorABC::GetMADPerAllele( std::size_t start_idx, std::size_t end_idx )
{
    
    std::vector<FLOAT_TYPE> allele_median_vec( _num_alleles, 0.0f);
    std::vector<FLOAT_TYPE> allele_mad_vec( _num_alleles, 0.0f);
    
    std::vector< boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > allele_accumulator_vec;
    
    std::vector< boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > distance_accumulator_vec;
    
    
    // First, calculate the median
    allele_accumulator_vec.resize(_num_alleles);
    for ( auto current_idx=start_idx; current_idx<end_idx; ++current_idx ) {
        
        for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
            
            boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col( _simulation_result_vector[current_idx].sim_data_matrix, current_allele );
            
            for ( auto current_freq : current_col ) {
                allele_accumulator_vec[current_allele]( current_freq );
            }
        }
    }
    
    // Store the allele medians
    for (auto current_allele = 0; current_allele < _num_alleles; ++current_allele) {
        allele_median_vec[current_allele] = boost::accumulators::median( allele_accumulator_vec[current_allele] );
    }
    
    // Now, accumulate the distances of values from the median
    distance_accumulator_vec.resize(_num_alleles);
    for ( auto current_idx=start_idx; current_idx<end_idx; ++current_idx ) {
        
        for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
            
            boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col( _simulation_result_vector[current_idx].sim_data_matrix, current_allele );
            
            for ( auto current_freq : current_col ) {
                
                auto current_allele_value = current_freq;
                auto current_allele_median = allele_median_vec[current_allele];
                
                auto tmpval = std::fabs( current_allele_value - current_allele_median );
                
                distance_accumulator_vec[current_allele](tmpval);
            }
        }
    }
    
    // Finally, return MAD for each allele
    for (auto current_allele = 0; current_allele < _num_alleles; ++current_allele) {
        allele_mad_vec[current_allele] = boost::accumulators::median( distance_accumulator_vec[current_allele] );
    }
    
    return allele_mad_vec;
}

void clsCMulatorABC::DivideEachAllele( std::size_t start_idx, std::size_t end_idx, std::vector<FLOAT_TYPE> value_vector )
{
    
    for ( auto current_idx=start_idx; current_idx<end_idx; ++current_idx ) {
        
        //std::cout << "before: " << _simulation_result_vector[current_idx].sim_data_matrix << std::endl;
        
        for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
            
            boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col( _simulation_result_vector[current_idx].sim_data_matrix, current_allele );
            
            current_col = current_col / value_vector[current_allele];
        }
        //std::cout << "after: " <<_simulation_result_vector[current_idx].sim_data_matrix << std::endl;
    }
    
}

void clsCMulatorABC::ScaleSimulationResults()
{
    auto sd_vector = GetSDPerAllele(0, _simulation_result_vector.size());
    
    std::cout << "Whole vector sd values: ";
    for ( auto tmpval : sd_vector ) { std::cout << "\t" << tmpval; }
    std::cout << std::endl;
    
    DivideEachAllele(0, _simulation_result_vector.size(), sd_vector);
}

void clsCMulatorABC::DoCoverageTest()
{
    std::cout << "Coverage" << std::endl;
    
    std::size_t acceptance_count = _zparams.GetFloat("_acceptance_rate", 0.01f) * _repeats;
    
    // instead of going through all
    std::cout << "acceptance count is " << acceptance_count << std::endl;
    // std::size_t num_sims = _simulation_result_vector.size() - acceptance_count;
    
    std::vector<std::size_t> fitness_underestimate_count_vec;
    std::vector<std::size_t> mutrate_underestimate_count_vec;
    std::size_t popsize_underestimate_count = 0;
    
    boost::numeric::ublas::matrix<FLOAT_TYPE> fitness_pvals(1,1);
    
    // Little sorting - make sure the first 100 results are the best of the best
    // these would be considered pseudo-data
    std::nth_element(_simulation_result_vector.begin(),
                     _simulation_result_vector.begin() + acceptance_count,
                     _simulation_result_vector.end());
    
    // for all the results not considerer pseudo-data
    // - compute the SD of frequencies for each allele
    // - scale the values accordingly
    // auto sd_vector = GetSDPerAllele(100, _simulation_result_vector.size());
    // DivideEachAllele(num_pseudo_datasets, _simulation_result_vector.size(), sd_vector);
    
    
    // scale the pseudo-data
    //for ( auto current_pseudo_dataset=0; current_pseudo_dataset<num_pseudo_datasets; ++current_pseudo_dataset ) {
    
    // auto pseudo_sd_vector = _simulation_result_vector[current_pseudo_dataset].GetSDForEachAllele();
    // _simulation_result_vector[current_pseudo_dataset].DivideEachAllele(pseudo_sd_vector);
    
    //    _simulation_result_vector[current_pseudo_dataset].DivideEachAllele(sd_vector);
    //}
    
    // write the pvals t file
    std::string filename = "fitness_pvals.txt";
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    std::string tmp_str = "";
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    outfile << "dataset";
    for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
        outfile << "\tallele" << current_allele;
    }
    outfile << std::endl;
    
    
    // now process each of the pseudo-datasets
    
    // TODO: create a vector of pairs{ index_of_sim_result, distance_from_current_pseudo
    // this paves the way for parallelization
    // the function in async should return the number of simulations for calculation of the pvalue
    //for ( auto current_pseudo_dataset=0; current_pseudo_dataset<acceptance_count; ++current_pseudo_dataset ) {
    
    // 2017-02-16 Because the coverage takes so much time, I'll take a smaller sample from
    // the posterior to check this. Haven't searched for justification yet.
    // NOTE: I won't take the best 10 but rather random set of datasets
    //for ( auto current_pseudo_dataset=0; current_pseudo_dataset<10; ++current_pseudo_dataset ) {
    auto to_be_sampled = ( acceptance_count > 100 ? 100 : acceptance_count );
    auto RandomPseudoDatasetsList = GetUniqueIndexSet( to_be_sampled );
    for ( auto current_pseudo_dataset : RandomPseudoDatasetsList ) {
        
        auto start_time_dataset = std::chrono::system_clock::now();
        
        // calculate the distance from the current dataset to the the simulation results not
        // taken as pseudo data
        for ( auto current_scaled_sim=acceptance_count;
              current_scaled_sim<_simulation_result_vector.size();
              ++current_scaled_sim ) {
            
            auto tmpdist = CMulatorToolbox::GetDistanceSimActual( _simulation_result_vector[current_pseudo_dataset].sim_data_matrix,
                                                                 _simulation_result_vector[current_scaled_sim].sim_data_matrix );
            
            _simulation_result_vector[current_scaled_sim].distance_from_actual = tmpdist;
        }
        
        // bound the simulations best fitting the current pseudo-dataset
        std::nth_element( _simulation_result_vector.begin() + acceptance_count,
                         _simulation_result_vector.begin() + acceptance_count + acceptance_count - 1,
                         _simulation_result_vector.end() );
        
        // go through the best fitting simulations and calculate the pvalue
        for ( auto current_scaled_sim=acceptance_count;
             current_scaled_sim<acceptance_count+acceptance_count;
             ++current_scaled_sim ) {
            
            switch (_factor_to_infer) {
                case Fitness:
                    
                    if ( fitness_underestimate_count_vec.empty() ) {
                        
                        fitness_underestimate_count_vec.resize( _num_alleles, 0 );
                        fitness_pvals.resize( acceptance_count, _num_alleles );
                    }
                    
                    for ( auto current_allele=0; current_allele<fitness_underestimate_count_vec.size(); ++current_allele ) {
                        if ( _simulation_result_vector[current_scaled_sim].fitness_values[current_allele] <
                            _simulation_result_vector[current_pseudo_dataset].fitness_values[current_allele] ) {
                            ++fitness_underestimate_count_vec[current_allele];
                        }
                    }
                    break;
                    
                case PopulationSize:
                    
                    break;
                    
                case MutationRate:
                    
                    break;
            }
        }
        
        // this is the pval for the current dataset
        for ( auto current_allele=0; current_allele<fitness_underestimate_count_vec.size(); ++current_allele ) {
            
            // TODO: why the corrections? (van der vaart did the +1 +2)
            auto count_float = static_cast<FLOAT_TYPE>(fitness_underestimate_count_vec[current_allele]) + 1.0f;
            auto total_float = static_cast<FLOAT_TYPE>(acceptance_count) + 2.0f;
            auto current_pval = count_float / total_float;
            
            fitness_pvals(current_pseudo_dataset, current_allele) = current_pval;
            
            fitness_underestimate_count_vec[current_allele] = 0;
        }
        
        // print the pvals for the current dataset
        outfile << current_pseudo_dataset;
        for ( auto current_allele=0; current_allele<fitness_underestimate_count_vec.size(); ++current_allele ) {
            outfile << "\t" << fitness_pvals(current_pseudo_dataset, current_allele);
        }
        outfile << std::endl;
        
        int remaining_sets = static_cast<int>(acceptance_count - current_pseudo_dataset - 1);
        auto end_time_dataset = std::chrono::system_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_dataset - start_time_dataset);
        float calculated_speed = 1.0f / static_cast<float>(elapsed_ms.count());
        calculated_speed *= 1000.0f; // 1/milisecond to 1/second
        calculated_speed *= 60.0f; // 1/second to 1/min
        
        std::cout << "remaining sets " << remaining_sets << std::endl;
        std::cout << "rcalculated_speed " << calculated_speed << std::endl;
        
        int minutes_remaining = remaining_sets / calculated_speed;
        
        auto duration_remaining = std::chrono::minutes(minutes_remaining);
        auto current_time = std::chrono::system_clock::now();
        
        auto completion_ETA = current_time + duration_remaining;
        
        auto completion_ETA_timet = std::chrono::system_clock::to_time_t(completion_ETA);
        auto completion_ETA_tm = *std::localtime(&completion_ETA_timet);
        
        std::cout << "Remaining " << remaining_sets
        << " sets (rate of "
        << std::round(calculated_speed)
        << "sets/min) - ETA is "
        << std::put_time(&completion_ETA_tm, "%c")
        << std::endl;
        
    }// pseudo datasets
    
    outfile.close();
    
}

void clsCMulatorABC::DoCoverageCook()
{
    std::cout << "Coverage Cook" << std::endl;
    
    std::size_t acceptance_count = _zparams.GetFloat("_acceptance_rate", 0.01f) * _repeats;
    
    std::cout << "acceptance count is " << acceptance_count << std::endl;
    
    std::vector<std::size_t> fitness_underestimate_count_vec;
    //std::vector<std::size_t> mutrate_underestimate_count_vec;
    //std::size_t popsize_underestimate_count = 0;
    
    boost::numeric::ublas::matrix<FLOAT_TYPE> fitness_pvals(1,1);
    
    // Little sorting - make sure the first 100 results are the best of the best
    // these would be considered pseudo-data
    std::nth_element(_simulation_result_vector.begin(),
                     _simulation_result_vector.begin() + acceptance_count,
                     _simulation_result_vector.end());
    
    // write the pvals output file
    std::string filename = "fitness_pvals.txt";
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    std::string tmp_str = "";
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    outfile << "dataset";
    for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
        outfile << "\tallele" << current_allele;
    }
    outfile << std::endl;
    
    // init
    switch (_factor_to_infer) {
        case Fitness:
            
            if ( fitness_underestimate_count_vec.empty() ) {
                
                fitness_underestimate_count_vec.resize( _num_alleles, 0 );
                fitness_pvals.resize( acceptance_count, _num_alleles );
            }
            
            break;
            
        case PopulationSize:
            
            break;
            
        case MutationRate:
            
            break;
    }
    
    
    // now process each of the pseudo-datasets
    /*std::vector< std::future< std::vector<std::size_t> > > fut_vec(100);
    for ( auto current_pseudo_dataset=0; current_pseudo_dataset<100; ++current_pseudo_dataset ) {
        fut_vec[current_pseudo_dataset] = std::async( std::launch::async,
                                                     &clsCMulatorABC::CoverageSingleDataset,
                                                     this,
                                                     current_pseudo_dataset,
                                                     100,
                                                     _simulation_result_vector.size() );
    }*/
    
    for ( auto current_pseudo_dataset=0; current_pseudo_dataset<100; ++current_pseudo_dataset ) {
        
        auto start_time_dataset = std::chrono::system_clock::now();
        
        std::cout << "Waiting for one..." << std::endl;
        //fitness_underestimate_count_vec = fut_vec[current_pseudo_dataset].get();
        std::cout << "Got one." << std::endl;
        
        fitness_underestimate_count_vec = CoverageSingleDataset( current_pseudo_dataset, 0, 100 );
        
        // this is the pval for the current dataset
        for ( auto current_allele=0; current_allele<fitness_underestimate_count_vec.size(); ++current_allele ) {
            
            // TODO: why the corrections? (van der vaart did the +1 +2)
            auto count_float = static_cast<FLOAT_TYPE>(fitness_underestimate_count_vec[current_allele]) + 1.0f;
            auto total_float = static_cast<FLOAT_TYPE>(acceptance_count) + 2.0f;
            auto current_pval = count_float / total_float;
            
            fitness_pvals(current_pseudo_dataset, current_allele) = current_pval;
            
            fitness_underestimate_count_vec[current_allele] = 0;
        }
        
        // print the pvals for the current dataset
        outfile << current_pseudo_dataset;
        for ( auto current_allele=0; current_allele<fitness_underestimate_count_vec.size(); ++current_allele ) {
            outfile << "\t" << fitness_pvals(current_pseudo_dataset, current_allele);
        }
        outfile << std::endl;
        
        int remaining_sets = static_cast<int>(acceptance_count - current_pseudo_dataset - 1);
        auto end_time_dataset = std::chrono::system_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_dataset - start_time_dataset);
        float calculated_speed = 1.0f / static_cast<float>(elapsed_ms.count());
        calculated_speed *= 1000.0f; // 1/milisecond to 1/second
        calculated_speed *= 60.0f; // 1/second to 1/min
        
        std::cout << "remaining sets " << remaining_sets << std::endl;
        std::cout << "rcalculated_speed " << calculated_speed << std::endl;
        
        int minutes_remaining = remaining_sets / calculated_speed;
        
        auto duration_remaining = std::chrono::minutes(minutes_remaining);
        auto current_time = std::chrono::system_clock::now();
        
        auto completion_ETA = current_time + duration_remaining;
        
        auto completion_ETA_timet = std::chrono::system_clock::to_time_t(completion_ETA);
        auto completion_ETA_tm = *std::localtime(&completion_ETA_timet);
        
        std::cout << "Remaining " << remaining_sets
        << " sets (rate of "
        << std::round(calculated_speed)
        << "sets/min) - ETA is "
        << std::put_time(&completion_ETA_tm, "%c")
        << std::endl;
        
    }// pseudo datasets
    
    outfile.close();
    
}

std::vector<std::size_t> clsCMulatorABC::CoverageSingleDataset( std::size_t dataset_idx, std::size_t start_idx, std::size_t end_idx )
{
    std::vector<std::size_t> local_fitness_underestimate_count_vec(_num_alleles, 0);
    
    // process
    for ( auto current_scaled_sim=start_idx; current_scaled_sim<end_idx; ++current_scaled_sim ) {
        
        for ( auto current_allele=0; current_allele<local_fitness_underestimate_count_vec.size(); ++current_allele ) {
            if ( _simulation_result_vector[current_scaled_sim].fitness_values[current_allele] <
                _simulation_result_vector[dataset_idx].fitness_values[current_allele] ) {
                ++local_fitness_underestimate_count_vec[current_allele];
            }
        }
    }

    return local_fitness_underestimate_count_vec;
}
