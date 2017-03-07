/**********************
 CMulator Core Methods
 **********************/

// #define DEBUG_VERBOSE

#include "clsCMulatorABC.h"


clsCMulatorABC::clsCMulatorABC( ZParams sim_params, ActualDataFile actual_data_file ) :
_zparams(sim_params),
_prior_type(PriorDistributionType::UNIFORM),
_actual_data_file(actual_data_file),
_simulation_result_vector(),
_prior_archive(),
_use_rejection_threshold(true)
{
    ResetRejectionThreshold();
    
    _repeats = _zparams.GetInt( "_num_repeats" );
    
    _num_alleles = _zparams.GetUnsignedInt( "_num_alleles", 0 );
    
    
    auto dummy = _zparams.GetInt( "_top_percent_to_keep", -7 );
    if (dummy>0) {
        std::cerr << "_top_percent_to_keep no longer supported, use _acceptance_rate" << std::endl;
        throw "_top_percent_to_keep no longer supported, use _acceptance_rate";
    }
    
    //_percent_to_keep = _zparams.GetFloat( "_acceptance_rate" ) * 100;
    _sims_to_keep = _zparams.GetFloat( "_acceptance_rate" ) * _repeats;
    
    std::cout << " acceptance " << _zparams.GetFloat( "_acceptance_rate" ) << std::endl;
    std::cout << " keeping " << _sims_to_keep << " simulations" << std::endl;
    
    _selected_actual_generations = actual_data_file.GetActualGenerations();
    _actual_data_raw_freqs = actual_data_file.GetActualFrequencies();
    
    auto _time_for_seeding = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    _boost_gen.seed(_time_for_seeding);
}


void clsCMulatorABC::RunABCInference( FactorToInfer factor, std::size_t number_of_batches )
{
    _simulation_result_vector.clear();
    
    // performace tracking
    boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::mean>> rate_stats;
    
    auto remaining_repeats = _repeats;
    auto repeats_in_batch = _repeats / number_of_batches;
    
    // reserve for 10%, which tends to be the final result count empirically
    _simulation_result_vector.reserve(10 * repeats_in_batch);
    
    std::size_t threshold_update_counter = 0;
    
    double total_running_time = 0.0;
    auto start_global = std::chrono::high_resolution_clock::now();
    
    
    while (remaining_repeats > 0) {
        
        if (repeats_in_batch > remaining_repeats) {
            repeats_in_batch = remaining_repeats;
        }
        remaining_repeats -= repeats_in_batch;
        
        std::vector<SimulationResult> tmp_result_vector;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        _factor_to_infer = factor;
        
        switch (_factor_to_infer) {
            case Fitness:
                tmp_result_vector = RunFitnessInferenceBatch(repeats_in_batch);
                break;
                
            case PopulationSize:
                tmp_result_vector = RunPopulationSizeInferenceBatch(repeats_in_batch);
                break;
                
            case MutationRate:
                tmp_result_vector = RunMutationInferenceBatch(repeats_in_batch);
                break;
        }
        auto end = std::chrono::high_resolution_clock::now();
        
        
        threshold_update_counter += tmp_result_vector.size();
        
        
        std::cout << repeats_in_batch << " simulations -> " << tmp_result_vector.size() << " results." << std::endl;
        
        
        
        // move results from chunk to the main result vector
        //std::size_t batch_start_idx = _simulation_result_vector.size();
        for ( auto &&tmp_result : tmp_result_vector ) {
            
            _simulation_result_vector.push_back(std::move(tmp_result));
            
            /*
            if ( threshold_update_counter > _sims_to_keep ) {
                
                threshold_update_counter = 0;
                
                std::cout << std::endl << "Updating rejection threshold... ";
                
                
                // in order to keep sims_to_keep simulations, we need to order items within [0,sims_to_keep)
                std::nth_element(   _simulation_result_vector.begin(),
                                 _simulation_result_vector.begin() + _sims_to_keep - 1,
                                 _simulation_result_vector.end() );
                
                
                _rejection_threshold = (_simulation_result_vector.begin() + _sims_to_keep - 1)->distance_from_actual;
                
                std::cout << "Done. Value set to " << _rejection_threshold << std::endl;
            }
             */
        }
        //std::size_t batch_end_idx = _simulation_result_vector.size();
        
        //_rejection_threshold = 0.0;
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        //auto batch_sd_vec = GetSDPerAllele(batch_start_idx, batch_end_idx);
        
        //std::cout << "Batch sd values: ";
        //for ( auto tmpval : batch_sd_vec ) { std::cout << "\t" << tmpval; }
        //std::cout << std::endl;
        
        auto sim_count = repeats_in_batch;
        
        auto calculated_speed = static_cast<double>(sim_count) / static_cast<double>(elapsed_ms.count());
        calculated_speed *= 1000; // 1/milisecond to 1/second
        rate_stats(calculated_speed);
        
        // estimate time for completion
        // in seconds
        auto seconds_remaining = remaining_repeats / static_cast<int>(std::round(calculated_speed));
        auto duration_remaining = std::chrono::seconds(seconds_remaining);
        auto current_time = std::chrono::system_clock::now();
        auto completion_ETA = current_time + duration_remaining;
        auto completion_ETA_timet = std::chrono::system_clock::to_time_t(completion_ETA);
        auto completion_ETA_tm = *std::localtime(&completion_ETA_timet);
        
        //std::locale::global(std::locale(""));
        
        std::cout << "Remaining " << remaining_repeats
        << " repeats (rate of "
        << std::round(calculated_speed)
        << "sim/sec) - ETA is "
        << std::put_time(&completion_ETA_tm, "%c")
        << std::endl;
    } // while (remaining_repeats > 0)
    
    if (_simulation_result_vector.size() > _sims_to_keep) {
        
        std::cout << "Total results count is " << _simulation_result_vector.size()
        << ", capacity is " << _simulation_result_vector.capacity() << std::endl;
        
        
        /*
         std::cout << "Deleting unwanted simulation results... ";
        
        std::nth_element(_simulation_result_vector.begin(),
                         _simulation_result_vector.begin() + _sims_to_keep - 1,
                         _simulation_result_vector.end());
        
        _simulation_result_vector.erase(_simulation_result_vector.begin() + _sims_to_keep, _simulation_result_vector.end());
        */
        
        // added 2016-09-04 maybe this will help to prevent page fault previouly caused by large vectors
        _simulation_result_vector.shrink_to_fit();
        //std::cout << "done." << std::endl;
        
        
        std::cout << "Final results count is " << _simulation_result_vector.size() << ", capacity is " << _simulation_result_vector.capacity() << std::endl;
    }
    
    // std::cout << "Sorting results... ";
    
    // std::sort(_simulation_result_vector.begin(), _simulation_result_vector.end());
    
    // std::cout << "Done." << std::endl << std::endl;
    
    auto end_global = std::chrono::high_resolution_clock::now();
    auto global_elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_global - start_global);
    total_running_time = static_cast<double>(global_elapsed_ms.count()) / 1000.0;
    
    std::cout << "Total running time: " << total_running_time << " seconds" << std::endl;
    std::cout << "Average simulation rate: " << boost::accumulators::mean(rate_stats) << " simulations/second" << std::endl;
}


FLOAT_TYPE clsCMulatorABC::ResetRejectionThreshold()
{
    return SetRejectionThreshold(THRESHOLD_RESET_VALUE);
}


FLOAT_TYPE clsCMulatorABC::GetRejectionThreshold()
{
    return _rejection_threshold;
}


FLOAT_TYPE clsCMulatorABC::SetRejectionThreshold(FLOAT_TYPE new_threshold)
{
    FLOAT_TYPE previous_threshold = _rejection_threshold;
    
    // atomic
    //_rejection_threshold.store(new_threshold);
    _rejection_threshold = new_threshold;
    
    return previous_threshold;
}


std::size_t clsCMulatorABC::GetRepeats()
{
    return _repeats;
}

std::vector<SimulationResult> clsCMulatorABC::GetResultsVector(bool only_accepted_results)
{
    if (only_accepted_results) {
        
        std::nth_element(_simulation_result_vector.begin(),
                         _simulation_result_vector.begin() + _sims_to_keep,
                         _simulation_result_vector.end());
        
        std::vector<SimulationResult> tmp_vec( _simulation_result_vector.begin(),
                                               _simulation_result_vector.begin() + _sims_to_keep );
        
        std::cout << "sims " << _sims_to_keep << std::endl;
        std::cout << "tmpvec " << tmp_vec.size() << std::endl;
        
        return tmp_vec;
    }
    
    return _simulation_result_vector;
}

void clsCMulatorABC::SetImmediateRejection(bool new_val)
{
    _use_rejection_threshold = new_val;
}

FLOAT_TYPE clsCMulatorABC::GetDistanceSimActual( ActualDataFile actual_data, SimulationResult sim_result )
{
    // get the generations actual data is available for
    auto generation_list = actual_data.GetActualGenerations();
    auto actual_freqs_matrix = actual_data.GetActualFreqsAsMatrix();
    
    // get matrix indeces equivalent of the desired generations by shifting
    std::vector<int> shifted_generations( generation_list.size(), sim_result.generation_shift );
    
    std::transform( generation_list.cbegin(), generation_list.cend(),
                   shifted_generations.cbegin(),
                   shifted_generations.begin(),
                [](int generation, int shift) { return generation - shift; });
    
    MATRIX_TYPE sim_freqs_matrix( generation_list.size(), sim_result.sim_data_matrix.size2() );
    
    for ( auto current_target_row=0; current_target_row<generation_list.size(); ++current_target_row ) {
        boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<float>> original_row( sim_result.sim_data_matrix, generation_list[current_target_row] );
        boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<float>> target_row( sim_freqs_matrix, current_target_row );
        
        std::cout << "copying row " << generation_list[current_target_row] << " as row " << current_target_row << std::endl;
        
        target_row = original_row;
    }
    
    if ( sim_freqs_matrix.size1() != actual_freqs_matrix.size1() ||
        sim_freqs_matrix.size2() != actual_freqs_matrix.size2() ) {
        std::cerr << "Actual and sim matrices don't match in size: actual("
        << actual_freqs_matrix.size1()
        << ","
        << actual_freqs_matrix.size2()
        << ") sim("
        << sim_freqs_matrix.size1()
        << ","
        << sim_freqs_matrix.size2()
        << ")"
        << std::endl;
        
        throw "Actual and sim matrices don't match in size";
    }
    
    // get the scaling factor
    auto scaling_factor_vec = GetSDPerAllele(0, _simulation_result_vector.size());
    
    // we now have compatible actual_data and simulated_data matrices
    std::vector<FLOAT_TYPE> distance_per_allele(scaling_factor_vec.size(), 0.0f);
    
    for ( auto current_allele=0; current_allele<scaling_factor_vec.size(); ++current_allele ) {
        
    }
    
    // calculate scaled distance
    FLOAT_TYPE ret_val = -0.1;
    return ret_val;
}
