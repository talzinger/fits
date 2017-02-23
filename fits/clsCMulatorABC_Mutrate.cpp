
#include "clsCMulatorABC.h"

std::vector<SimulationResult> clsCMulatorABC::RunMutationInferenceBatch( std::size_t num_simulations )
{
    // initialization
    CMulator local_sim_object(_zparams);
    
    if ( !local_sim_object.IsAbleToInferMutationRate() ) {
        std::cerr << "Not enough parameters to infer mutation rate" << std::endl;
        throw "Not enough parameters to infer mutation rate";
    }
    
    local_sim_object.SetGenerationShift(_actual_data_file.GetFirstGeneration());
    
    // set initial frequencies from actual data
    auto init_freq_vec = _actual_data_file.GetInitFreqs();
    for (auto i = 0; i < init_freq_vec.size(); ++i) {
        local_sim_object.SetAlleleInitFreq(i, init_freq_vec[i]);
    }
    
    // identify wt
    auto wt_allele_it = std::max_element(init_freq_vec.begin(), init_freq_vec.end());
    auto wt_allele_idx = static_cast<unsigned int>(std::distance(init_freq_vec.begin(), wt_allele_it));
    // std::cout << " Warning: initializing allele frequency in code and not in parameters. To be fixed. Min=0.0, Max=2.0" << std::endl;
    local_sim_object.SetWTAllele(wt_allele_idx, 0.0f, 2.0f);
    //std::cout << " WT allele set to " << wt_allele_idx << std::endl;
    
    
    std::vector<FLOAT_TYPE> min_mutrates;
    std::vector<FLOAT_TYPE> max_mutrates;
    auto min_matrix = local_sim_object.GetMinMutationRateMatrix();
    auto max_matrix = local_sim_object.GetMaxMutationRateMatrix();
    
    for ( auto row=0; row<local_sim_object.GetAlleleNumber(); ++row ) {
        for ( auto col=0; col<local_sim_object.GetAlleleNumber(); ++col ) {
            
            min_mutrates.push_back( min_matrix(row,col) );
            max_mutrates.push_back( max_matrix(row,col) );
        }
    }
    
    PriorSampler<FLOAT_TYPE> sampler( min_mutrates, max_mutrates, PriorDistributionType::UNIFORM);
    
    auto mutrate_vector_list = sampler.SamplePrior(num_simulations, min_mutrates.size() );
    
    
    std::vector<SimulationResult> tmp_res_vector;
    
    // std::cout << "Prior archive size : " << _prior_archive.size() << std::endl;
    // simulation for each set of parameters
    // 2017-02-07 changed to regular for loop to be able to record index
    //for (auto current_mutrate_vector : mutrate_vector_list) {
    for (auto current_mutrate_idx=0; current_mutrate_idx<mutrate_vector_list.size(); ++current_mutrate_idx ) {
        
        auto current_mutrate_vector = mutrate_vector_list[current_mutrate_idx];
        
        local_sim_object.Reset_Soft();
        
        MATRIX_TYPE tmp_mutrates_matrix( local_sim_object.GetAlleleNumber(), local_sim_object.GetAlleleNumber() );
        
        if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
            std::cout << "Putting sampled mutation rates in matrix." << std::endl;
        }
        
        std::vector<FLOAT_TYPE> tmp_line_sum( tmp_mutrates_matrix.size1(), 0.0f );
        for ( auto i=0; i<current_mutrate_vector.size(); ++i ) {
            
            auto row = i / local_sim_object.GetAlleleNumber();
            auto col = i % local_sim_object.GetAlleleNumber();
            
            // power of the mutation rate
            tmp_mutrates_matrix(row,col) = std::pow( 10, current_mutrate_vector[i] );
         
            // to normalize such that row sums up to 1.0
            if ( row != col ) {
                tmp_line_sum[row] += tmp_mutrates_matrix(row, col);
            }
        }
        
        // normalize
        for ( auto row=0; row<tmp_mutrates_matrix.size1(); ++row ) {
            tmp_mutrates_matrix(row,row) = 1.0f - tmp_line_sum[row];
        }
        
        if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
            std::cout << "Mutation matrix processing finished:" << std::endl;
            std::cout << tmp_mutrates_matrix << std::endl;
        }
        
        local_sim_object.SetMutationRateMatrix(tmp_mutrates_matrix);
        
        local_sim_object.EvolveAllGenerations();
        
        _prior_archive.push_back( current_mutrate_vector );
        
        
        SimulationResult sim_result(local_sim_object);
        
        //sim_result.prior_sample_index = _prior_archive.size() - 1;
        
        
        auto raw_frequency_data = local_sim_object.GetRawFrequencyData(_selected_actual_generations);
        
        if ( raw_frequency_data.empty() ) {
            throw std::range_error("Simulated data vector is empty");
        }
        
        
        // didn't get the correct number of simulated generations/alleles
        if (_actual_data_raw_freqs.size() != raw_frequency_data.size()) {
            
            std::string tmp_str = "Actual data size is "
            + std::to_string(_actual_data_raw_freqs.size())
            + " and simulated data size is "
            + std::to_string(raw_frequency_data.size());
            
            throw std::range_error(tmp_str);
        }
        
        sim_result.distance_from_actual = CMulatorToolbox::GetDistanceSimActual(_actual_data_raw_freqs, raw_frequency_data);
        
        // just to verify the results is the same for debug
        /*sim_result.individual_distance = CMulatorToolbox::GetDistanceVectorSimActual(   _actual_data_raw_freqs,
         raw_frequency_data,
         local_sim_object.GetAlleleNumber() );
         
         sim_result.distance_from_actual = std::accumulate( sim_result.individual_distance.cbegin(), sim_result.individual_distance.cend(), 0.0f );*/
        //assert( tmp == sim_result.distance_from_actual );
        
        if ( _rejection_threshold > 0.0 && _use_rejection_threshold ) {
            
            if ( sim_result.distance_from_actual < _rejection_threshold ) {
                tmp_res_vector.push_back(std::move(sim_result));
            }
        }
        else {
            tmp_res_vector.push_back(std::move(sim_result));
        }
    } // fitness
    
    
    return tmp_res_vector;
    
} // RunFitnessInferenceBatch


