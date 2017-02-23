
#include "clsCMulatorABC.h"

std::vector<SimulationResult> clsCMulatorABC::RunFitnessInferenceBatch( std::size_t num_simulations )
{
    // initialization
    CMulator local_sim_object(_zparams);
    
    if ( !local_sim_object.IsAbleToInferFitness() ) {
        std::cerr << "Not enough parameters to infer fitness" << std::endl;
        throw "Not enough parameters to infer fitness";
    }
    
    local_sim_object.SetGenerationShift(_actual_data_file.GetFirstGeneration());
    
    // set initial frequencies from actual data
    auto init_freq_vec = _actual_data_file.GetInitFreqs();
    for (auto i = 0; i < init_freq_vec.size(); ++i) {
        local_sim_object.SetAlleleInitFreq(i, init_freq_vec[i]);
    }
    
    // identify wt
   // auto wt_allele_it = std::max_element(init_freq_vec.begin(), init_freq_vec.end());
    //auto wt_allele_idx = static_cast<unsigned int>(std::distance(init_freq_vec.begin(), wt_allele_it));
    // std::cout << " Warning: initializing allele frequency in code and not in parameters. To be fixed. Min=0.0, Max=2.0" << std::endl;
    //local_sim_object.SetWTAllele(wt_allele_idx, 0.0f, 2.0f);
    //std::cout << " WT allele set to " << wt_allele_idx << std::endl;
    
    
    if ( _zparams.GetString("_prior_distribution", "empty") == "uniform" ) {
        _prior_type = UNIFORM;
    }
    
    if ( _zparams.GetString("_prior_distribution", "empty") == "fitness_composite" ) {
        _prior_type = FITNESS_COMPOSITE;
    }
    
    PriorSampler<float> sampler(local_sim_object.GetAlleleMinFitnessValues(),
                                local_sim_object.GetAlleleMaxFitnessValues(),
                                _prior_type);
    sampler.SetManualCategoryProportions(_zparams);
    
    
    //auto fitness_vector_list = _fitness_range.GetRandomCombinations(num_simulations, false);
    auto fitness_vector_list = sampler.SamplePrior(num_simulations, local_sim_object.GetAlleleNumber());
    
    if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
        std::cout << "BEGIN Debug: Prior distribution" << std::endl;
        std::cout << "================================" << std::endl;
        
        for ( auto current_fitness_vector : fitness_vector_list ) {
            
            for ( auto current_fitness_val : current_fitness_vector ) {
                
                std::cout << current_fitness_val << "\t";
            }
            
            std::cout << std::endl;
        }
        
        std::cout << "END Debug: Prior distribution" << std::endl;
        std::cout << "==============================" << std::endl;
    }
    
    std::vector<SimulationResult> tmp_res_vector;
    
    
    // simulation for each set of parameters
    for (auto current_fitness_vector : fitness_vector_list) {
        
        
        local_sim_object.Reset_Soft();
        local_sim_object.SetFitnessValues(current_fitness_vector);
        local_sim_object.EvolveAllGenerations();
        
        _prior_archive.push_back( current_fitness_vector );
        
        SimulationResult sim_result(local_sim_object);
        
        
        
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
        
        
        if ( _rejection_threshold > 0.0  && _use_rejection_threshold ) {
            
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


