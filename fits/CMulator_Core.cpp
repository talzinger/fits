//
//  CMulator_Core.cpp
//  fits
//
//  Created by Tal Zinger on 15/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "CMulator.h"

int CMulator::EvolveToGeneration( int target_generation )
{
    
    // this distribution is unique - it takes a closed range [inclusive]
    if ( _debug_mode ) {
        std::cout << "Initializing local int distribution for alt N to max generation " << target_generation-1 << std::endl;
    }
    boost::random::uniform_int_distribution<int> local_int_distrib(1, target_generation-1);

    if ( _alt_N > 0 && _alt_generation < 0) {
        _alt_generation = local_int_distrib(_boost_gen);
    }
    
    if (!_initialized_with_parameters) {
        std::cerr << " EvolveToGeneration: object not initialized with parameters." << std::endl;
        throw " EvolveToGeneration: object not initialized with parameters.";
    }
    
    if ( _allele_init_freqs.empty() ) {
        std::cerr << " EvolveToGeneration: initial frequencies not provided." << std::endl;
        throw " EvolveToGeneration: initial frequencies not provided.";
    }
    
    if ( _sim_data.empty() ) {
        std::cerr <<	" EvolveToGeneration: sim_data uninitialized." << std::endl;
        throw " EvolveToGeneration: sim_data uninitialized.";
    }
    
    if ( _allele_fitness.empty() ) {
        if ( _allele_max_fitness.empty() || _allele_min_fitness.empty() ) {
            std::cerr <<	" EvolveToGeneration: fitness data empty but so are min or max." << std::endl;
            throw " EvolveToGeneration: fitness data empty but so are min and max.";
        }
    }
    
    // just to inform the calling function what was the initial generation
    int init_generation = _current_generation;
    
    
    // MDOUBLE currentP = 0.0;
    //	MDOUBLE newvals[4] = {0,0,0,0};
    std::vector<FLOAT_TYPE> newvals(_num_alleles);
    
    // 2016-05-16 don't see why this is int and not MDOUBLE
    // int currentSum = 0;
    FLOAT_TYPE currentSum = 0.0;
    
    // should we revert to original population size
    bool alt_popsize_flag = false;
    
    // generation 0 is given, evolve until desired generation
    for (_current_generation=1; _current_generation<=target_generation; _current_generation++) {
        
        if ( _alt_generation == _current_generation ) {
            _old_N = GetPopulationSize();
            SetPopulationSize(_alt_N);
            alt_popsize_flag = true;
            
            if (_debug_mode) {
                std::cout << "Generation "
                << _current_generation
                << ": changed to alt popsize "
                << _N << " instead of "
                << _old_N << std::endl;
            }
        }
        
        if ( alt_popsize_flag ) {
            SetPopulationSize(_old_N);
            alt_popsize_flag = false;
            
            if (_debug_mode) {
                std::cout << "Generation "
                << _current_generation
                << ": reverted popsize to "
                << _N << std::endl;
            }
        }
        
        // normalize fitness
        FLOAT_TYPE wBar = 0.0;
        for (auto i=0; i<_num_alleles; i++ ) {
            
            _allele_fitness_bar[i] = _allele_fitness[i] * _sim_data[_current_generation - 1][i];
            wBar+= _allele_fitness_bar[i];
        }
        //FLOAT_TYPE wBar = 1;
        for (auto i=0; i<_num_alleles; i++ ) {
            _allele_fitness_adjusted[i] = _allele_fitness[i] / wBar;
        }
        
        
        // TODO: replace these nested for loops with matrix multiplication
        // vector of last generation X mutation_rate_matrix -> tmp1_vector
        // tmp1_vector X fitness_vector -> new_generation
        for ( auto currenly_evolving_allele=0; currenly_evolving_allele<_num_alleles; currenly_evolving_allele++ ) {
            
            _sim_data[_current_generation][currenly_evolving_allele] = 0;
            _sim_data_matrix(_current_generation, currenly_evolving_allele) = 0.0f;
            
            
            for ( auto currently_other_allele=0; currently_other_allele<_num_alleles; currently_other_allele++ ) {
                
                // side by side migration to ublast
                /*auto using_old_mutation_rate_matrix =
                _mutation_rates[currently_other_allele][currenly_evolving_allele]
                * _sim_data[_current_generation-1][currently_other_allele]
                * _allele_fitness_adjusted[currently_other_allele];*/
                
                auto using_new_mutation_rate_matrix =
                _mutation_rates_matrix(currently_other_allele,currenly_evolving_allele)
                * _sim_data[_current_generation-1][currently_other_allele]
                * _allele_fitness_adjusted[currently_other_allele];
                

                // make sure mutation rates are the same
                /*if ( using_old_mutation_rate_matrix != using_new_mutation_rate_matrix ) {
                    std::cerr << "from " << currenly_evolving_allele << " to" << currently_other_allele << std::endl;
                    throw "Mutation rates old and new do not agree: " + std::to_string(using_old_mutation_rate_matrix) + " " + std::to_string(using_new_mutation_rate_matrix) ;
                }*/
                
                _sim_data[_current_generation][currenly_evolving_allele] += using_new_mutation_rate_matrix;
                
                /*_sim_data[_current_generation][currenly_evolving_allele] +=
                 _mutation_rates[currently_other_allele][currenly_evolving_allele]
                 * _sim_data[_current_generation-1][currently_other_allele]
                 * _allele_fitness_adjusted[currently_other_allele]; */
            }
        }
        
        // normalize probabilities
        FLOAT_TYPE sumall = 0.0;
        
        for ( auto j=0; j<_num_alleles; j++ ) {
            sumall += _sim_data[_current_generation][j];
        }
        for ( auto j=0; j<_num_alleles; j++ ) {
            _sim_data[_current_generation][j] = _sim_data[_current_generation][j] / sumall;
        }
        
        
        // logistic growth
        if (_logistic_growth) {
            _N = static_cast<int>( _logistic_growth_K *(_N0/(_N0+(exp(-_logistic_growth_r*static_cast<FLOAT_TYPE>(_logistic_growth_t))*(_logistic_growth_K-_N0)))) );
            _logistic_growth_t++;
            
        }
        
        
        // now apply stochastic step
        currentSum = 0;
        
        // randomally select individuals possessing each allele
        if (!_skip_stochastic_step)
        {
            for (auto allele_i = 0; allele_i<_num_alleles; allele_i++) {
                
                // from float to MDOUBLE
                auto currentP = _sim_data[_current_generation][allele_i];
                //currentP = (MDOUBLE)_sim_data[_current_generation][allele_i];
                
                // this returns number of actual individuals chosen - practically an integer
                //newvals[allele_i] = static_cast<int>(bnldev(currentP, _N, time_for_seeding));
                
                boost::random::binomial_distribution<int> local_bin_distrib( _N, currentP );
                newvals[allele_i] = local_bin_distrib(_boost_gen);
                
                // boost::random::binomial_distribution<int> dist( );
                
                
                // actually selected individuals, as each is sampled separately, total samples will probably be != N
                currentSum += newvals[allele_i];
            }
            
            // Normalize frequencies
            for ( auto allele_i = 0; allele_i<_num_alleles; allele_i++ ) {
                _sim_data[_current_generation][allele_i] = static_cast<FLOAT_TYPE>(newvals[allele_i] / currentSum);
            }
        }
        
        
        // apply bottleneck if relevant
        if ( (_bottleneck_interval) > 0 &&
            (_current_generation % _bottleneck_interval == 0) ) {
            
            
            _logistic_growth_t = 0;
            //_N = _bottleneck_size; // commented 2016-03-08 - this caused the population to irreversibly shrink. no wonder drift seemed so strong.
            _N0 = _bottleneck_size;	// this should be like that as this is used when logistic growth is in use.
            
            
            // like stochastic step but with other population size
            // brackets remains from multithreaded period where mutex was used via lock guard (has to be scoped)
            currentSum = 0;
            {
                // randomally select individuals possessing each allele
                for ( auto allele_i = 0; allele_i < _num_alleles; allele_i++ ) {
                    
                    // from float to MDOUBLE
                    //currentP = static_cast<FLOAT_TYPE>(_sim_data[_current_generation][allele_i]);
                    //currentP = (MDOUBLE)_sim_data[_current_generation][allele_i];
                    auto currentP = _sim_data[_current_generation][allele_i];
                    
                    
                    // this returns number of actual individuals chosen - practically an integer
                    //newvals[allele_i] = static_cast<int>(bnldev(currentP, _bottleneck_size, time_for_seeding));
                    
                    boost::random::binomial_distribution<int> local_bin_distrib( _bottleneck_size, currentP );
                    newvals[allele_i] = local_bin_distrib(_boost_gen);
                    
                    // actually selected individuals, as each is sampled separately, total samples will probably be != N
                    currentSum += newvals[allele_i];
                }
            }
            
            // Normalize frequencies
            for ( auto allele_i=0; allele_i<_num_alleles; allele_i++ ) {
                _sim_data[_current_generation][allele_i] = static_cast<FLOAT_TYPE>(newvals[allele_i] / currentSum);
            }
            
            // finished processing for this generation, compare to real data if applicable
            
        }
    }
    
    // at least one simulation performed
    _first_sim_in_batch = 0;
    
    return init_generation;
}


int CMulator::EvolveAllGenerations()
{
    if (!_initialized_with_parameters) {
        throw " EvolveAllGenerations: object not initialized with parameters.";
    }
    
    return EvolveToGeneration(_num_generations);
}


// TODO: if parameter does not exist getint returns 0. how to overcome this when getting parameter to say if we want indeed to evolve until fixation??
int CMulator::EvolveUntilFixation( int tested_allele )
{
    
    if (!_initialized_with_parameters) {
        throw " EvolveUntilFixation: object not initialized with parameters.";
    }
    
    int _consecutive_exceeding_generations = 0;
    int processed_generations = 0;
    
    if ( (tested_allele < 0) || (tested_allele >= _num_alleles) ) {
        std::cerr << "EvolveUntilFixation - invalid number of alleles" << std::endl;
        throw "EvolveUntilFixation - invalid number of alleles";
    }
    
    
    // we need space for two generations
    if ( _num_generations < 2 ) {
        _num_generations = 2;
        _sim_data.resize(boost::extents[2][_num_alleles]);
        Reset_Soft(); // this will update init values
    }
    
    while ( _consecutive_exceeding_generations < _fixation_consecutive_generations ) {
        
        if ( _sim_data[_current_generation][tested_allele] < _fixation_lower_threshold
            || _sim_data[_current_generation][tested_allele] > _fixation_upper_threshold) {
            _consecutive_exceeding_generations++;
        }
        else {
            _consecutive_exceeding_generations = 0;
        }
        
        if ( _current_generation >= MAX_GENERATION ) {
            return -1;
        }
        
        EvolveToGeneration(_current_generation+1);
        
        for ( auto j=0; j<_num_alleles; j++ ) {
            _sim_data[0][j] = _sim_data[1][j];
            _allele_init_freqs[j] = _sim_data[1][j];
        }
        
        _current_generation = 0;
        processed_generations++;
        
        if ( processed_generations > _max_generations_till_fixation ) {
            std::cerr << "Exceeded max generations till fixation: " << processed_generations << " frequency is " << _sim_data[0][tested_allele] << std::endl;
            return -1;
        }
    }
    
    return processed_generations;
}

