//
//  CMulator_GetSet.cpp
//  CMulatorLib
//
//  Created by Tal Zinger on 31/12/2015.
//  Copyright © 2015 Tal Zinger. All rights reserved.
//

// #define DEBUG_VERBOSE

#include "CMulator.h"


FLOAT_TYPE CMulator::GetSimulatedFrequency( int generation, int allele ) const
{
    if (  !IsValid_Generation(generation) ) {
        std::cerr << "ERROR: illegal generation: " << generation << std::endl;
        throw std::out_of_range("ERROR: illegal generation.");
    }
    
    if ( !IsValid_Allele(allele) ) {
        std::cerr << "ERROR: illegal allele: " << generation << std::endl;
        throw std::out_of_range("ERROR: illegal allele.");
    }
    
    return _sim_data[generation][allele];
}

void CMulator::SetSimulatedFrequency(  int generation, int allele, FLOAT_TYPE frequency )
{
    if (  !IsValid_Generation(generation) ) {
        std::cerr << "ERROR: illegal generation: " << generation << std::endl;
        throw std::out_of_range("ERROR: illegal generation.");
    }
    
    if ( !IsValid_Allele(allele) ) {
        std::cerr << "ERROR: illegal allele: " << allele << std::endl;
        throw std::out_of_range("ERROR: illegal allele.");
    }
    
    if ( !IsValid_Frequency(frequency) ) {
        std::cerr << "ERROR: illegal allele frequency: " << frequency << std::endl;
        throw std::out_of_range("ERROR: illegal allele frequency.");
    }
    
    _sim_data[generation][allele] = frequency;
}

std::vector<FLOAT_TYPE> CMulator::GetAlleleFitnessValues() const
{
	return _allele_fitness;
}


int CMulator::GetNumOfGenerations() const
{
	return _num_generations;
}


void CMulator::SetNumOfGeneration(int generations)
{
	if ( !IsValid_Generation(generations) ) {
		throw " SetNumOfGeneration: invalid value (" + std::to_string(generations) + ").";
	}
	
	_num_generations = generations;
	// UpdateParameters();
	Reset_Soft();
}


std::string CMulator::GetSimUID() const
{
	return _uid;
}


void CMulator::SetSimUID( std::string new_sim_uid )
{
	_uid = new_sim_uid;
	// UpdateParameters();
}


void CMulator::SetBottleneckSize( int size )
{
	if ( !IsValid_BottleneckSize(size) ) {
		throw " SetBottleneckSize: invalid value (" + std::to_string(size) + ").";
	}
	
	_bottleneck_size = size;
}


int CMulator::GetBottleneckSize() const
{
	return _bottleneck_size;
}


std::string CMulator::GetRFileName() const
{
	return _R_filename;
}


void CMulator::SetRFileName( const std::string& filename )
{
	_R_filename = filename;
	// UpdateParameters();
}


void CMulator::SetAlleleFitnessValue( int allele, FLOAT_TYPE fitness )
{
	if ( !IsValid_Fitness(fitness) ) {
		throw " SetAlleleFitnessValue: invalid fitness value (" + std::to_string(fitness) + ").";
	}
	
	if ( !IsValid_Allele(allele) ) {
		throw " SetAlleleFitnessValue: invalid allele (" + std::to_string(allele) + ").";
	}
	
	_allele_fitness[allele] = fitness;
}


void CMulator::SetAlleleInitFreqs( std::vector<FLOAT_TYPE> freqs )
{
	// todo: make sure sum to 1 but with epsilon
	// FLOAT_TYPE tmp_sum = 0.0;
	
	_allele_init_freqs.resize( freqs.size() );
	_allele_init_freqs = freqs;
	
	// UpdateParameters();
	
	// for consistency reasons, we must reset the simulated data
	// this will also update the init freqs
	Reset_Soft();
}

std::vector<FLOAT_TYPE> CMulator::GetAlleleInitFreqs() const
{
	return _allele_init_freqs;	
}


std::vector<FLOAT_TYPE> CMulator::GetAlleleMaxFitnessValues() const
{
	return _allele_max_fitness;
}


std::vector<FLOAT_TYPE> CMulator::GetAlleleMinFitnessValues() const
{
	return _allele_min_fitness;
}


FLOAT_TYPE CMulator::GetAlleleFreq( int generation, int allele ) const
{
	if ( !IsValid_Generation(generation) ) {
		throw "GetAlleleFreq: illegal generation " + std::to_string(generation);
	}
	
	if ( !IsValid_Allele(allele) ) {
		throw "GetAlleleFreq: illegal allele " + std::to_string(allele);
	}
	
	return _sim_data[generation][allele];
}


std::vector<FLOAT_TYPE> CMulator::GetRawFrequencyData( bool resample )
{
	std::vector<FLOAT_TYPE> tmp_freqs;
    
    //std::cout << "Raw begin1...";
	for ( auto cur_generation=0; cur_generation<_num_generations; cur_generation++ ) {
		
        std::vector<FLOAT_TYPE> resampled_freqs(_num_alleles, 0.0f);
        
        FLOAT_TYPE currentSum = 0.0f;
        
        for ( auto allele_i = 0; allele_i<_num_alleles; allele_i++ ) {
            
            
            auto currentP = _sim_data[_current_generation][allele_i];
            
            boost::random::binomial_distribution<int> local_bin_distrib( _sample_size, currentP );
            resampled_freqs[allele_i] = local_bin_distrib(_boost_gen);
            
            currentSum += resampled_freqs[allele_i];
        }
        
        // Normalize frequencies
        for ( auto allele_i = 0; allele_i<_num_alleles; allele_i++ ) {
            resampled_freqs[allele_i] = resampled_freqs[allele_i] / currentSum;
        }
		
		for (int cur_allele=0; cur_allele<_num_alleles; cur_allele++) {
            
            // TODO: this is indeed inefficient, move the if such that the resampling won't happen always
            if ( resample ) {
                tmp_freqs.push_back( resampled_freqs[cur_allele] );
            }
            else {
                tmp_freqs.push_back( _sim_data[cur_generation][cur_allele] );
            }
		}
	}
	
    //std::cout << "Raw end1.";
	return tmp_freqs;
}


std::vector<FLOAT_TYPE> CMulator::GetRawFrequencyData( std::vector<int> selected_generations, bool resample )
{
    std::vector<FLOAT_TYPE> tmp_freqs;
	
    //std::cout << "Raw begin2...";
    
	// sort so the results would be given in a predictable order
    std::sort( selected_generations.begin(), selected_generations.end() );
	
	for ( auto cur_generation : selected_generations ) {
		
		auto tmp_cur_generation = cur_generation - _generation_shift;
		
        std::vector<FLOAT_TYPE> resampled_freqs(_num_alleles, 0.0f);
        
        FLOAT_TYPE currentSum = 0.0f;
        
        for ( auto allele_i = 0; allele_i<_num_alleles; allele_i++ ) {
            
            //auto currentP = _sim_data[tmp_cur_generation][allele_i];
            auto currentP = GetAlleleFreq(tmp_cur_generation, allele_i);
            
            boost::random::binomial_distribution<int> local_bin_distrib( _sample_size, currentP );
            resampled_freqs[allele_i] = local_bin_distrib(_boost_gen);
            
            currentSum += resampled_freqs[allele_i];
        }
        
        // Normalize frequencies
        for ( auto allele_i = 0; allele_i<_num_alleles; allele_i++ ) {
            resampled_freqs[allele_i] = resampled_freqs[allele_i] / currentSum;
        }
        
		for ( auto cur_allele=0; cur_allele<_num_alleles; cur_allele++ ) {
			
			if ( !IsValid_Generation(cur_generation) ) {
                std::cerr << " GetRawFrequencyData: generation out of range (" + std::to_string(tmp_cur_generation) + "). number of generations is " << _num_generations << std::endl;
                throw " GetRawFrequencyData: generation out of range (" + std::to_string(tmp_cur_generation) + ").";
			}
			
			// cout << "getting generation " << tmp_cur_generation << " with freq " << _sim_data[tmp_cur_generation][cur_allele] << " size queue before is " << tmp_freqs.size() << endl;
            
            // TODO: this is indeed inefficient, move the if such that the resampling won't happen always
            if ( resample ) {
                tmp_freqs.push_back( resampled_freqs[cur_allele] );
            }
            else {
                tmp_freqs.push_back( _sim_data[tmp_cur_generation][cur_allele] );
            }
            //tmp_freqs.push_back( _sim_data[tmp_cur_generation][cur_allele] );
		}
	}
	
    //std::cout << "Raw end2.";
	return tmp_freqs;
}


FLOAT_TYPE CMulator::GetMutationRate( int from, int to ) const
{
	if ( from < 0 || from >= _num_alleles ) {
		throw " GetMutationRate: from allele out of range (" + std::to_string(from) + ").";
	}
	
	if ( to < 0 || to >= _num_alleles ) {
		throw " GetMutationRate: from allele out of range (" + std::to_string(from) + ").";
	}
	
    return _mutation_rates_matrix(from,to);
	//return _mutation_rates[from][to];
}


void CMulator::SetMutationRate( int from, int to, FLOAT_TYPE rate )
{
	if ( from < 0 || from >= _num_alleles ) {
		throw " SetMutationRate: from allele out of range (" + std::to_string(from) + ").";
	}
	
	if ( to < 0 || to >= _num_alleles ) {
		throw " SetMutationRate: from allele out of range (" + std::to_string(from) + ").";
	}
	
	if ( rate < 0 || rate > 1.0 ) {
		throw " SetMutationRate: rate of mutation not valid (" + std::to_string(rate) + ").";
	}
	
	// update both parameter and internal variable
	//string param_name = PARAM_MUTATION_RATE + to_string(from) + "_" + to_string(to);
	//Parameters::updateParameter( param_name, to_string(rate).c_str() );
	
	_mutation_rates[from][to] = rate;
    _mutation_rates_matrix(from,to) = rate;
	// UpdateParameters();
}


void CMulator::SetAlleleInitFreq( int allele, FLOAT_TYPE freq )
{
	if ( !IsValid_Allele(allele) ) {
		throw " SetAlleleInitFreq: allele out of range (" + std::to_string(allele) + ").";
	}
	
	if ( !IsValid_Frequency(freq) && freq != PARAM_DEFAULT_VAL_FLOAT ) {
		throw " SetAlleleInitFreq: frequency out of range (" + std::to_string(freq) + ").";
	}
	
	_sim_data[0][allele] = freq;
	_allele_init_freqs[allele] = freq;
	
	// UpdateParameters();
	Reset_Soft();
}


void CMulator::SetAlleleMinFitness( int allele, FLOAT_TYPE fitness )
{
	if ( !IsValid_Allele(allele) ) {
		throw " SetAlleleMinFitness: allele out of range (" + std::to_string(allele) + ").";
	}
	
	if ( !IsValid_Fitness(fitness  && fitness != PARAM_DEFAULT_VAL_FLOAT ) ) {
		throw " SetAlleleMinFitness: fitness invalid (" + std::to_string(fitness) + ").";
	}
	
	_allele_min_fitness[allele] = fitness;
	// UpdateParameters();
}


void CMulator::SetAlleleMaxFitness( int allele, FLOAT_TYPE fitness )
{
	if ( !IsValid_Allele(allele) ) {
		throw " SetAlleleMaxFitness: allele out of range (" + std::to_string(allele) + ").";
	}
	
	if ( !IsValid_Fitness(fitness)  && fitness != PARAM_DEFAULT_VAL_FLOAT ) {
		throw " SetAlleleMaxFitness: fitness invalid (" + std::to_string(fitness) + ").";
	}
	
	_allele_max_fitness[allele] = fitness;
	// UpdateParameters();
}


void CMulator::SetSamplingSize( int sampling_size )
{
	if ( !IsValid_SamplingSize( sampling_size ) ) {
		throw " SetSamplingSize: invalid value (" + std::to_string(sampling_size) + ").";
	}
	
	_sample_size = sampling_size;
	// UpdateParameters();
}


int CMulator::GetSamplingSize() const
{
	return _sample_size;
}


void CMulator::SetBottleneckInterval( int interval )
{
	if ( !IsValid_BottleneckInterval(interval) ) {
		throw " SetBottleneckInterval: invalid value (" + std::to_string(interval) + ").";
	}
	
	_bottleneck_interval = interval;
	// UpdateParameters();
}


int CMulator::GetBottleneckInterval() const
{
	return _bottleneck_interval;
}


/*
void CMulator::SetPopulationSize( int N )
{
	if ( !IsValid_PopulationSize(N) ) {
		throw " SetPopulationSize: invalid value (" + to_string(N) + ").";
	}
	
	_N = N;
	// UpdateParameters();
}
*/

int CMulator::GetPopulationSize() const
{
	return _N;
}


int CMulator::GetGenerationShift() const
{
	if (!_initialized_with_parameters) {
		throw " GetGenerationShift: object not initialized with parameters.";
	}
	
	return _generation_shift;
}


void CMulator::SetGenerationShift(int shift)
{
	if ( shift < 0 ) {
		throw " SetGenerationShift: invalid value (" + std::to_string(shift) + ").";
	}
	
	_generation_shift = shift;
}


int CMulator::GetAlleleNumber() const
{
	if (!_initialized_with_parameters) {
		throw " GetAlleleNumber: object not initialized with parameters.";
	}
	
	return _num_alleles;
}


//void CMulator::SetAlleleNumber(int alleles)
//{
	// This requires a massive reset in the data, so update params and reload
    //_num_alleles = alleles;
    
	//string tmp_str = to_string( alleles );
	//Parameters::updateParameter( PARAM_NUM_ALLELES, tmp_str.c_str() );
	// Reset_Hard();
//}


FLOAT_TYPE CMulator::GetInitAlleleFreq(int allele) const
{
	if (!_initialized_with_parameters) {
		throw " GetInitAlleleFreq: object not initialized with parameters.";
	}
	
	return _allele_init_freqs[allele];
}


int CMulator::GetWTAllele() const
{
	return _wt_allele_index;
}


void CMulator::SetWTAllele(int wt_index, FLOAT_TYPE min_fitness_nonwt, FLOAT_TYPE max_fitness_nonwt)
{
	for (int i = 0; i < _num_alleles; ++i) {
		SetAlleleMaxFitness(i, max_fitness_nonwt);
		SetAlleleMinFitness(i, min_fitness_nonwt);
	}

	SetAlleleMaxFitness(wt_index, 1.0f);
	SetAlleleMinFitness(wt_index, 1.0f);

	_wt_allele_index = wt_index;
}


void CMulator::SetFitnessValues( std::vector<FLOAT_TYPE> _given_allele_fitness )
{
    if ( _given_allele_fitness.size() != _num_alleles ) {
		std::cerr << " SetFitnessValues: num of alleles is " << _num_alleles << " but trying to set " << _given_allele_fitness.size() << " values " << std::endl;
		throw " SetFitnessValues: trying to set fitness with different size than number of alleles";
	}
	
	_allele_fitness = _given_allele_fitness;
    _available_fitness = true;
}


unsigned int CMulator::GetRandomSeed() const
{
	return _time_for_seeding;
}


void CMulator::SetRandomSeed(unsigned int new_seed)
{
	_time_for_seeding = new_seed;
}


unsigned int CMulator::GetTopPercentToKeep() const
{
	return _top_percent_to_keep;
}


void CMulator::SetTopPercentToKeep(unsigned int new_percent)
{
	if ( new_percent > 100 ) {
		std::cerr << "Invalid value for percent to keep: " << std::to_string(new_percent) << std::endl;
		throw "Invalid value for percent to keep: " + std::to_string(new_percent);
	}

    std::cerr << "SetTopPercentToKeep" << std::endl;
	_top_sims_to_keep = -1;
	_top_percent_to_keep = new_percent;
}


unsigned int CMulator::GetTopSimsToKeep() const
{
	return _top_sims_to_keep;
}


void CMulator::SetTopSimsToKeep(unsigned int new_sims)
{
	// can't actually get illegal value as this is unsigend

	_top_percent_to_keep = -1;
	_top_sims_to_keep = new_sims;
}


void CMulator::SetSkipStochasticStep(bool skip)
{
	_skip_stochastic_step = skip;
}


bool CMulator::GetSkipStochasticStep() const
{
	return _skip_stochastic_step;
}


void CMulator::SetMutationRateMatrix( MATRIX_TYPE new_mutation_matrix )
{
    // TODO: check that rates are summed to 1
    // TODO: update some dirty bit
    _mutation_rates_matrix = new_mutation_matrix;
}


MATRIX_TYPE CMulator::GetMutationRateMatrix() const
{
    return _mutation_rates_matrix;
}

MATRIX_TYPE CMulator::GetMinMutationRateMatrix() const
{
    return _min_mutation_rate_matrix;
}

MATRIX_TYPE CMulator::GetMaxMutationRateMatrix() const
{
    return _max_mutation_rate_matrix;
}


int CMulator::GetRepeats() const
{
    return _repeats;
}

void CMulator::SetPopulationSize( int N )
{
    _N = N;
    _available_popsize = true;
}
