//
//  CMulator_Initialization.cpp
//  fits
//
//  Created by Tal Zinger on 15/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "CMulator.h"

/*****************************************************************
 Initialization of simulator variables
 
 Each function throws an exception if its parameters are missing.
 *****************************************************************/

void CMulator::InitMutationRates( ZParams zparams, bool available_mutrate_range )
{
    for (auto from_allele=0; from_allele<_num_alleles; from_allele++) {
        
        FLOAT_TYPE sanity_mutation_rate_sum = 0.0;
        
        for ( auto to_allele=0; to_allele<_num_alleles; to_allele++ ) {
            
            std::string current_mutation_rate_str = PARAM_MUTATION_RATE + std::to_string(from_allele) + "_" + std::to_string(to_allele);
            
            auto tmp_mutation_rate = zparams.GetFloat( current_mutation_rate_str, PARAM_DEFAULT_VAL_FLOAT );
            
            if ( !IsValid_Frequency( tmp_mutation_rate ) && !available_mutrate_range ) {
                std::cerr << _uid << ": Mutation rate from " << from_allele << " to " << to_allele <<
                " is not valid or missing : " << tmp_mutation_rate;
                
                throw _uid + ": Mutation rate from " +
                std::to_string(from_allele) + " to " +
                std::to_string(to_allele) + " is not valid or missing : " + std::to_string( tmp_mutation_rate );
            }
            
            _mutation_rates[from_allele][to_allele] = tmp_mutation_rate;
            
            // side by side while moving to ublast
            _mutation_rates_matrix.at_element(from_allele, to_allele) = tmp_mutation_rate;
            
            sanity_mutation_rate_sum += tmp_mutation_rate;
        }
        
        
        if ( std::fabs(1.0 - sanity_mutation_rate_sum ) > _epsilon_float_compare ) {
            std::cerr << _uid << ": Mutation rates don't sum up to 1.0: " << sanity_mutation_rate_sum
            << ". Delta is " << 1.0 - sanity_mutation_rate_sum << " and epsilon is " << _epsilon_float_compare << std::endl;
            throw( _uid + ": Mutation rates don't sum up to 1.0: " + std::to_string(sanity_mutation_rate_sum) );
        }
    }
}


// 2016-12-22 - Note that now this is log-uniform - i.e. the data read is the power x (as in 10^x)
// 20170116 why should we mention mutation 0_0 in this format? it turns out quite tricky
void CMulator::InitMutationInferenceVariables( ZParams zparams )
{
    
    for (auto from_allele=0; from_allele<_num_alleles; from_allele++) {
        
        for (auto to_allele=0; to_allele<_num_alleles; to_allele++) {
            
            std::string current_min_mutation_rate_str = PARAM_MIN_LOG_MUTATION_RATE + std::to_string(from_allele) + "_" + std::to_string(to_allele);
            
            //std::cout << "attempting to read rate: " << current_min_mutation_rate_str << std::endl;
            auto current_min_mutation_rate_flt = zparams.GetFloat( current_min_mutation_rate_str );
            /*if ( current_min_mutation_rate_flt < 0 ) {
                std::cerr << "Mutation rate inference: missing min mutation rate from " << from_allele << " to " << to_allele << std::endl;
                throw "Mutation rate inference: missing min mutation rate from " + std::to_string(from_allele) + " to " + std::to_string(to_allele);
            }*/
            _min_mutation_rate_matrix(from_allele, to_allele) = current_min_mutation_rate_flt;
            
            
            std::string current_max_mutation_rate_str = PARAM_MAX_LOG_MUTATION_RATE + std::to_string(from_allele) + "_" + std::to_string(to_allele);
            
            auto current_max_mutation_rate_flt = zparams.GetFloat( current_max_mutation_rate_str );
            /*if ( current_max_mutation_rate_flt < 0 ) {
                std::cerr << "Mutation rate inference: missing max mutation rate from " << from_allele << " to " << to_allele << std::endl;
                throw "Mutation rate inference: missing max mutation rate from " + std::to_string(from_allele) + " to " + std::to_string(to_allele);
            }
             */
            _max_mutation_rate_matrix(from_allele, to_allele) = current_max_mutation_rate_flt;
            
        }
    }
}


void CMulator::InitBasicVariables( ZParams zparams )
{
    /* Mandatory Parameters - No Defaults */
    
    try {
        if ( zparams.GetInt( "Debug", 0 ) > 0 ) {
            _debug_mode = true;
        }
        
        // _N = Parameters::getInt( PARAM_POPULATION_SIZE, PARAM_DEFAULT_VAL_INT );
        _N = zparams.GetInt(PARAM_POPULATION_SIZE);
        
        if ( !IsValid_PopulationSize(_N)  ) {
            throw "Parameters: Missing or invalid N (must be positive).";
        }
        
        _alt_N = zparams.GetInt( PARAM_ALT_POPULATION_SIZE, 0 );
        _alt_generation = zparams.GetInt( PARAM_ALT_GENERATION, -1 );
        
        _repeats = zparams.GetInt( PARAM_SIM_REPEATS, 1 );
        
        _N0 = _N;
        
        // _num_alleles = Parameters::getInt( PARAM_NUM_ALLELES, PARAM_DEFAULT_VAL_INT );
        _num_alleles = zparams.GetInt(PARAM_NUM_ALLELES);
        
        if (  !IsValid_NumAlleles(_num_alleles) ) {
            throw "Parameters: Missing or invalid number of alleles (must be >=2)";
        }
        
        //_num_generations = Parameters::getInt( PARAM_NUM_GENERATIONS, PARAM_DEFAULT_VAL_INT );
        _num_generations = zparams.GetInt(PARAM_NUM_GENERATIONS);
        
        //_sample_size = Parameters::getInt( PARAM_SAMPLE_SIZE, PARAM_DEFAULT_VAL_INT );
        _sample_size = zparams.GetInt( PARAM_SAMPLE_SIZE, PARAM_DEFAULT_VAL_INT );
        
        //_bottleneck_interval = Parameters::getInt( PARAM_BOTTLENECK_INTERVAL, PARAM_DEFAULT_VAL_INT );
        _bottleneck_interval = zparams.GetInt( PARAM_BOTTLENECK_INTERVAL, PARAM_DEFAULT_VAL_INT );
        
        //_bottleneck_size = Parameters::getInt( PARAM_BOTTLENECK_SIZE, PARAM_DEFAULT_VAL_INT );
        _bottleneck_size = zparams.GetInt( PARAM_BOTTLENECK_SIZE, PARAM_DEFAULT_VAL_INT );
    }
    catch (const char* exp_txt) {
        std::cerr << "Exception while reading parameters: " << exp_txt << std::endl;
        throw "Can't continue initializtion.";
    }
    catch (std::exception& exp) {
        std::cerr << "Exception while reading parameters: " << exp.what() << std::endl;
        throw "Can't continue initializtion.";
    }
    catch (...) {
        std::cerr << "Unknown exception while reading essential parameters." << std::endl;
        throw "Can't continue initializtion.";
    }
    
    /* Not Mandatory, General Parameters */
    _epsilon_float_compare = zparams.GetFloat( PARAM_EPSILON_FLOAT_COMPARE, PARAM_EPSILON_FLOAT_COMPARE_DEFAULT );
    
    _uid = zparams.GetString( PARAM_SIM_ID, std::string("SimX") );
    
    _logistic_growth = ( zparams.GetInt( PARAM_LOGISTIC_GROWTH, PARAM_DEFAULT_VAL_INT ) > 0 );
    _logistic_growth_K = zparams.GetInt( PARAM_LOGISTIC_GROWTH_K, PARAM_DEFAULT_VAL_INT );
    _logistic_growth_r = zparams.GetFloat( PARAM_LOGISTIC_GROWTH_r, PARAM_DEFAULT_VAL_FLOAT );
    _logistic_growth_t = 0;
    
    _skip_stochastic_step = ( zparams.GetInt(PARAM_SKIP_STOCHASTIC_STEP, PARAM_SKIP_STOCHASTIC_STEP_DEFAULT) > 0);
    
    // seed already initialized, replace only if manual seed is given
    auto original_seed = _time_for_seeding;
    _time_for_seeding = zparams.GetUnsignedInt(PARAM_MANUAL_SEED, original_seed);
    
    
}



void CMulator::InitGeneralInferenceVariables( ZParams zparams )
{
    _generation_shift = zparams.GetInt( PARAM_GENERATION_SHIFT, PARAM_DEFAULT_VAL_INT );
    
    //_fitness_increment = zparams.GetFloat( PARAM_ALLELE_FITNESS_INCRENEMT, PARAM_ALLELE_FITNESS_INCRENEMT_DEFAULT);
    
    _top_percent_to_keep = zparams.GetFloat( PARAM_ACCEPTANCE_RATE, 0.0 ) * 100.0f;
    
    // 2017-01-15 can't allow this to happen that we don't have the percent.
    /* TODO: return to this thought when the CMulator object knows what is its purpose - simulate or inference.
     till then - this is passed to inference section */
    /*try {
        
    }
    catch (...) {
        std::cerr << "Error: cannot get top percent to keep in parameter file." << std::endl;
        throw "Error: cannot get top percent to keep in parameter file.";
    }
    */
    
    //_top_sims_to_keep = zparams.GetUnsignedInt( PARAM_TOP_SIMS_TO_KEEP, TOP_SIMS_TO_KEEP_DEFAULT );
    
    _top_sims_to_keep = _top_percent_to_keep * _repeats / 100;
    
    _max_generations_till_fixation = 1000000;
    
    
    // TODO: make this obsolete, convert to actual file to take this data from
    // auto tmp_no_init_freqs_as_parameters = Parameters::getInt( PARAM_ALLELE_NO_INIT_FREQS, 0 );
    auto tmp_no_init_freqs_as_parameters = zparams.GetInt( PARAM_ALLELE_NO_INIT_FREQS, 1 );
    _no_init_freqs_as_parameters = ( tmp_no_init_freqs_as_parameters == 1 );
    
}


void CMulator::InitInitialAlleleFreqs( ZParams zparams )
{
    FLOAT_TYPE sanity_allele_init_freq_sum = 0.0;
    
    FLOAT_TYPE highest_freq = -1.0;
    
    for (auto current_allele_num=0; current_allele_num<_num_alleles; current_allele_num++) {
        
        std::string current_allele_freq_str = PARAM_ALLELE_INIT_FREQ + std::to_string(current_allele_num);
        
        _allele_init_freqs[current_allele_num] = zparams.GetFloat(current_allele_freq_str, PARAM_DEFAULT_VAL_FLOAT);
        
        // this parameter was not found, were we expecting it to be there?
        if ( _allele_init_freqs[current_allele_num] == PARAM_DEFAULT_VAL_FLOAT ) {
            if (!_no_init_freqs_as_parameters) {
                std::cerr << "CMulator Init Params: no init freqs defined for allele " << current_allele_num << ". To avoid init freqs in parameter file set _no_init_freqs_as_parameters to 1." << std::endl;
                std::string tmp_str = std::string("CCMulator Init Params: no init freqs defined for allele ") +
                std::to_string(current_allele_num) + std::string(". To avoid init freqs in parameter file set _no_init_freqs_as_parameters to 1.");
                throw tmp_str.c_str();
            }
        }
        sanity_allele_init_freq_sum += _allele_init_freqs[current_allele_num];
        
        _sim_data[0][current_allele_num] = _allele_init_freqs[current_allele_num];
        _sim_data_matrix(0, current_allele_num) = _allele_init_freqs[current_allele_num];
        
        if ( highest_freq < _sim_data[0][current_allele_num] && _infer_wt_automatically ) {
            highest_freq = _sim_data[0][current_allele_num];
            _wt_allele_index = current_allele_num;
        }
        
    } // end fitness and init freq
    
    if ( _infer_wt_automatically && _wt_allele_index < 0 ) {
        std::cerr << _uid << ": WT allele not found. Highest frequency is " << highest_freq << std::endl;
        throw( _uid + ": WT allele not found. Highest frequency is " + std::to_string(highest_freq) );
    }
    
    if ( std::fabs( 1.0 - sanity_allele_init_freq_sum ) > _epsilon_float_compare ) {
        std::cerr << _uid << ": Allele initial frequencies don't sum up to 1.0: " << sanity_allele_init_freq_sum << std::endl;
        throw( _uid + ": Allele initial frequencies don't sum up to 1.0: " + std::to_string(sanity_allele_init_freq_sum) );
    }
}


void CMulator::InitFitnessValues( ZParams zparams )
{
    
    for (auto current_allele_num=0; current_allele_num<_num_alleles; current_allele_num++) {
        
        std::string current_allele_fitness_str = PARAM_ALLELE_FITNESS + std::to_string(current_allele_num);
        std::string current_allele_min_fitness = PARAM_ALLELE_MIN_FITNESS + std::to_string(current_allele_num);
        std::string current_allele_max_fitness = PARAM_ALLELE_MAX_FITNESS + std::to_string(current_allele_num);
        
        _allele_fitness[current_allele_num] = zparams.GetFloat(current_allele_fitness_str, -1.0);
        
        if ( _allele_fitness[current_allele_num] < 0 ) {
            //std::cerr << "Missing fitness value for allele " << current_allele_num << std::endl;
            throw "Missing fitness value for allele " + std::to_string(current_allele_num);
        }
    }
    
}


void CMulator::InitFitnessInference( ZParams zparams )
{

    for (auto current_allele_num=0; current_allele_num<_num_alleles; current_allele_num++) {
        
        std::string current_allele_min_fitness = PARAM_ALLELE_MIN_FITNESS + std::to_string(current_allele_num);
        std::string current_allele_max_fitness = PARAM_ALLELE_MAX_FITNESS + std::to_string(current_allele_num);
        
        _allele_min_fitness[current_allele_num] = zparams.GetFloat(current_allele_min_fitness, -1.0);
        if ( _allele_min_fitness[current_allele_num] < 0.0 ) {
            //std::cerr << "Missing minimum fitness value for allele " << current_allele_num << std::endl;
            throw "Missing minimum fitness value for allele " + std::to_string(current_allele_num);
        }
        
        
        _allele_max_fitness[current_allele_num] = zparams.GetFloat(current_allele_max_fitness, -1.0);
        if ( _allele_max_fitness[current_allele_num] < 0.0 ) {
            //std::cerr << "Missing maximum fitness value for allele " << current_allele_num << std::endl;
            throw "Missing maximum fitness value for allele " + std::to_string(current_allele_num);
        }
        
    }
}


void CMulator::InitFixationVariables( ZParams zparams )
{
    _fixation_upper_threshold = zparams.GetFloat("_fixation_upper_threshold", -1.0 );
    
    _fixation_lower_threshold = zparams.GetFloat("_fixation_lower_threshold", -1.0 );
    
    _fixation_consecutive_generations = zparams.GetFloat("_fixation_consecutive_generations", -1.0 );
}

void CMulator::InitPopulationSizeInference( ZParams zparams )
{
    _Nlog_min = zparams.GetInt("_Nlog_min");
    _Nlog_max = zparams.GetInt("_Nlog_max");

}


void CMulator::InitMemberVariables( ZParams zparams )
{
    _available_mutrate = true;
    _available_fitness = true;
    _available_popsize = true;
    _available_essential = true;
    _available_inference = true;
    _available_initfreqs = true;
    _available_fitness_range = true;
    _available_mutrate_range = true;
    _available_popsize_range = true;
    
    // exception here will crash the simulator, this is sort of fine with me
    //std::cout << "\tBasic" << std::endl;
    InitBasicVariables(zparams);
    
    // I want to make sure these are initialized
    // only possible after basic data was read
    _mutation_rates.resize(boost::extents[_num_alleles][_num_alleles]);
    _mutation_rates_matrix.resize(_num_alleles, _num_alleles );
    _wt_allele_index = -1;
    
    // generations + 1 because 0 is given and is not the first processed generation.
    _sim_data.resize(boost::extents[_num_generations+1][_num_alleles]);
    
    _sim_data_matrix.resize( _num_generations+1, _num_alleles );
    _allele_init_freqs.resize(_num_alleles);
    _allele_fitness.resize(_num_alleles);
    _allele_fitness_bar.resize(_num_alleles);
    _allele_fitness_adjusted.resize(_num_alleles);
    
    _min_mutation_rate_matrix.resize(_num_alleles, _num_alleles);
    _max_mutation_rate_matrix.resize(_num_alleles, _num_alleles);
    
    _allele_fitness.resize(_num_alleles);
    _allele_max_fitness.resize(_num_alleles);
    _allele_min_fitness.resize(_num_alleles);
    
    _infer_wt_automatically = (zparams.GetFloat("_infer_wt_automatically", -1.0) > 0);
    _infer_mutation_rate = (zparams.GetFloat(PARAM_INFER_MUTATION_RATE_FLAG, -1.0) > 0);
    
    
    try {
        InitGeneralInferenceVariables(zparams);
    }
    catch (...) {
        _available_essential = false;
        std::cerr << "Error: essential inference parameter missing." << std::endl;
        throw "Error: essential simulation parameter missing.";
    }
    
    
    try {
        InitMutationInferenceVariables(zparams);
    }
    catch (...) {
        _available_mutrate_range = false;
    }
    
    try {
        InitMutationRates(zparams, _available_mutrate_range);
    }
    catch (...) {
        _available_mutrate = false;
    }
    
    
    try {
        InitPopulationSizeInference(zparams);
    }
    catch (...) {
        _available_popsize_range = false;
    }
    
    if ( !_no_init_freqs_as_parameters ) {
        try {
            //std::cout << "\tInitInitialAlleleFreqs" << std::endl;
            InitInitialAlleleFreqs(zparams);
        }
        catch (...) {
            _available_initfreqs = false;
            std::cerr << "Error: Initial allele frequencies expected but not found." << std::endl;
            throw "Error: Initial allele frequencies expected but not found.";
        }
    }
    
    
    try {
        InitFitnessValues(zparams);
    }
    catch (...) {
        _available_fitness = false;
    }
    
    try {
        InitFitnessInference(zparams);
    }
    catch (...) {
        _available_fitness_range = false;
    }
    
    
    /* Internal State */
    _current_generation = 0;
    
    _initialized_with_parameters = true;
}


