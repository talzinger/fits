//
//  clsFieldRange.cpp
//  NumFieldRange
//
//  Created by Tal Zinger on 13/10/2015.
//  Copyright © 2015 Tal Zinger. All rights reserved.
//

// #define DEBUG_VERBOSE

#include "clsFieldRange.hpp"

// using namespace std;

clsFieldRange::clsFieldRange()
	: _initialized(false)
{}

void InitializeMemberVariables()
{}

// copy constructor
// TODO: only copy members if original is initialized!
clsFieldRange::clsFieldRange( const clsFieldRange &original ) :
_initialized(true),
_rnd_engine(),
_increment(original._increment),
_max_reached(original._max_reached),
_num_of_combinations(original._num_of_combinations),
_num_of_steps_passed(original._num_of_steps_passed),
_range_type(original._range_type)
{
	_rnd_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    _rnd_int_seed = -static_cast<int>(_rnd_seed);
    auto _rnd_unsinged_int_seed = static_cast<unsigned int>(_rnd_seed);
    
	_rnd_engine.seed(_rnd_seed);
    _boost_gen.seed(_rnd_unsinged_int_seed);

	_min_vector = original._min_vector;
	_max_vector = original._max_vector;
	
	_data_vector = original._data_vector;
	
	_is_data_reset = original._is_data_reset;
	
	
	_all_combinations = original._all_combinations;

	// can't copy this, have to reinitialize
	// unless I find a way to properly copy this
	_combination_iterator = _all_combinations.begin();
		
	
	_subrange_min_vector = original._subrange_min_vector;
	_subrange_max_vector = original._subrange_max_vector;
	_subrange_increment = original._subrange_increment;
	
	_use_subrange = original._use_subrange;
	_combinations_with_subrange = original._combinations_with_subrange;

	_subrange_increment = original._subrange_increment;
	_repeats = original._repeats;

	_prior_dist_type = original._prior_dist_type;

	
}

// ctor matrix
clsFieldRange::clsFieldRange( boost::numeric::ublas::matrix<FLOAT_TYPE> min_values, boost::numeric::ublas::matrix<FLOAT_TYPE> max_values, FLOAT_TYPE increment, int repeats, PriorDistributionType prior_dist)
	: _subrange_min_vector(_min_vector.size()),
	_subrange_max_vector(_min_vector.size()),
	_repeats(repeats),
	_initialized(true),
	_all_combinations(),
	_prior_dist_type(prior_dist),
	_combinations_with_subrange(false),
_range_type(RangeType::MATRIX_RANGE)
{
	if (min_values.size1() != max_values.size1()) {
		throw "clsFieldRange: min-max matrices do not match in size1.";
	}
	if (min_values.size2() != max_values.size2()) {
		throw "clsFieldRange: min-max matrices do not match in size2.";
	}
	if (min_values.size1() != min_values.size2()) {
		throw "clsFieldRange: matrices are not square.";
	}

	if (increment <= 0) {
		throw "clsFieldRange: increment must be positive.";
	}

    _rnd_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    _rnd_int_seed = -static_cast<int>(_rnd_seed);
    auto _rnd_unsinged_int_seed = static_cast<unsigned int>(_rnd_seed);
    
    _rnd_engine.seed(_rnd_seed);
    _boost_gen.seed(_rnd_unsinged_int_seed);
    
	// need to make sure sizes are correct, this is why not initializing earlier
	_min_matrix = min_values;
	_max_matrix = max_values;

	//_data_vector = min_values;
	_is_data_reset = true;

	_increment = increment;
	_max_reached = false;
	_num_of_steps_passed = 0;
	_num_of_combinations = 1;
}


//clsFieldRange::clsFieldRange( vector<FLOAT_TYPE> min_values, vector<FLOAT_TYPE> max_values, FLOAT_TYPE increment, int repeats ) 
clsFieldRange::clsFieldRange( std::vector<FLOAT_TYPE> min_values, std::vector<FLOAT_TYPE> max_values, FLOAT_TYPE increment, int repeats, PriorDistributionType prior_dist) :
_subrange_min_vector(_min_vector.size()),
_subrange_max_vector(_min_vector.size()),
_repeats(repeats),
_initialized(true),
_all_combinations(),
_prior_dist_type(prior_dist),
_rnd_seed(0),
_combinations_with_subrange(false),
_range_type(RangeType::VECTOR_RANGE)
{
	if ( min_values.size() != max_values.size() ) {
		throw "clsFieldRange: min-max vectors do not match in size.";
	}
	if ( increment <=0 ) {
		throw "clsFieldRange: increment must be positive.";
	}
	
	// need to make sure sizes are correct, this is why not initializing earlier
	_min_vector = min_values;
	_max_vector = max_values;
	
	_data_vector = min_values;
	_is_data_reset = true;
	
	_increment = increment;
	_max_reached = false;
	_num_of_steps_passed = 0;
	_num_of_combinations = 1;
	
	
	try {
		for ( auto i = 0; i < min_values.size(); i++ ) {
			_num_of_combinations *= static_cast<int>((max_values[i] - min_values[i]) / increment + 1.0f);
		}
	}
	catch (...) {
        std::cerr << "unknown exception in ctor of clsFieldRange" << std::endl;
		std::exit(1);
	}

	_subrange_increment = _increment;
	_use_subrange = false;

	// generate all recombinations
	GenerateAllRecombinations();
	_combination_iterator = _all_combinations.begin();
}


std::vector<std::vector<FLOAT_TYPE>> clsFieldRange::GetAllCombinations()
{
	// don't get dirty data while its being generated
#ifdef ENABLE_PARALLEL
	lock_guard<mutex> combi_lock(clsFieldRange::_COMBI_GENERATION_MUTEX);
#endif

	if (!_initialized) {
        std::cerr << "clsFieldRange: cannot get combinations of uninitialized object" << std::endl;
		throw "clsFieldRange: cannot get combinations of uninitialized object";
	}
	
	// already generated, simply return
	if (_all_combinations.size() > 0) {
		
		return _all_combinations;
	}

	if (_all_combinations.size() == 0) {
		GenerateAllRecombinations();
	}

	// either this was generated by this instance, or the condition failed and reached here (generated by another instance)
	return _all_combinations;
}


bool clsFieldRange::GenerateAllRecombinations()
{
	_all_combinations.clear();

	// could be shorted, wanted to keep ability to catch errors
	while (!_max_reached) {

		bool res = GenerateNextCombination();

		if (!res) {
			return false;
		}
	}

	// cout << "setting iterator" << endl;
	_combination_iterator = _all_combinations.begin();
	return true;
}


bool clsFieldRange::GenerateNextCombination()
{
	bool add_flag = true; // will cause the first digit to increase
	
	FLOAT_TYPE eps = _increment / 2;
	
	// all range covered
	if (_max_reached) {
		return false;
	}
	
		
	for ( auto i=0; i<_data_vector.size(); i++ ) {
		
		if (add_flag) {
			
			if ( _use_subrange ) {
				if ( _data_vector[i] >= _subrange_min_vector[i] && _data_vector[i] < _subrange_max_vector[i] ) {
					_data_vector[i]+=_subrange_increment;
				}
				else {
					_data_vector[i]+=_increment;
				}
			} 
			else {
				_data_vector[i]+=_increment;
			}
			add_flag = false;
		}
		
		//cout << _data_vector[0] << endl;
		
		if (_data_vector[i]>_max_vector[i]+eps) {
			
			// this is the last digit, nowhere to go
			if (i==_data_vector.size()-1) {
				_data_vector = _max_vector;
				_max_reached = true;
				return false;
			}
			
			_data_vector[i] = _min_vector[i];
			add_flag = true;
		}
		
	}
	
	// will be positive if we could not actually increase value
	//_max_reached = add_flag;
	
	// check if we have reached the maximum
	bool tmpflag_maxreached = true;

	for (auto i = 0; i < _min_vector.size(); i++) {
		tmpflag_maxreached = (tmpflag_maxreached && (_data_vector[i] >= _max_vector[i]));
	}

	_max_reached = tmpflag_maxreached;

	_all_combinations.push_back(_data_vector);

	_num_of_steps_passed++;

	return true;
}


std::size_t clsFieldRange::GetPossibleCombinations()
{
	return _num_of_combinations;
}

std::size_t clsFieldRange::GetCurrentPositionInCombinations()
{
	return _num_of_steps_passed;
}


void clsFieldRange::SetRepeats(int repeats)
{
	_repeats = repeats;
}


std::size_t clsFieldRange::GetRepeats()
{
	return _repeats;

}

bool clsFieldRange::IsInitialized()
{
	return _initialized;

}

std::vector<std::vector<FLOAT_TYPE>> clsFieldRange::SampleRecombinations(int num)
{
#ifdef ENABLE_PARALLEL
	lock_guard<mutex> combi_lock(clsFieldRange::_COMBI_GENERATION_MUTEX);
#endif

	if (_all_combinations.size() == 0) {
        std::cout << "empty, generate combis" << std::endl;
		GenerateAllRecombinations();
	}
	
	std::uniform_int_distribution<int> dis(0, static_cast<int>(_all_combinations.size()));

	std::vector<std::vector<FLOAT_TYPE>> tmp_vec;

	for (auto i = 0; i < num; i++) {
		auto rnd_indx = dis(_rnd_engine);
		tmp_vec.push_back(_all_combinations[rnd_indx]);
	}

	return tmp_vec;
}


clsFieldRange& clsFieldRange::operator=(clsFieldRange other)
{
	swap(other);

	return *this;
}


void clsFieldRange::swap(clsFieldRange& other)
{
	_all_combinations.swap(other._all_combinations);
	_combination_iterator = _all_combinations.begin();
	other._combination_iterator = other._all_combinations.begin();

	_data_vector.swap(other._data_vector);
	_min_vector.swap(other._min_vector);
	_max_vector.swap(other._max_vector);

	std::swap(_increment, other._increment);

	// std::swap(_is_data_reset, other._is_data_reset);
	bool tmp_isdatareset = other._is_data_reset;
	other._is_data_reset = _is_data_reset;
	_is_data_reset = tmp_isdatareset;

	// std::swap(_max_reached, other._max_reached);
	bool tmp_maxreached = other._max_reached;
	other._max_reached = _max_reached;
	_max_reached = tmp_maxreached;
	

	std::swap(_num_of_combinations, other._num_of_combinations);
	std::swap(_num_of_steps_passed, other._num_of_steps_passed);

	// sub range - use alternative increments in this range
	_subrange_min_vector.swap(other._subrange_min_vector);
	_subrange_max_vector.swap(other._subrange_max_vector);

	std::swap(_subrange_increment, other._subrange_increment);

	// std::swap(_use_subrange, other._use_subrange);
	bool tmp_usesubrange = other._use_subrange;
	other._use_subrange = _use_subrange;
	_use_subrange = tmp_usesubrange;

	
	std::swap(_combinations_with_subrange, other._combinations_with_subrange);

	std::swap(_repeats, other._repeats);
	std::swap(_initialized, other._initialized);
}


std::vector<std::vector<FLOAT_TYPE>> clsFieldRange::GetRandomCombinations( std::size_t num, bool normalize = false)
{

    std::vector<std::vector<FLOAT_TYPE>> tmp_all_vec;
	
	
	// std::gamma_distribution<double> gamma_dist();

	
    std::vector<FLOAT_TYPE> delta_vector(_min_vector.size());

	// to generalize this generation of numbers we create a delta for each allele
	// then this will be multiplied by a random coefficient within [0,1) in order to achieve a
	// uniform distribution for each allele
	for ( auto current_allele = 0; current_allele < _min_vector.size(); ++current_allele ) {
		delta_vector[current_allele] = _max_vector[current_allele] - _min_vector[current_allele];
	}
    
    //const int FITNESS_LTH_IDX = 0;
    const FLOAT_TYPE FITNESS_LTH_PROB = 0.1;
    //const int FITNESS_DEL_IDX = 1;
    const FLOAT_TYPE FITNESS_DEL_PROB = 0.79;
   // const int FITNESSS_NEU_IDX = 2;
    const FLOAT_TYPE FITNESS_NEU_PROB = 0.1;
    //const FLOAT_TYPE FITNESS_NEU_VAL = 1.0;
    //const int FITNESS_ADV_IDX = 3;
    const FLOAT_TYPE FITNESS_ADV_PROB = 0.01;
    
    // Sanjuan parameters 1.2, 0.092
    std::gamma_distribution<FLOAT_TYPE> dist_gamma(0.5f, 1.0f);
    std::lognormal_distribution<FLOAT_TYPE> dist_lognorm(0.092, 6.0f);
    std::normal_distribution<FLOAT_TYPE> dist_norm(1.0f,0.2f);
    std::discrete_distribution<int> dist_disc( { FITNESS_LTH_PROB, FITNESS_DEL_PROB, FITNESS_NEU_PROB, FITNESS_ADV_PROB});
    //std::uniform_real_distribution<FLOAT_TYPE> dist_uniform( 0.0f, 1.0f );
    
    boost::random::uniform_01<FLOAT_TYPE> dist_uniform;
    
	for ( auto i = 0; i < num; ++i ) {

        std::vector<FLOAT_TYPE> tmp_current_vec(_min_vector.size());

		FLOAT_TYPE tmp_sum = 0.0;

		for (auto current_allele = 0; current_allele < _min_vector.size(); ++current_allele) {
	
            // no sampling is needed
            if ( _min_vector[current_allele] == _max_vector[current_allele] ) {
                tmp_current_vec[current_allele] = _min_vector[current_allele];
                tmp_sum += tmp_current_vec[current_allele];
                continue;
            }
            
#ifdef COMPOSITE_DISTRIBUTION
            
            // only so it won't be uninitialized
            auto tmp_fitness = _min_vector[current_allele];
            
            do {
                // first step - what category
                auto chosen_category = dist_disc(_rnd_engine);
                
                switch (chosen_category) {
                        
                    case FITNESS_ADV_IDX:
                        tmp_fitness = FITNESS_NEU_VAL + dist_gamma(_rnd_engine);
                        break;
                        
                    case FITNESSS_NEU_IDX:
                        tmp_fitness = FITNESS_NEU_VAL;
                        // tmp_fitness = 0.0 - std::fabs(dist_norm(_rnd_engine));
                        break;
                    
                    case FITNESS_DEL_IDX:
                        tmp_fitness = 1.0f - dist_lognorm(_rnd_engine) * dist_norm(_rnd_engine);
                        break;
                        
                    case FITNESS_LTH_IDX:
                        tmp_fitness = 0.0f;
                        break;
                        
                    default:
                        break;
                }
                
                
            } while ( tmp_fitness > _max_vector[current_allele] ||
                     tmp_fitness
                     < _min_vector[current_allele] );
#else
            // uniform distribution
            auto delta_random_factor = dist_uniform(_boost_gen);
            
            auto tmp_fitness = _min_vector[current_allele] + delta_vector[current_allele] * delta_random_factor;
            
            // Quantization
            
            // round to adhere to increment value (binning/quantization)
            //auto tmp_quantize = tmp_current_vec[current_allele] / _increment;
            
            //tmp_quantize = static_cast<FLOAT_TYPE>(std::floor(tmp_quantize));
            
            //tmp_current_vec[current_allele] = _increment * tmp_quantize;
#endif
            
            tmp_current_vec[current_allele] = tmp_fitness;
            tmp_sum += tmp_current_vec[current_allele];
		}


		// normalize if needed
		if (normalize) {
			for ( auto current_allele = 0; current_allele < _min_vector.size(); ++current_allele ) {
				tmp_current_vec[current_allele] = tmp_current_vec[current_allele] / tmp_sum;
			}
		}
		

		// tmp_all_vec.push_back(std::move(tmp_current_vec));
		tmp_all_vec.push_back(tmp_current_vec);
	}
		
	return tmp_all_vec;
}


