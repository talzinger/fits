//
//  clsFieldRange.hpp
//  NumFieldRange
//
//  Created by Tal Zinger on 13/10/2015.
//  Copyright © 2015 Tal Zinger. All rights reserved.
//

#ifndef clsFieldRange_hpp
#define clsFieldRange_hpp

#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/gamma.hpp>

#include "CMulator.h" // just for the FLOAT_TYPE , MATRIX_TYPE type definitions

// using namespace std;

#define FIRST_DIGIT 0

class clsFieldRange {
public:

	enum PriorDistributionType {
		UNIFORM,
		GAMMA
	};

	enum RangeType {
		VECTOR_RANGE,
		MATRIX_RANGE
	};

private:

	std::vector<std::vector<FLOAT_TYPE>> _all_combinations;
	std::vector<std::vector<FLOAT_TYPE>>::iterator _combination_iterator;

	std::vector<FLOAT_TYPE> _data_vector;
	std::vector<FLOAT_TYPE> _min_vector;
	std::vector<FLOAT_TYPE> _max_vector;

	MATRIX_TYPE _min_matrix;
	MATRIX_TYPE _max_matrix;
	
	FLOAT_TYPE _increment;
	
	std::mt19937_64 _rnd_engine;
    boost::mt19937 _boost_gen;
    
    // boost::random::mt19937 _rnd_engine;
	long long _rnd_seed;
    int _rnd_int_seed;

	bool _is_data_reset;
	bool _max_reached;
	
	int _num_of_combinations;
	int _num_of_steps_passed;
	
	// sub range - use alternative increments in this range
	std::vector<FLOAT_TYPE> _subrange_min_vector;
	std::vector<FLOAT_TYPE> _subrange_max_vector;
	FLOAT_TYPE _subrange_increment;
	bool _use_subrange;
	int _combinations_with_subrange;
	
    std::size_t _repeats;
	
	bool _initialized;

	RangeType _range_type;
	

	// to be used to generate the combination vector
	// exact copy of the old Next. Kept so it would be easy to go back to old version
	bool GenerateNextCombination();

	bool GenerateAllRecombinations();

	void InitializeMemberVariables();

public:
	PriorDistributionType _prior_dist_type;
	FLOAT_TYPE _prior_dist_gamma_alpha;
	FLOAT_TYPE _prior_dist_gamma_beta;
	
    bool IsInitialized();
	
	void SetRepeats( int repeats );
	std::size_t GetRepeats();
	
	std::vector<std::vector<FLOAT_TYPE>> GetAllCombinations();
	std::vector<std::vector<FLOAT_TYPE>> SampleRecombinations(int num);

    std::vector<std::vector<FLOAT_TYPE>> GetRandomCombinations(std::size_t num, bool normalize); // normalize - such that sum of values in combination is 1

    std::vector<MATRIX_TYPE> GetRandomCombinationMatrices( std::size_t num );

	
	clsFieldRange();
    
	clsFieldRange( std::vector<FLOAT_TYPE> min_values,
                  std::vector<FLOAT_TYPE> max_values,
                  FLOAT_TYPE increment,
                  int repeats = 1,
                  PriorDistributionType prior_dist = PriorDistributionType::UNIFORM );
    
	clsFieldRange( MATRIX_TYPE min_values,
                  MATRIX_TYPE max_values,
                  FLOAT_TYPE increment,
                  int repeats = 1,
                  PriorDistributionType prior_dist = PriorDistributionType::UNIFORM);

	// clsFieldRange( std::vector<FLOAT_TYPE> min_values, std::vector<FLOAT_TYPE> max_values, FLOAT_TYPE increment );
	clsFieldRange( const clsFieldRange &original );
	
	std::size_t GetPossibleCombinations();
	std::size_t GetCurrentPositionInCombinations();
	

	void swap(clsFieldRange& other);
	clsFieldRange& operator=(clsFieldRange other);
};


#endif /* clsFieldRange_hpp */
