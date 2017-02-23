//
//  CMulator.h
//  simFourAllelesV1
//
//  Created by Tal Zinger on 6/9/15.
//  Copyright (c) 2015 Tal Zinger. All rights reserved.
//


//#define _SCL_SECURE_NO_WARNINGS

//#ifndef __simFourAllelesV1__CMulator__
//#define __simFourAllelesV1__CMulator__
#ifndef __CMulator__
#define __CMulator__

#include <cmath>
#include <chrono>
#include <random>
#include <iostream>
#include <string>

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "ZParams.h"

#define FLOAT_TYPE float
#define MATRIX_TYPE boost::numeric::ublas::matrix<FLOAT_TYPE>


class CMulator {

private:
    
    bool _debug_mode;
    
    // No need to create a member of this object.
    // create this object only when initializing and then use only member variables, not parameters
    // ZParams _zparams;
	
	/* Basic Parameters For Simulation */
	// TODO: turn all int to unsigend int
	int _N;		// current population size
	int _N0; 	// founding population size
    
    int _alt_N; // alternative population size to introduce noise to simulation
    int _old_N; // the original value, to revert from altenative value
    int _alt_generation; // if given, in this generation the popsize will alter; otherwise it would be random.
    
	int _num_generations;
	int _num_alleles;
	int _bottleneck_interval;
	int _bottleneck_size;
	int _sample_size;
	int _generation_shift;
    
	std::string _uid;
    
	bool _no_init_freqs_as_parameters; // init freqs will be taken later
	bool _logistic_growth;
	int _logistic_growth_K;		// population capacity
	int _logistic_growth_t;		// time; couples to number of generations but may be reset following bottleneck
	FLOAT_TYPE _logistic_growth_r;	// growth rate
	int _repeats;

	int _wt_allele_index;
	bool _infer_wt_automatically;		// don't expect min and max values for each alleles. infer the wt, assign 1 to min and max, and assign arbitrary min and max for rest of alleles

    bool _infer_mutation_rate;
    
	std::string _R_filename;
	int _max_generations_till_fixation;
	unsigned int _top_percent_to_keep; // mutually exclusive with _top_sims_to_keep
	unsigned int _top_sims_to_keep; // mutually exclusive with _top_percent_to_keep
	bool _skip_stochastic_step;		// use only deterministic selection
	
    /* Fitness */
	//FLOAT_TYPE _fitness_increment;
	std::vector<FLOAT_TYPE> _allele_fitness;
	std::vector<FLOAT_TYPE> _allele_init_freqs;
	std::vector<FLOAT_TYPE> _allele_max_fitness;
	std::vector<FLOAT_TYPE> _allele_min_fitness;
	
    // 2016-11-15
    // What data is available to the simulator - the reset would have to be supplied before actually running simulation
    bool _available_mutrate;
    bool _available_mutrate_range;
    bool _available_fitness;
    bool _available_fitness_range;
    bool _available_popsize;
    bool _available_popsize_range;
    bool _available_essential;
    bool _available_inference;
    bool _available_initfreqs;
	
	// store a single simulation - column for allele, row for generation
	// TODO: maybe switch to matrix data type to enable matrix multiplications?
	typedef boost::multi_array<FLOAT_TYPE,2> allele_dataArray;
	typedef allele_dataArray::index _allele_index;
	allele_dataArray _sim_data;
    // slowly develop this matrix side-by-side with the above _sim_data
    // until fully operational and then replace
    boost::numeric::ublas::matrix<FLOAT_TYPE> _sim_data_matrix;
    //boost::numeric::ublas::vector<float> _current_generation_vector;
	
	/* Mutation rates */
	typedef boost::multi_array<FLOAT_TYPE,2> allele_mutationArray;
	typedef allele_mutationArray::index _mutation_rate_index;
	allele_mutationArray _mutation_rates;
    
    MATRIX_TYPE _mutation_rates_matrix;
    MATRIX_TYPE _min_mutation_rate_matrix;
    MATRIX_TYPE _max_mutation_rate_matrix;
	/* END Basic Parameters For Simulation */
	
    unsigned int _Nlog_min;
    unsigned int _Nlog_max;
	
	/* Runtime Parameters */
	int _current_generation;
    
    // TODO: why keeping these as member variables?
	std::vector<FLOAT_TYPE> _allele_fitness_bar;
	std::vector<FLOAT_TYPE> _allele_fitness_adjusted;
	
	FLOAT_TYPE _fixation_upper_threshold;
	FLOAT_TYPE _fixation_lower_threshold;
	FLOAT_TYPE _fixation_consecutive_generations; // for at least this num of generation threthold is exceeded to count as fixation
	/* END Runtime Parameters */
	
	
	
	/* Assisting Functions */
	
	
    void InitBasicVariables( ZParams zparams );
    void InitInitialAlleleFreqs( ZParams zparams );
    
    void InitGeneralInferenceVariables( ZParams zparams );
    
    void InitPopulationSizeInference( ZParams zparams );
    
    // boolean flag to know whether we should be expecting the mutation rate at all (or the range was given)
    void InitMutationRates( ZParams zparams, bool available_mutrate_range );
    
    void InitMutationInferenceVariables( ZParams zparams );
    
    void InitFitnessValues( ZParams zparams );
    void InitFitnessInference( ZParams zparams );
    
    void InitFixationVariables( ZParams zparams );
	/* END Assisting Functions */
	
	
	/* Technical Parameters */
	// for stochastic step
	//std::mt19937_64 _my_random_generator;
	unsigned long long _my_random_seed;
    boost::mt19937 _boost_gen;

	unsigned int _time_for_seeding;
	
	int _first_sim_in_batch;
	
	// if the object was created with no parameter file, make sure it is later initialized with parameters
	bool _initialized_with_parameters;
	
	// to mitigate floating point rounding errors
	FLOAT_TYPE _epsilon_float_compare;
	/* END Technical Parameters */
	
	
	
	/* Internal Sanity Checks */
	bool IsValid_Allele( int allele ) const { return ( allele >= 0 && allele < _num_alleles ); }
    bool IsValid_NumAlleles( int num_alleles ) const { return num_alleles >= 2; }
	
    bool IsValid_Generation( int generation ) const
    {
        if ( _num_generations>0 ) return ( generation >=0 && generation <= _num_generations);
        else return ( generation >=0  );
    }
    
	bool IsValid_Frequency( FLOAT_TYPE freq ) const { return freq >= 0 && freq <= 1.0; }
	bool IsValid_Fitness( FLOAT_TYPE fitness ) const { return fitness >= 0; }
	bool IsValid_PopulationSize( int N ) const { return N>0; }
	bool IsValid_SamplingSize( int sample_size ) const { return sample_size > 0 && sample_size <=_N; }
	bool IsValid_BottleneckInterval( int interval ) const { return interval < _num_generations; } // non-positive treated as no bottleneck
	bool IsValid_BottleneckSize( int size ) const { return size > 0 && size <= _N; }
	/* End Internal Sanity Checks */
	
// set functions
// used to be public
    // but many of them require reset to make the simulated data in sync with the parameters
    // I don't like this and don't like the static param, library
    
public:
    // using the available data, is the simulator able to do that
    bool IsAbleToInferFitness() const {
        if (!_available_mutrate) std::cerr << "mutrate not available" << std::endl;
        if (!_available_popsize) std::cerr << "popsize not available" << std::endl;
        if (!_available_essential) std::cerr << "essential not available" << std::endl;
        if (!_available_fitness_range) std::cerr << "fitness range not available" << std::endl;
        return _available_mutrate && _available_popsize && _available_essential && _available_fitness_range;
    }
    bool IsAbleToInferPopulationSize() const {
        if (!_available_mutrate) std::cerr << "mutrate not available" << std::endl;
        if (!_available_fitness) std::cerr << "fitness not available" << std::endl;
        if (!_available_essential) std::cerr << "essential not available" << std::endl;
        if (!_available_popsize_range) std::cerr << "popsize range not available" << std::endl;
        return _available_mutrate && _available_fitness && _available_essential && _available_popsize_range;
    }
    bool IsAbleToInferMutationRate() const {
        if (!_available_mutrate_range) std::cerr << "mutrate range not available" << std::endl;
        if (!_available_popsize) std::cerr << "popsize not available" << std::endl;
        if (!_available_essential) std::cerr << "essential not available" << std::endl;
        if (!_available_fitness) std::cerr << "fitness not available" << std::endl;
        return _available_fitness && _available_popsize && _available_essential && _available_mutrate_range;
    }
    bool IsAbleToSimulate() const {
        if (!_available_mutrate) std::cerr << "mutrate rates not available" << std::endl;
        if (!_available_popsize) std::cerr << "popsize not available" << std::endl;
        if (!_available_essential) std::cerr << "essential not available" << std::endl;
        if (!_available_fitness) std::cerr << "fitness not available" << std::endl;
        return _available_mutrate && _available_popsize && _available_essential && _available_fitness;
    }
private:
    // Past Set functions
    void SetNumOfGeneration( int generations );
    void SetAlleleNumber( int alleles );
    void SetBottleneckInterval( int interval );
    void SetBottleneckSize( int size );
    
    void SetSamplingSize( int sampling_size );
    void SetSkipStochasticStep(bool skip);
    
    
    void SetRepeats( int repeats );
    void SetTopPercentToKeep(unsigned int new_percent);
    void SetTopSimsToKeep(unsigned int new_sims);
    
    void EnableLogisticGrowth();
    void EnableLogisticGrowth( int k, FLOAT_TYPE r );
    void DisableLogisticGrowth();
    bool GetLogisticGrowth( int &k, FLOAT_TYPE &r ) const;
    
    void SetAlleleMinFitness( int allele, FLOAT_TYPE fitness );
    void SetAlleleMaxFitness( int allele, FLOAT_TYPE fitness );
    void SetAlleleFitnessValue( int allele, FLOAT_TYPE fitness );
    
    
    void SetMutationRate( int from, int to, FLOAT_TYPE rate );
    
public:
	/* Constructors */
	CMulator();
	//CMulator( std::string param_filename );
    
	CMulator( const CMulator &original );
    CMulator( const ZParams &zparams );
    
    //TODO: move ctor
	/* END Constructors */
	
	/* Get/Set Functions */
    void SetGenerationShift( int shift );
    void SetAlleleInitFreq( int allele, FLOAT_TYPE freq );
    void SetWTAllele( int wt_index, FLOAT_TYPE min_fitness_nonwt, FLOAT_TYPE max_fitness_nonwt );
    void SetFitnessValues( std::vector<FLOAT_TYPE> _allele_fitness );
    
    void InitMemberVariables( ZParams zparams );
    
	// to be used only when the object is still uninitialized (default ctor used)
	//void InitializeFromParamFile(std::string param_filename);
    //void ReadParametersFromFile( std::string filename );
    
	void SetSimUID( std::string new_sim_uid );
	int GetPopulationSize() const;
	void SetPopulationSize( int N );
    
	int GetNumOfGenerations() const;
	
	int GetAlleleNumber() const;
	
	int GetBottleneckInterval() const;
		
	int GetBottleneckSize() const;
	
	int GetSamplingSize() const;
	
	int GetGenerationShift() const;
	
	std::string GetSimUID() const;

	bool GetSkipStochasticStep() const;
	
	//FLOAT_TYPE GetFitnessIncrement() const;
    
    // to avoid direct access to array - allow to raise exceptions
    FLOAT_TYPE GetSimulatedFrequency( int generation, int allele ) const;
    void SetSimulatedFrequency(  int generation, int allele, FLOAT_TYPE frequency );

	unsigned int GetTopPercentToKeep() const;
	
	unsigned int GetTopSimsToKeep() const;
	
	// this shouldn't be here.. but rather in a wrapping class
	// but I can't afford any other thread calling parameters class
	int GetRepeats() const;
	
	FLOAT_TYPE GetInitAlleleFreq(int allele) const;
	
	std::vector<FLOAT_TYPE> GetAlleleInitFreqs() const;
	void SetAlleleInitFreqs( std::vector<FLOAT_TYPE> freqs );

	std::string GetRFileName() const;
	void SetRFileName( const std::string& filename );
	
    std::vector<FLOAT_TYPE> GetAlleleFitnessValues() const;
	std::vector<FLOAT_TYPE> GetAlleleMaxFitnessValues() const;
	std::vector<FLOAT_TYPE> GetAlleleMinFitnessValues() const;

	int GetWTAllele() const;
	
	FLOAT_TYPE GetMutationRate( int from, int to ) const;
    
    void SetMutationRateMatrix( MATRIX_TYPE new_mutation_matrix );
    MATRIX_TYPE GetMutationRateMatrix() const;
    MATRIX_TYPE GetMinMutationRateMatrix() const;
    MATRIX_TYPE GetMaxMutationRateMatrix() const;
    
	FLOAT_TYPE GetAlleleFreq( int generation, int allele ) const;
    

	/* END Get/Set Functions */
	
	
	
	/* Other Information Retreival */
	// int IsFirstInBatch();
	// functions for debugging, should be ommitted in future.
	unsigned int GetRandomSeed() const;
	void SetRandomSeed(unsigned int new_seed);
	/* END Other Information  Retrieval */
	
	
	
	/* Operational Functions */
	int EvolveToGeneration( int target_generation );
	int EvolveAllGenerations();
	
	/* This does not retain actual data, justt last generation.
	   We assume this could run for many generation, don't want to
	   overload the system.
	   Evolve up till a given hard limit or until fixation
	 Return the generation where reached fixation.
	 Return negative value if hard limit reached
	 */
	int EvolveUntilFixation( int tested_allele );
	
	//void UpdateParameters(); // make sure Parameters hold the values stored in the member variables
	/* END Operational Functions */
	
	
	
	/* Output Functions */
	//void PrintDataForR( const std::string& filename, bool append_file, bool zero_based_alleles, bool print_header, bool extended_header, int genshift, int resample );

	//void PrintDataAsMatrix( std::string seperator ) const;
	
	//void PrintSingleAllele( std::string filename, int allele, bool use_newline );
	
	//enum ReportType { FREQS_ONLY, FITNESS_ONLY, FULL };
	
	//void PrintReport( ReportType report_type );
	
    std::string GetAllOutputAsTextForR( bool header = true ) const;
	std::string GetAllOutputAsText( bool header = true ) const;
	std::string GetAllOutputAsCSV() const;
    
    MATRIX_TYPE GetAllOutputAsMatrix() const;
	
    // MATRIX_TYPE SamplePopulationForSequencing();
    
	// std::string GetAllOutputForR( bool zero_based_alleles, bool print_header, bool extended_header, int genshift, int resample );

	// void GetSimulationData( std::string job_uid, CMulatorData &sim_data );

    // resample - randomally select indifviduals according to _sample_size
    // can't be const because the random engine changes state
    std::vector<FLOAT_TYPE> GetRawFrequencyData( bool resample = false );
	
	// export only subset of generation, if actual data is not continuous (i.e missing generations)
    // resample - randomally select indifviduals according to _sample_size
    std::vector<FLOAT_TYPE> GetRawFrequencyData( std::vector<int> selected_generations, bool resample = false );
	/* END Output Functions */
	
	/* Reset */
	// delete simulated data, set current generation to 0.
	void Reset_Soft();
	void Reset_Soft( std::string new_sim_uid );
	
	// initialize member variables from loaded parameters
    // becomes irrelevant if parameters are only refered to in ctor
    // and not constantly available as static
	// void Reset_Hard();
	/* END Reset */	
	
public:	
	// not supported any longer
    // const std::string PARAM_ACTUAL_DATA_FILENAME = "_actual_data_file";
    // const std::string PARAM_FIRST_SIM_IN_BATCH = "_first_sim_in_batch";
    // const std::string PARAM_ALLELE_FITNESS_INCRENEMT = "_fitness_increment";
    // const FLOAT_TYPE PARAM_ALLELE_FITNESS_INCRENEMT_DEFAULT = 0.05;
    // const std::string PARAM_R_FILENAME = "_out_R_file";
    
    
    /* Constants */
	const int MAX_GENERATION = 1000;
	
	const std::string PARAM_POPULATION_SIZE = "_N";
	const int PARAM_POPULATION_SIZE_DEFAULT = 0;
	
    const std::string PARAM_ALT_POPULATION_SIZE = "_alt_N";
    const std::string PARAM_ALT_GENERATION = "_alt_generation";

	const std::string PARAM_GENERATION_SHIFT = "_generation_shift";
	const int PARAM_GENERATION_SHIFT_DEFAULT = 0;
	
	
	const std::string PARAM_BOTTLENECK_SIZE = "_bottleneck_size";
	const std::string PARAM_BOTTLENECK_INTERVAL = "_bottleneck_interval";
	
	const std::string PARAM_SAMPLE_SIZE = "_sample_size";
    
	const std::string PARAM_NUM_GENERATIONS = "_num_generations";
    
	const std::string PARAM_NUM_ALLELES = "_num_alleles";
    
	const std::string PARAM_ALLELE_INIT_FREQ = "_allele_init_freq";
	const std::string PARAM_ALLELE_NO_INIT_FREQS = "_no_init_freqs_as_parameters"; // assume it will be provided later
    
    const std::string PARAM_INFER_MUTATION_RATE_FLAG = "_infer_mutation_rate";
	const std::string PARAM_MUTATION_RATE = "_mutation_rate";
    const std::string PARAM_MIN_LOG_MUTATION_RATE = "_min_log_mutation_rate";
	const std::string PARAM_MAX_LOG_MUTATION_RATE = "_max_log_mutation_rate";
    
	const std::string PARAM_ALLELE_FITNESS = "_allele_fitness";
	const std::string PARAM_ALLELE_MAX_FITNESS = "_allele_max_fitness";
	const std::string PARAM_ALLELE_MIN_FITNESS = "_allele_min_fitness";
	
	const std::string PARAM_SIM_ID = "_sim_uid";
	const std::string PARAM_SIM_REPEATS = "_num_repeats";
    
    const std::string PARAM_MANUAL_SEED = "_manual_seed";
	
	const std::string PARAM_LOGISTIC_GROWTH = "_logistic_growth";
	const std::string PARAM_LOGISTIC_GROWTH_K = "_logistic_growth_K";
	const std::string PARAM_LOGISTIC_GROWTH_r = "_logistic_growth_r";
	const FLOAT_TYPE ALLELE_FITNESS_DEFAULT_MIN = 0.0;
	const FLOAT_TYPE ALLELE_FITNESS_DEFAULT_MAX = 2.0;

	const std::string PARAM_FIXED_UPPER_THRESHOLD = "_fixation_upper_threshold";
	const std::string PARAM_FIXED_LOWER_THRESHOLD = "_fixation_lower_threshold";
	const std::string PARAM_FIXED_CONSECUTIVE_GENERATIONS = "_fixation_consecutive_generations";
	
	const std::string PARAM_EPSILON_FLOAT_COMPARE = "_epsilon_float_compare";
	const FLOAT_TYPE PARAM_EPSILON_FLOAT_COMPARE_DEFAULT = 0.00001f; // within this range float values are considered equal
	const FLOAT_TYPE PARAM_EPSILON_SIM_DEFAULT = 0.001f; // default range within actual data to accept simulations
	
	const int PARAM_DEFAULT_VAL_INT = -1;
	const FLOAT_TYPE PARAM_DEFAULT_VAL_FLOAT = -1.0;
	const std::string PARAM_DEFAULT_VAL_STRING = "NA";
	
	//const std::string PARAM_TOP_PERCENT_TO_KEEP = "_top_percent_to_keep";
	//const unsigned int TOP_PERCENT_TO_KEEP_DEFAULT = 0;
	
	//const std::string PARAM_TOP_SIMS_TO_KEEP = "_top_sims_to_keep";
	//const unsigned int TOP_SIMS_TO_KEEP_DEFAULT = 0;
    
    const std::string PARAM_ACCEPTANCE_RATE = "_acceptance_rate";
    const FLOAT_TYPE ACCEPTANCE_RATE_DEFAULT = 0.01f;
    

	const std::string PARAM_SKIP_STOCHASTIC_STEP = "_skip_stochastic_step";
	const int PARAM_SKIP_STOCHASTIC_STEP_DEFAULT = 0;

    // composite prior for fitness inference
    // if one is manually chosen - all rest must be also
    const std::string PARAM_PRIOR_FRACTION_DEL = "_prior_fraction_del";
    const std::string PARAM_PRIOR_FRACTION_NEU = "_prior_fraction_neu";
    const std::string PARAM_PRIOR_FRACTION_ADV = "_prior_fraction_adv";
    
	/* END Constants */
};

#endif /* defined(__CMulator__) */
