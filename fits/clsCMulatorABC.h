#ifndef clsCMulatorABC_hpp
#define clsCMulatorABC_hpp

#include <future>
#include <vector>
#include <chrono>

// performance
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>

#include <cmath> //for pow
#include "ZParams.h"
#include "CMulator.h"
#include "CMulatorToolbox.hpp"
//#include "clsFieldRange.hpp"
#include "SimulationResult.hpp"
#include "PriorSampler.hpp"
#include "ActualDataFile.hpp"

enum FactorToInfer {
    Fitness, MutationRate, PopulationSize
};

class clsCMulatorABC 
{
private:
    
    boost::mt19937 _boost_gen;
    
    const FLOAT_TYPE THRESHOLD_RESET_VALUE = -1.0;
    
    std::size_t _repeats;
    // std::size_t _percent_to_keep;
    std::size_t _sims_to_keep;
    
    FactorToInfer _factor_to_infer;
    
    unsigned int _num_alleles;
    
	FLOAT_TYPE _rejection_threshold;
    bool _use_rejection_threshold; // immediately reject simulations, don't store
    
    ZParams _zparams;
    
    PriorDistributionType _prior_type;

	std::vector<SimulationResult> _simulation_result_vector;
    
    ActualDataFile _actual_data_file;
    
    std::vector<ActualDataEntry> _actual_data_vector;
    std::vector<FLOAT_TYPE> _actual_data_raw_freqs;
    std::vector<int> _selected_actual_generations;
    

    
    // stores samples from the prior
    std::vector< std::vector<FLOAT_TYPE> > _prior_archive;
    
public:
    
    clsCMulatorABC( ZParams sim_params, ActualDataFile actual_data_file );

    std::vector<int> GetUniqueIndexSet( int num_items );

    FLOAT_TYPE ResetRejectionThreshold();
	FLOAT_TYPE SetRejectionThreshold(FLOAT_TYPE new_threshold);
	FLOAT_TYPE GetRejectionThreshold();
    
    void SetImmediateRejection(bool new_val);

    std::size_t GetRepeats();
    
    std::vector<SimulationResult> GetResultsVector(bool only_accepted_results=false);
    
    void RunABCInference( FactorToInfer factor, std::size_t number_of_batches );
    std::vector<SimulationResult> RunFitnessInferenceBatch( std::size_t num_simulations );
    std::vector<SimulationResult> RunPopulationSizeInferenceBatch( std::size_t num_simulations );
    std::vector<SimulationResult> RunMutationInferenceBatch( std::size_t num_simulations );
    
    std::vector<FLOAT_TYPE> GetSDPerAllele( std::size_t start_idx, std::size_t end_idx );
    std::vector<FLOAT_TYPE> GetMADPerAllele( std::size_t start_idx, std::size_t end_idx );
    void DivideEachAllele( std::size_t start_idx, std::size_t end_idx, std::vector<FLOAT_TYPE> value_vector );
    
    FLOAT_TYPE GetDistanceSimActual( MATRIX_TYPE actual_data, MATRIX_TYPE sim_data );
    
    void ScaleSimulationResults();
    void DoCoverageTest();
    
    void DoCoverageCook();
    std::vector<std::size_t> CoverageSingleDataset( std::size_t dataset_idx, std::size_t start_idx, std::size_t end_idx );
    
    // std::vector<SimulationResult> RunGenerationInferenceBatch( std::size_t num_simulations );
};

#endif
