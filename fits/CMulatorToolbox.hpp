//
//  CMulatorToolbox.hpp
//  CMulatorLib
//
//  Created by Tal Zinger on 03/01/2016.
//  Copyright © 2016 Tal Zinger. All rights reserved.
//

#ifndef CMulatorToolbox_hpp
#define CMulatorToolbox_hpp

#include <iostream>
// #include <cmath>
#include <fstream>
#include <iomanip> // string set percision
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>


#include <boost/numeric/odeint/util/ublas_wrapper.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "ZParams.h"
#include "CMulator.h"
#include "SimulationResult.hpp"


/* Helper functions and data types for CMulator */
class CMulatorToolbox {

public: /* Constants */

	
	// perecent of the results within bounds of definition to be considered as such
	/*static const int CATEGORY_INCLUSION_RELAXED_THRESHOLD = 50;
    static const int CATEGORY_INCLUSION_STRICT_THRESHOLD = 95;*/

	const FLOAT_TYPE _DEFAULT_DISTANCE_SIM_ACTUAL = -1.0;

	// define fitness thresholds [min,max)
	/*static const FLOAT_TYPE FITNESS_LETHAL_MIN;
	static const FLOAT_TYPE FITNESS_LETHAL_MAX;
	static const FLOAT_TYPE FITNESS_DELETERIOUS_MIN;
	static const FLOAT_TYPE FITNESS_DELETERIOUS_MAX;
	static const FLOAT_TYPE FITNESS_NEUTRAL_MIN;
	static const FLOAT_TYPE FITNESS_NEUTRAL_MAX;
	static const FLOAT_TYPE FITNESS_ADVANTAGEOUS_MIN;
	static const FLOAT_TYPE FITNESS_ADVANTAGEOUS_MAX;*/
    

	// fields for sim data output if a new method is to be written
	static const int SIM_DATA_COLUMN_UID = 0;
	static const int SIM_DATA_COLUMN_GENERATION = 1;
	static const int SIM_DATA_COLUMN_ALLELE = 2;
	static const int SIM_DATA_COLUMN_FREQ = 3;
	static const int SIM_DATA_COLUMN_INIT_FREQ = 4;
	static const int SIM_DATA_COLUMN_ALLELE_FITNESS = 5;
	static const int SIM_DATA_COLUMN_ALLELE_POS = 6;
	static const int SIM_DATA_COLUMNS = 7;

    static const int REPORT_FLOAT_PRECISION = 3;

	/* Data Types */
public:

	// TODO: move code to cpp
    //static std::vector<FLOAT_TYPE> GetHeuristicMaxFitness( std::vector<ActualDataEntry> actual_data );
    //static std::vector<FLOAT_TYPE> GetHeuristicMinFitness( std::vector<ActualDataEntry> actual_data );
	
	
	
	/* Functions */
	static FLOAT_TYPE GetDistanceSimActual( std::vector<FLOAT_TYPE> actual_data, std::vector<FLOAT_TYPE> sim_data );
    
    static FLOAT_TYPE GetDistanceSimActual( MATRIX_TYPE actual_data, MATRIX_TYPE sim_data );
    
	//static FLOAT_TYPE GetDistanceSimActual(ActualDataEntry actual_data, SimulationResult sim_data);
    static std::vector<FLOAT_TYPE> GetDistanceVectorSimActual( std::vector<FLOAT_TYPE> actual_data, std::vector<FLOAT_TYPE> sim_data, std::size_t num_alleles  );
    
    static FLOAT_TYPE CalculateDistKernel( std::vector<FLOAT_TYPE> actual_data, std::vector<FLOAT_TYPE> sim_data );

	//static FLOAT_TYPE OrderResultsAndUpdateThreshold(std::vector<SimulationResult>& result_vector, int max_idx_to_keep, bool erase_over_threshold);
	
	static void WriteFitnessDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename);
    static void WriteMutRateDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename);
    static void WritePopSizeDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename);
    
    static void WriteStringToFile( std::string filename, std::string str );
    
	//static std::string ReportStatsFitness(const std::vector<SimulationResult>& result_vector, std::string filename);
    

	static std::vector<FLOAT_TYPE> GetMinFitness(const std::vector<SimulationResult>& result_vector);
	static std::vector<FLOAT_TYPE> GetMaxFitness(const std::vector<SimulationResult>& result_vector);

    
    //static ZParams ReadParametersFromFile( std::string param_filename );
    
    static FLOAT_TYPE FreqsSameMagnitude( FLOAT_TYPE original_val );
    
    // building the posteriors from simulation results
    static std::vector<std::vector<FLOAT_TYPE>> BuildPosteriorFitness(const std::vector<SimulationResult>& result_vector);
    
    // static std::vector<float> GetPosteriorPval
};

#endif /* CMulatorToolbox_hpp */
