// CMulatorMain.cpp : Defines the entry point for the console application.
//

#define DEBUG_VERBOSE

#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>

#include <boost/format.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/bernoulli_distribution.hpp>

#include "CMulator.h"
#include "clsCMulatorABC.h"
#include "CMulatorToolbox.hpp"
#include "PriorSampler.hpp"
#include "SimulationResult.hpp"
#include "ResultsStats.hpp"
#include "ActualDataFile.hpp"



int InferABC( FactorToInfer factor,
                std::string param_filename, std::string actual_data_filename,
                std::string distrib_output_filename, std::string summary_output_filename )
{
    std::cout << "Parameter file: " << param_filename << std::endl;
    std::cout << "Actual data file: " << actual_data_filename << std::endl;
    std::cout << "Posterior distribution output: " << distrib_output_filename << std::endl;
    std::cout << "Summary data output: " << summary_output_filename << std::endl;
 
    
    ActualDataFile actual_data_file;
    try {
        std::cout << "Reading actual data... ";
        actual_data_file.LoadActualData(actual_data_filename);
        std::cout << "Done" << std::endl;
        
        //std::cout << "num of alleles: " << actual_data_file.GetNumberOfAlleles() << std::endl;
        //std::cout << "sd for actual data: " << std::endl;;
        //auto resvec = actual_data_file.GetSDPerAllele();
        
        //for ( auto val : resvec ) {
        //    std::cout << "\t" << val;
       // }
        //std::cout << std::endl;
        //return 0;
    }
    catch (std::exception& e) {
        std::cerr << "Exception whille loading actual data: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown exception whille loading actual data." << std::endl;
        return 1;
    }
    
    
    ZParams my_zparams;
    try {
        std::cout << "Reading parameters... ";
        my_zparams.ReadParameters(param_filename, true);
        std::cout << "Done" << std::endl;
    }
    catch (std::exception& e) {
        std::cerr << "Exception whille loading parameters: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown exception whille loading parameters." << std::endl;
        return 1;
    }
    
    
    clsCMulatorABC abc_object_sim( my_zparams, actual_data_file );
    try {
        std::cout << "Starting ABC." << std::endl;
        
        abc_object_sim.GetUniqueIndexSet(10);
        abc_object_sim.SetImmediateRejection(false);
        abc_object_sim.RunABCInference(factor, 100);
        abc_object_sim.ScaleSimulationResults();
        
        if (my_zparams.GetInt("_coverage", 0) > 0) {
            std::cout << "Starting coverage." << std::endl;
            //abc_object_sim.DoCoverageTest();
            abc_object_sim.DoCoverageCook();
        }
        
        //auto start_coverage_time = std::chrono::high_resolution_clock::now();
        //DoCoverageTest();
        //auto end_coverage_time = std::chrono::high_resolution_clock::now();
        //auto coverage_elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_coverage_time - start_coverage_time);
        //total_running_time = static_cast<double>(coverage_elapsed_ms.count()) / 1000.0;
        //std::cout << "Coverage running time: " << total_running_time << " seconds" << std::endl;
        
        std::cout << "ABC Finished." << std::endl;
    }
    catch (std::exception& e) {
        std::cerr << "exception in runnig abc: " << e.what() << std::endl;
        return 1;
    }
    catch (const char *str) {
        std::cerr << "exception in runnig abc: " << str << std::endl;
        return 1;
    }
    catch (std::string str) {
        std::cerr << "exception in runnig abc: " << str << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "unkown exception while runnig abc." << std::endl;
        return 1;
    }
    
    
    
    std::cout << std::endl << "Writing posterior distribution and summary files... " << std::endl;
    try {
        ResultsStats result_stats;
        
        std::cout << " getting results" << std::endl;
        auto results_vector = abc_object_sim.GetResultsVector(true);
        
        result_stats.SetRejectionThreshold( abc_object_sim.GetRejectionThreshold() );
        
        switch (factor) {
            case Fitness: {
                result_stats.CalculateStatsFitness(results_vector);
                std::string tmp_summary_str = result_stats.GetSummaryFitness();
                
                std::cout << tmp_summary_str << std::endl;
                
                CMulatorToolbox::WriteFitnessDistribToFile(results_vector, distrib_output_filename);
                CMulatorToolbox::WriteStringToFile(summary_output_filename, tmp_summary_str);
                break;
            }
                
            case PopulationSize: {
                result_stats.CalculateStatsPopulationSize(results_vector);
                std::string tmp_summary_str = result_stats.GetSummaryPopSize();
                
                std::cout << tmp_summary_str << std::endl;
                
                CMulatorToolbox::WritePopSizeDistribToFile(results_vector, distrib_output_filename);
                CMulatorToolbox::WriteStringToFile(summary_output_filename, tmp_summary_str);
                break;
            }
                
            case MutationRate: {
                result_stats.CalculateStatsMutation(results_vector);
                std::string tmp_summary_str = result_stats.GetSummaryMutRate();
                
                std::cout << tmp_summary_str << std::endl;
                
                CMulatorToolbox::WriteMutRateDistribToFile(results_vector, distrib_output_filename);
                CMulatorToolbox::WriteStringToFile(summary_output_filename, tmp_summary_str);
                break;
            }
        }
    }
    catch (std::exception& e) {
        std::cerr << "Exception caught while attempting to report stats: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown exception while attempting to report stats." << std::endl;
        return 1;
    }
    
    
    return 0;
}


int RunSingleSimulation(std::string param_filename, std::string output_filename)
{
    // simulator object
    CMulator sim_object;
    ZParams my_zparams;
    
    try {
        my_zparams.ReadParameters(param_filename, true);
        
        if ( my_zparams.GetInt(sim_object.PARAM_ALLELE_NO_INIT_FREQS, 1) > 0 ) {
            throw "Simulating - allele initial frequencies must be provided, and you need to set _no_init_freqs_as_parameters=0";
        }
        
        sim_object.InitMemberVariables(my_zparams);
        
        if ( !sim_object.IsAbleToSimulate() ) {
            throw "Missing parameters for simulation.";
        }
        
        //sim_object.ReadParametersFromFile(param_filename);
        // sim_object.InitializeFromParamFile(param_filename);
    }
    catch (const char* txt) {
        std::cerr << "Exception caught while initializing simulator object:" << std::endl;
        std::cerr << txt << std::endl;
        return 1;
    }
    catch (std::exception& e) {
        std::cerr << "Exception caught while initializing simulator object:" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown exception while initializing simulator object." << std::endl;
        throw;
        return 1;
    }
    
    std::cout << "Simulator object initialized." << std::endl;
    std::cout << "Random seed set to " << sim_object.GetRandomSeed() << std::endl;
    std::cout << "Running simulation..." << std::endl;
    try {
        sim_object.EvolveAllGenerations();
    }
    catch (const char* txt) {
        std::cerr << "Exception: " << txt << std::endl;
    }
    std::cout << "Done." << std::endl;
    
    //std::string sim_output = sim_object.GetAllOutputAsText();
    std::string sim_output = sim_object.GetAllOutputAsTextForR(true);
    
    // also output to screen
    // TODO: make this via parameter
    std::cout << sim_output << std::endl;
    
    std::cout << std::endl;
    
    std::cout << sim_object.GetAllOutputAsMatrix() << std::endl;
    
    std::ofstream outfile(output_filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << output_filename << std::endl;
        return 1;
    }
    std::cout << "Writing to file...";
    outfile << sim_output;
    outfile.close();
    std::cout << " Done." << std::endl;
    
    return 0;
}


void TestMutationRates()
{
    boost::numeric::ublas::matrix<float> min(2, 2);
    boost::numeric::ublas::matrix<float> max(2, 2);
    
    for (auto i = 0; i < min.size1(); ++i) {
        for (auto j = 0; j < min.size2(); ++j) {
            
            min(i, j) = 0;
            max(i, j) = 0.01;
        }
    }    
}


void test_parameters( std::string filename )
{
    ZParams my_zparams(filename, true);
    
    CMulator sim(my_zparams);
    
    // todo: use flag to dump the object as it is created, i don't retain the params along with member variables
    //std::cout << sim.GetAllZParams() << std::endl;
}


void test_range()
{
    
    
    /*
    clsFieldRange range_object( min, max, 0.05, 100);
    
    auto range_output = range_object.GetRandomCombinations(1000000, false);
    
    for (auto current_set : range_output) {
        
        for (auto current_value : current_set) {
            std::cout << current_value << "\t";
        }
        std::cout << std::endl;
    }
     */
    
    
    //PriorSampler<float>::PriorDistributionType dist_type = PriorSampler<float>::PriorDistributionType::UNIFORM;
    //template class PriorSampler<float>;
    std::vector<float> min {0.0};
    std::vector<float> max {2.0};
    
    std::vector<unsigned int> min_int { 40, 40, 40 };
    std::vector<unsigned int> max_int { 100, 100, 100 };
    
    PriorSampler<float> sampler( min, max, PriorDistributionType::FITNESS_COMPOSITE );
    
    auto res_vec = sampler.SamplePrior(10000, min.size());
    
    for ( auto vec : res_vec ) {
        for ( auto val : vec ) {
            
            std::cout << val << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
}


void print_syntaxes(std::string exec_name)
{
    std::cout << "\t" << exec_name << " -fitness <param_file> <actual_data_file> <distribution_output_file> <summary_output_file" << std::endl << std::endl;
    std::cout << "\t" << exec_name << " -mutation <param_file> <output_file>" << std::endl << std::endl;
    std::cout << "\t" << exec_name << " -simulate <param_file> <output_file>" << std::endl << std::endl;
}

void print_welcome()
{
    std::cout << std::endl << "================================================";
    std::cout << std::endl << "                  FITS v0.87 ";
    std::cout << std::endl << "  Flexible Inference from Time-Series data";
    std::cout << std::endl << "================================================";
    std::cout << std::endl;
}


int main(int argc, char* argv[])
{
    
    if (argc <= 1) {
        print_welcome();
        print_syntaxes(argv[0]);
        return 1;
    }
    
    std::string tmp_first_param = argv[1];
    
    
    
    
    if ( tmp_first_param == "-fitness") {
        
        print_welcome();
        
        if (argc != 6) {
            std::cout << "illegal number of arguments (" << argc-1 << "). Syntax is:" << std::endl;
            print_syntaxes(argv[0]);
            return 1;
        }
        
        std::cout << "Inferring fitness" << std::endl;
        
        std::string param_filename = argv[2];
        std::string actual_data_filename = argv[3];
        std::string distrib_output_filename = argv[4];
        std::string summary_output_filename = argv[5];
        
        return InferABC( FactorToInfer::Fitness, param_filename, actual_data_filename,
                            distrib_output_filename, summary_output_filename );
    }
    
    if ( tmp_first_param == "-mutation") {
        
        print_welcome();
        
        if (argc != 6) {
            std::cout << "illegal number of arguments (" << argc-1 << "). Syntax is:" << std::endl;
            print_syntaxes(argv[0]);
            return 1;
        }
        
        std::cout << "Inferring mutation rate" << std::endl;
        
        std::string param_filename = argv[2];
        std::string actual_data_filename = argv[3];
        std::string distrib_output_filename = argv[4];
        std::string summary_output_filename = argv[5];
        
        return InferABC( FactorToInfer::MutationRate, param_filename, actual_data_filename,
                        distrib_output_filename, summary_output_filename );
    }
    
    if ( tmp_first_param == "-popsize") {
        
        print_welcome();
        
        if (argc != 6) {
            std::cout << "illegal number of arguments (" << argc-1 << "). Syntax is:" << std::endl;
            print_syntaxes(argv[0]);
            return 1;
        }
        
        std::cout << "Inferring population size" << std::endl;
        
        std::string param_filename = argv[2];
        std::string actual_data_filename = argv[3];
        std::string distrib_output_filename = argv[4];
        std::string summary_output_filename = argv[5];
        
        return InferABC( FactorToInfer::PopulationSize, param_filename, actual_data_filename,
                        distrib_output_filename, summary_output_filename );
    }
    
    if (tmp_first_param == "-simulate") {
        
        print_welcome();
        
        if (argc != 4) {
            std::cout << "illegal number of arguments (" << argc-1 << "). Syntax is:" << std::endl;
            print_syntaxes(argv[0]);
            return 1;
        }
        
        std::cout << "Running a single simulation" << std::endl;
        
        std::string param_filename = argv[2];
        std::string output_filename = argv[3];
        
        return RunSingleSimulation(param_filename, output_filename);
    }
    
    
    
    if (tmp_first_param == "-test_range") {
        try {
            test_range();
        }
        catch( const char* txt ) {
            std::cerr << "exception in test_range: " << txt << std::endl;
        }
        
        return 1;
    }
 
        
    // we should never reach here...
    print_welcome();
    std::cout << "Invalid suntax (option chosen was " << tmp_first_param << "). Use only the following:" << std::endl;
    print_syntaxes(argv[0]);
    return 1;
}

