//
//  CMulatorToolbox_Report.cpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "CMulatorToolbox.hpp"



void CMulatorToolbox::WriteFitnessDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename)
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    // header
    outfile << "sim_id" << "\t" << "distance";
    for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
        outfile << "\t" << "allele" << i;
    }
    outfile << std::endl;
    
    
    // 2016-09-04 adding sanity check to report problem
    if ( result_vector.empty() ) {
        std::cerr << "Error in WriteToFile: result_vector is empty. Filename: " << filename << std::endl;
        std::string my_err = "";
        my_err = "Error in WriteToFile: result_vector is empty. Filename: ";
        my_err += filename;
        throw my_err.c_str();
    }
    
    for ( auto tmp_entry : result_vector ) {
        
        outfile << tmp_entry.sim_id << "\t" << tmp_entry.distance_from_actual;
        
        for ( auto tmpval : tmp_entry.fitness_values ) {
            outfile << "\t" << tmpval;
        }
        
        outfile << std::endl;
    }
    
    outfile.close();
}

void CMulatorToolbox::WriteMutRateDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename)
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    // header
    outfile << "sim_id" << "\t" << "distance";
    
    for (auto i = 0; i < result_vector[0].fitness_values.size(); ++i) {
        for (auto j = 0; j < result_vector[0].fitness_values.size(); ++j) {
            outfile << "\t" << "allele" << i << "_" << j;
        }
    }
    outfile << std::endl;
    
    
    // 2016-09-04 adding sanity check to report problem
    if ( result_vector.empty() ) {
        std::cerr << "Error in WriteToFile: result_vector is empty. Filename: " << filename << std::endl;
        std::string my_err = "";
        my_err = "Error in WriteToFile: result_vector is empty. Filename: ";
        my_err += filename;
        throw my_err.c_str();
    }
    
    for ( auto tmp_entry : result_vector ) {
        
        outfile << tmp_entry.sim_id << "\t" << tmp_entry.distance_from_actual;
        
        for (auto row = 0; row < result_vector[0].mutation_rates.size1(); ++row) {
            for (auto col = 0; col < result_vector[0].mutation_rates.size2(); ++col) {
                outfile << "\t" << boost::format("%.3d") % tmp_entry.mutation_rates(row,col);
            }
        }
        outfile << std::endl;
    }
    
    outfile.close();
}

void CMulatorToolbox::WritePopSizeDistribToFile(const std::vector<SimulationResult>& result_vector, std::string filename)
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    // header
    outfile << "sim_id" << "\t" << "distance" << "\t" << "N" << std::endl;
    
    
    // 2016-09-04 adding sanity check to report problem
    if ( result_vector.empty() ) {
        std::cerr << "Error in WriteToFile: result_vector is empty. Filename: " << filename << std::endl;
        std::string my_err = "";
        my_err = "Error in WriteToFile: result_vector is empty. Filename: ";
        my_err += filename;
        throw my_err.c_str();
    }
    
    for ( auto tmp_entry : result_vector ) {
        outfile << tmp_entry.sim_id << "\t"
                << tmp_entry.distance_from_actual << "\t"
                << boost::format("%e") % tmp_entry.N << std::endl;
    }
    
    outfile.close();
}


void CMulatorToolbox::WriteStringToFile( std::string filename, std::string str )
{
    std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
    
    if (!outfile.is_open()) {
        std::cerr << "unable to open file for writing: " << filename << std::endl;
        throw "unable to open file for writing: " + filename;
    }
    
    outfile << str;
    
    outfile.close();
}
