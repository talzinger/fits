//
//  CMulator_IO.cpp
//  CMulatorLib
//
//  Created by Tal Zinger on 31/12/2015.
//  Copyright © 2015 Tal Zinger. All rights reserved.
//

// #define DEBUG_VERBOSE

//#include <stdio.h>
#include <iostream>
#include <fstream>

#include "CMulator.h"
//#include "numR.h"


MATRIX_TYPE CMulator::GetAllOutputAsMatrix() const
{
    MATRIX_TYPE tmp_matrix(_current_generation, _num_alleles);
    
    for ( auto gen=0; gen<_current_generation; ++gen ) {
        
        for ( auto allele=0; allele<_num_alleles; ++allele ) {
            
            tmp_matrix(gen,allele) = _sim_data[gen][allele];
        }
    }
    
    return tmp_matrix;
}


std::string CMulator::GetAllOutputAsCSV() const
{
	if (!_initialized_with_parameters) {
		throw " GetAllOutputAsCSV: object not initialized with parameters.";
	}
	
    std::string tmp_str = "";
	if (_current_generation==0) {
		return tmp_str;
	}
	
	// header 
	for (int myAllele=0; myAllele<_num_alleles; myAllele++) {
		tmp_str +=  "Allele " + std::to_string(myAllele) + ",";
	}
	
	tmp_str += "\n";
	
	
	if (_current_generation>1) {
		
		for (auto myGeneration=0; myGeneration<_current_generation; myGeneration++) {
			
			for (auto myAllele=0; myAllele<_num_alleles; myAllele++) {
				
				tmp_str += std::to_string(_sim_data[myGeneration][myAllele]) + ",";
			}
			
			tmp_str += "\n";
		}
	}
	
	return tmp_str;	
}


/*
void CMulator::ReadParametersFromFile( std::string param_filename )
{
    ZParams my_zparams;
    
    try {
        // parameters should be read only
        std::cout << "Reading parameters file: " << param_filename << std::endl;
        my_zparams.ReadParameters(param_filename, true);
    }
    catch( const char* txt ) {
        std::cerr << "Error in zparams: " << txt << std::endl;
    }
    catch (...) {
        std::cerr << "Error while reading parameter file: " << param_filename << std::endl;
    }
	
    if (!my_zparams.IsEmpty()) {
        std::cout << "Zparams:" << std::endl;
        std::cout << my_zparams.GetAllParameters() << std::endl;
    }
    
    _initialized_with_parameters = true;
    
    InitMemberVariables(my_zparams);
}
*/


std::string CMulator::GetAllOutputAsText(bool header) const
{
	if (!_initialized_with_parameters) {
		throw " GetAllOutputAsText: object not initialized with parameters.";
	}
	
    std::string tmp_str = "";
	if (_current_generation==0) {
		return tmp_str;
	}
	
	// header
    if (header) {
        for (auto myAllele=0; myAllele<_num_alleles; myAllele++) {
            tmp_str +=  "allele" + std::to_string(myAllele) + "\t";
        }
        tmp_str += "\n";
    }
	
	if (_current_generation>1) {
		
		for (auto myGeneration=0; myGeneration<_current_generation; myGeneration++) {
			
			for (auto myAllele=0; myAllele<_num_alleles; myAllele++) {
				
				tmp_str += std::to_string(_sim_data[myGeneration][myAllele]) + "\t";
			}
			
			tmp_str += "\n";
		}
	}
	
	return tmp_str;
}

std::string CMulator::GetAllOutputAsTextForR( bool header ) const
{
    if (!_initialized_with_parameters) {
        throw " GetAllOutputAsText: object not initialized with parameters.";
    }
    
    std::string tmp_str = "";
    if (_current_generation==0) {
        return tmp_str;
    }
    
    // header
    if (header) {
        tmp_str += "gen\tbase\tfreq\n";
    }
    
    if (_current_generation>1) {
        
        for (auto myGeneration=0; myGeneration<_current_generation; myGeneration++) {
            
            for (auto myAllele=0; myAllele<_num_alleles; myAllele++) {
                
                tmp_str += std::to_string(myGeneration) + "\t";
                tmp_str += std::to_string(myAllele) + "\t";
                tmp_str += std::to_string(_sim_data[myGeneration][myAllele]) + "\n";
            }
        }
    }
    
    return tmp_str;
}
