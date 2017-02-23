//
//  CMulator_Utilities.cpp
//  fits
//
//  Created by Tal Zinger on 15/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "CMulator.h"

void CMulator::Reset_Soft()
{
    if ( _allele_init_freqs.size() != _num_alleles ) {
        throw " Reset soft: number of alleles changed unexpectedly, consider hard reset.";
    }
    
    if ( _allele_init_freqs.size() != _num_alleles ) {
        throw " Reset soft: number of alleles changed unexpectedly, consider hard reset.";
    }
    
    
    // for consistency, update initial frequencies
    for (auto j=0; j<_num_alleles; j++ ) {
        
        _sim_data[0][j] = _allele_init_freqs[j];
        _sim_data_matrix(0,j) = _allele_init_freqs[j];
    }
    
    // now zero all others
    for (auto i=1; i<_num_generations; i++ ) {
        
        for (auto j=0; j<_num_alleles; j++ ) {
            
            _sim_data[i][j] = 0;
            _sim_data_matrix(i,j) = 0.0f;
        }
    }
    
    _current_generation = 0;
}


void CMulator::Reset_Soft( std::string new_sim_uid )
{
    Reset_Soft();
    SetSimUID( new_sim_uid );
}
