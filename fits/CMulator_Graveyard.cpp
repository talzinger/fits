//
//  CMulator_Graveyard.cpp
//  fits
//
//  Created by Tal Zinger on 19/10/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//


/*
 CMulator::CMulator( std::string param_filename ) :
 _sim_data_matrix()
 {
 auto zparams = CMulatorToolbox::ReadParametersFromFile(param_filename);
 InitMemberVariables(zparams);
 _initialized_with_parameters = true;
	
	// InitMemberVariables();
 _time_for_seeding = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
 _boost_gen.seed(_time_for_seeding);
 }
 */


/*
 void CMulator::Reset_Hard()
 {
	if ( Parameters::empty() ) {
 throw " CMulator Hard Reset unsuccessful - no loaded parameters detected.";
	}
	
	
	// if parameters were missing, an exception would have prevented us from getting here
	_initialized_with_parameters = true;
	
	InitMemberVariables();
 }
 */





/*
 void CMulator::Reset_Hard()
 {
	if ( Parameters::empty() ) {
 throw " CMulator Hard Reset unsuccessful - no loaded parameters detected.";
	}
	
	
	// if parameters were missing, an exception would have prevented us from getting here
	_initialized_with_parameters = true;
	
	InitMemberVariables();
 }
 */


/*
 void CMulator::InitializeFromParamFile(string param_filename)
 {
	if (_initialized_with_parameters) {
 cerr << "CMulator InitializeFromParamFile: object already initialized" << endl;
 throw "CMulator InitializeFromParamFile: object already initialized";
	}
 
	// essentially ctor with param file
	ReadParametersFromFile(param_filename);
	_initialized_with_parameters = true;
 
	InitMemberVariables();
 }
 */


/*
string CMulator::GetAllOutputForR( bool zero_based_alleles, bool print_header, bool extended_header, int genshift, int resample )
{
    stringstream R_string;
    
    if (!_initialized_with_parameters) {
        throw " PrintDataForR: object not initialized with parameters.";
    }
    
    // print information regarding the simulation
    // this would be ignored by R as it is commented
    if (extended_header) {
        R_string << "# Simulator Report" << std::endl;
        R_string << "# ====================" << std::endl;
        R_string << "# N = " << _N << std::endl;
        R_string << "# Alleles = " << _num_alleles << std::endl;
        R_string << "# Generations = " << _num_generations << std::endl;
        R_string << "# Sample Size = " << _sample_size << std::endl;
        
        R_string << "# Initial frequencies: (" << _allele_init_freqs.size() << ")" << endl;
        for (auto tmpi=0; tmpi<_allele_init_freqs.size(); tmpi++) {
            R_string << "# \t" << tmpi << ": " << _allele_init_freqs[tmpi] << endl;
        }
        
        R_string << "# Fitness values:" << endl;
        for (auto tmpi=0; tmpi<_allele_fitness.size(); tmpi++) {
            R_string << "# \t" << tmpi << ": " << _allele_fitness[tmpi] << endl;
        }
        
        R_string << "# Mutation Rates:" << endl;
        for (auto from_allele=0; from_allele<_num_alleles; from_allele++) {
            
            R_string << "# ";
            
            for (auto to_allele=0; to_allele<_num_alleles; to_allele++) {
                
                R_string << _mutation_rates[from_allele][to_allele] << "\t";
            }
            
            R_string << endl;
        }
    }
    
    
    // header
    if (print_header) {
        
        // yes, this separation of tabs and titles is redundant, but is more readable in a sense
        R_string << "sim_uid"
        << "\t" << "gen"
        << "\t" << "base"
        << "\t" << "freq"
        << "\t" << "init_freq"
        << "\t" << "allele_fitness";
        
        
 
        R_string << std::endl;
        
    }
    
    // dump data
    for (int gen_count=0; gen_count<=_num_generations; gen_count++) {
        
        // non-invasive sampling
        MDOUBLE currentSum=0, currentP=0;
        vector<MDOUBLE> newvals(_num_alleles);
        //		MDOUBLE newvals[4] = {0,0,0,0};
        
        if (resample) {
            // randomally select individuals possessing each allele
            for ( int allele_i=0; allele_i<_num_alleles; allele_i++ ) {
                
                // from float to MDOUBLE
                currentP = (MDOUBLE)_sim_data[gen_count][allele_i];
                
                // this returns number of actual individuals chosen - practically an integer
                //newvals[allele_i] = (int)clsBinomialDistribution::SampleFromBinomialDistribution(currentP, _N, time_for_seeding);
                newvals[allele_i] = (int)bnldev(currentP, _sample_size, time_for_seeding);
                
                // actually selected individuals, as each is sampled separately, total samples will probably be != N
                currentSum += newvals[allele_i];
            }
            
            // Normalize frequencies
            for ( auto allele_i=0; allele_i<_num_alleles; allele_i++ ) {
                newvals[allele_i] = newvals[allele_i] / currentSum;
            }
        } else {
            for ( auto allele_i=0; allele_i<_num_alleles; allele_i++ ) {
                newvals[allele_i] = (MDOUBLE)_sim_data[gen_count][allele_i];
            }
        }
        
        
        
        
        // now write data for each alleles in this repeat and generation
        for (auto my_allele=0; my_allele<_num_alleles; my_allele++) {
            
            R_string << _uid
            << "\t" << gen_count+genshift
            << "\t" << (zero_based_alleles ? my_allele : my_allele+1 )
            << "\t" << newvals[my_allele]
            //				<< "\t" << (gen_count==0 ? _sim_data[gen_count][my_allele] : newvals[my_allele] ) // do not random-sample, so that 1st gen would be the same
            << "\t" << _sim_data[0][my_allele]
            << "\t" << _allele_fitness[my_allele];
            
            
            R_string << std::endl;
        }
        
        // some files used to have some random numbers in the last line of the table
        // I thought this would cure the problem, it seems like it has
        //R_string.flush();
        
    }
    
    return R_string.str();
}
*/

/*
void CMulator::PrintDataAsMatrix(string seperator) const
{
    if (!_initialized_with_parameters) {
        throw " PrintDataAsMatrix: object not initialized with parameters.";
    }
    
    // header
    for (int myAllele=0; myAllele<_num_alleles; myAllele++) {
        cout << "Allele " << myAllele << "\t";
    }
    
    cout << endl;
    
    
    if (_current_generation>1) {
        
        for (int myGeneration=0; myGeneration<_current_generation; myGeneration++) {
            
            for (int myAllele=0; myAllele<_num_alleles; myAllele++) {
                
                cout << _sim_data[myGeneration][myAllele] << "\t";
            }
            
            cout << endl;
        }
    }
    
}
*/

/*
void CMulator::PrintSingleAllele( string filename, int allele, bool use_newline )
{
    if (!_initialized_with_parameters) {
        throw " PrintSingleAllele: object not initialized with parameters.";
    }
    
    std::ofstream outfile (filename, std::ofstream::out | std::ofstream::app);
    
    if ( outfile.is_open()) {
        
        for (int i=0; i<=_num_generations; i++) {
            
            outfile << _sim_data[i][allele] << (use_newline ? "\n" : "\t");
        }
    }
    outfile.close();
}
*/

/*
 void CMulator::PrintDataForR( const std::string& filename, bool append_file, bool zero_based_alleles, bool print_header, bool extended_header, int genshift, int resample )
 {
	if (!_initialized_with_parameters) {
 throw " PrintDataForR: object not initialized with parameters.";
	}
	
	std::ofstream R_outfile;
	
	if (append_file) {
 R_outfile.open( filename, ofstream::out | ofstream::app);
	}
	else {
 // erase previous content
 R_outfile.open( filename, ofstream::out | ofstream::trunc);
	}
	
	if ( R_outfile.is_open()) {
 
 
 // print information regarding the simulation
 // this would be ignored by R as it is commented
 if (extended_header) {
 R_outfile << "# Simulator Report" << std::endl;
 R_outfile << "# ====================" << std::endl;
 R_outfile << "# N = " << _N << std::endl;
 R_outfile << "# Alleles = " << _num_alleles << std::endl;
 R_outfile << "# Generations = " << _num_generations << std::endl;
 R_outfile << "# Sample Size = " << _sample_size << std::endl;
 
 R_outfile << "# Initial frequencies: (" << _allele_init_freqs.size() << ")" << endl;
 for (int tmpi=0; tmpi<_allele_init_freqs.size(); tmpi++) {
 R_outfile << "# \t" << tmpi << ": " << _allele_init_freqs[tmpi] << endl;
 }
 
 R_outfile << "# Fitness values:" << endl;
 for (int tmpi=0; tmpi<_allele_fitness.size(); tmpi++) {
 R_outfile << "# \t" << tmpi << ": " << _allele_fitness[tmpi] << endl;
 }
 
 R_outfile << "# Mutation Rates:" << endl;
 for (int from_allele=0; from_allele<_num_alleles; from_allele++) {
 
 R_outfile << "# ";
 
 for (int to_allele=0; to_allele<_num_alleles; to_allele++) {
 
 R_outfile << _mutation_rates[from_allele][to_allele] << "\t";
 }
 
 R_outfile << endl;
 }
 }
 
 
 // header
 if (print_header) {
 
 // yes, this separation of tabs and titles is redundant, but is more readable in a sense
 R_outfile << "sim_uid"
 << "\t" << "gen"
 << "\t" << "base"
 << "\t" << "freq"
 << "\t" << "init_freq"
 << "\t" << "allele_fitness";
 
 

R_outfile << std::endl;

}
*/
/*
void CMulator::PrintDataForR( const std::string& filename, bool append_file, bool zero_based_alleles, bool print_header, bool extended_header, int genshift, int resample )
{
    if (!_initialized_with_parameters) {
        throw " PrintDataForR: object not initialized with parameters.";
    }
    
    std::ofstream R_outfile;
    
    if (append_file) {
        R_outfile.open( filename, ofstream::out | ofstream::app);
    }
    else {
        // erase previous content
        R_outfile.open( filename, ofstream::out | ofstream::trunc);
    }
    
    if ( R_outfile.is_open()) {
        
        
        // print information regarding the simulation
        // this would be ignored by R as it is commented
        if (extended_header) {
            R_outfile << "# Simulator Report" << std::endl;
            R_outfile << "# ====================" << std::endl;
            R_outfile << "# N = " << _N << std::endl;
            R_outfile << "# Alleles = " << _num_alleles << std::endl;
            R_outfile << "# Generations = " << _num_generations << std::endl;
            R_outfile << "# Sample Size = " << _sample_size << std::endl;
            
            R_outfile << "# Initial frequencies: (" << _allele_init_freqs.size() << ")" << endl;
            for (int tmpi=0; tmpi<_allele_init_freqs.size(); tmpi++) {
                R_outfile << "# \t" << tmpi << ": " << _allele_init_freqs[tmpi] << endl;
            }
            
            R_outfile << "# Fitness values:" << endl;
            for (int tmpi=0; tmpi<_allele_fitness.size(); tmpi++) {
                R_outfile << "# \t" << tmpi << ": " << _allele_fitness[tmpi] << endl;
            }
            
            R_outfile << "# Mutation Rates:" << endl;
            for (int from_allele=0; from_allele<_num_alleles; from_allele++) {
                
                R_outfile << "# ";
                
                for (int to_allele=0; to_allele<_num_alleles; to_allele++) {
                    
                    R_outfile << _mutation_rates[from_allele][to_allele] << "\t";
                }
                
                R_outfile << endl;
            }
        }
        
        
        // header
        if (print_header) {
            
            // yes, this separation of tabs and titles is redundant, but is more readable in a sense
            R_outfile << "sim_uid"
            << "\t" << "gen"
            << "\t" << "base"
            << "\t" << "freq"
            << "\t" << "init_freq"
            << "\t" << "allele_fitness";
            
 
            
            R_outfile << std::endl;
            
        }
        
        // dump data
        for (int gen_count=0; gen_count<=_num_generations; gen_count++) {
            
            // non-invasive sampling
            MDOUBLE currentSum=0, currentP=0;
            //			MDOUBLE newvals[4] = {0,0,0,0};
            vector<MDOUBLE> newvals(_num_alleles);
            
            if (resample) {
                // randomally select individuals possessing each allele
                for ( int allele_i=0; allele_i<_num_alleles; allele_i++ ) {
                    
                    // from float to MDOUBLE
                    currentP = (MDOUBLE)_sim_data[gen_count][allele_i];
                    
                    // this returns number of actual individuals chosen - practically an integer
                    newvals[allele_i] = (int)bnldev(currentP, _sample_size, time_for_seeding);
                    //newvals[allele_i] = (int)clsBinomialDistribution::SampleFromBinomialDistribution(currentP, _sample_size, time_for_seeding);
                    
                    // actually selected individuals, as each is sampled separately, total samples will probably be != N
                    currentSum += newvals[allele_i];
                }
                
                // Normalize frequencies
                for ( int allele_i=0; allele_i<_num_alleles; allele_i++ ) {
                    newvals[allele_i] = newvals[allele_i] / currentSum;
                }
            } else {
                for ( int allele_i=0; allele_i<_num_alleles; allele_i++ ) {
                    newvals[allele_i] = (MDOUBLE)_sim_data[gen_count][allele_i];
                }
            }
            
            
            
            // now write data for each alleles in this repeat and generation
            for (int my_allele=0; my_allele<_num_alleles; my_allele++) {
                
                R_outfile << _uid
                << "\t" << gen_count+genshift
                << "\t" << (zero_based_alleles ? my_allele : my_allele+1 )
                << "\t" << newvals[my_allele]
                //				<< "\t" << (gen_count==0 ? _sim_data[gen_count][my_allele] : newvals[my_allele] ) // do not random-sample, so that 1st gen would be the same
                << "\t" << _sim_data[0][my_allele]
                << "\t" << _allele_fitness[my_allele];
 
                R_outfile << std::endl;
            }
            
            // some files used to have some random numbers in the last line of the table
            // I thought this would cure the problem, it seems like it has
            R_outfile.flush();
            
        }
        
        R_outfile.close();
    }
    
}
 */

/*
void CMulator::PrintReport( ReportType report_type )
{
    if (!_initialized_with_parameters) {
        throw " PrintReport: object not initialized with parameters.";
    }
    
    switch (report_type) {
        case FULL:
            cout << "Simulator Report" << std::endl;
            cout << "====================" << std::endl;
            cout << "N = " << _N << std::endl;
            cout << "Alleles = " << _num_alleles << std::endl;
            cout << "Generations = " << _num_generations << std::endl;
            cout << "Sample Size = " << _sample_size << std::endl;
            
            cout << "Initial frequencies: (" << _allele_init_freqs.size() << ")" << endl;
            for (int tmpi=0; tmpi<_allele_init_freqs.size(); tmpi++) {
                cout << "\t" << tmpi << ": " << _allele_init_freqs[tmpi] << endl;
            }
            
            cout << "Fitness values:" << endl;
            for (int tmpi=0; tmpi<_allele_fitness.size(); tmpi++) {
                cout << "\t" << tmpi << ": " << _allele_fitness[tmpi] << endl;
            }
            
            cout << "Mutation Rates:" << endl;
            for (int from_allele=0; from_allele<_num_alleles; from_allele++) {
                
                for (int to_allele=0; to_allele<_num_alleles; to_allele++) {
                    
                    cout << _mutation_rates[from_allele][to_allele] << "\t";
                }
                
                cout << endl;
            }
            
            // intentionally leave out break so that allele frequencies will be printed
            // break;
            
            
        case FREQS_ONLY:
            if (_current_generation>1) {
                
                for (int myAllele=0; myAllele<_num_alleles; myAllele++) {
                    cout << "Allele " << myAllele << "\t";
                }
                
                cout << endl;
                
                for (auto myGeneration=0; myGeneration<_current_generation; myGeneration++) {
                    
                    for (auto myAllele=0; myAllele<_num_alleles; myAllele++) {
                        
                        cout << _sim_data[myGeneration][myAllele] << "\t";
                    }
                    
                    cout << endl;
                }
            }
            break;
            
        case FITNESS_ONLY:
            cout << "Simulation " << _uid;
            
            for (int i=0; i<_num_alleles; i++) {
                cout << "\tallele " << i << ": " << _allele_fitness[i];
            }
            cout << endl;
            break;
    }
}
*/

/*
// dump data
for (int gen_count=0; gen_count<=_num_generations; gen_count++) {
    
    // non-invasive sampling
    MDOUBLE currentSum=0, currentP=0;
    //			MDOUBLE newvals[4] = {0,0,0,0};
    vector<MDOUBLE> newvals(_num_alleles);
    
    if (resample) {
        // randomally select individuals possessing each allele
        for ( int allele_i=0; allele_i<_num_alleles; allele_i++ ) {
            
            // from float to MDOUBLE
            currentP = (MDOUBLE)_sim_data[gen_count][allele_i];
            
            // this returns number of actual individuals chosen - practically an integer
            newvals[allele_i] = (int)bnldev(currentP, _sample_size, time_for_seeding);
            //newvals[allele_i] = (int)clsBinomialDistribution::SampleFromBinomialDistribution(currentP, _sample_size, time_for_seeding);
            
            // actually selected individuals, as each is sampled separately, total samples will probably be != N
            currentSum += newvals[allele_i];
        }
        
        // Normalize frequencies
        for ( int allele_i=0; allele_i<_num_alleles; allele_i++ ) {
            newvals[allele_i] = newvals[allele_i] / currentSum;
        }
    } else {
        for ( int allele_i=0; allele_i<_num_alleles; allele_i++ ) {
            newvals[allele_i] = (MDOUBLE)_sim_data[gen_count][allele_i];
        }
    }
    
    
    
    // now write data for each alleles in this repeat and generation
    for (int my_allele=0; my_allele<_num_alleles; my_allele++) {
        
        R_outfile << _uid
        << "\t" << gen_count+genshift
        << "\t" << (zero_based_alleles ? my_allele : my_allele+1 )
        << "\t" << newvals[my_allele]
        //				<< "\t" << (gen_count==0 ? _sim_data[gen_count][my_allele] : newvals[my_allele] ) // do not random-sample, so that 1st gen would be the same
        << "\t" << _sim_data[0][my_allele]
        << "\t" << _allele_fitness[my_allele];
        
        
        R_outfile << std::endl;
    }
    
    // some files used to have some random numbers in the last line of the table
    // I thought this would cure the problem, it seems like it has
    R_outfile.flush();
    
}

R_outfile.close();
}

}
*/
/*
void CMulator::PrintReport( ReportType report_type )
{
    if (!_initialized_with_parameters) {
        throw " PrintReport: object not initialized with parameters.";
    }
    
    switch (report_type) {
        case FULL:
            cout << "Simulator Report" << std::endl;
            cout << "====================" << std::endl;
            cout << "N = " << _N << std::endl;
            cout << "Alleles = " << _num_alleles << std::endl;
            cout << "Generations = " << _num_generations << std::endl;
            cout << "Sample Size = " << _sample_size << std::endl;
            
            cout << "Initial frequencies: (" << _allele_init_freqs.size() << ")" << endl;
            for (int tmpi=0; tmpi<_allele_init_freqs.size(); tmpi++) {
                cout << "\t" << tmpi << ": " << _allele_init_freqs[tmpi] << endl;
            }
            
            cout << "Fitness values:" << endl;
            for (int tmpi=0; tmpi<_allele_fitness.size(); tmpi++) {
                cout << "\t" << tmpi << ": " << _allele_fitness[tmpi] << endl;
            }
            
            cout << "Mutation Rates:" << endl;
            for (int from_allele=0; from_allele<_num_alleles; from_allele++) {
                
                for (int to_allele=0; to_allele<_num_alleles; to_allele++) {
                    
                    cout << _mutation_rates[from_allele][to_allele] << "\t";
                }
                
                cout << endl;
            }
            
            // intentionally leave out break so that allele frequencies will be printed
            // break;
            
            
        case FREQS_ONLY:
            if (_current_generation>1) {
                
                for (int myAllele=0; myAllele<_num_alleles; myAllele++) {
                    cout << "Allele " << myAllele << "\t";
                }
                
                cout << endl;
                
                for (auto myGeneration=0; myGeneration<_current_generation; myGeneration++) {
                    
                    for (auto myAllele=0; myAllele<_num_alleles; myAllele++) {
                        
                        cout << _sim_data[myGeneration][myAllele] << "\t";
                    }
                    
                    cout << endl;
                }
            }
            break;
            
        case FITNESS_ONLY:
            cout << "Simulation " << _uid;
            
            for (int i=0; i<_num_alleles; i++) {
                cout << "\tallele " << i << ": " << _allele_fitness[i];
            }
            cout << endl;
            break;
    }
}
*/

/*
 int CMulator::IsFirstInBatch()
 {
	if (!_initialized_with_parameters) {
 throw " IsFirstInBatch: object not initialized with parameters.";
	}
	
	return _first_sim_in_batch;
 }
 */

/*
 void CMulator::UpdateParameters()
 {
	// There is only one place where this could be good - the prototype simulator object.
	//Because many instances actually access this, updates may corrupt everything.
	//That is why I'm commenting it all. 2016-05-26
 
 
	Parameters::updateParameter( PARAM_R_FILENAME, 		_R_filename.c_str() );
	
	Parameters::updateParameter( PARAM_POPULATION_SIZE, 		to_string(_N).c_str() );
	
	Parameters::updateParameter( PARAM_EPSILON_FLOAT_COMPARE, 	to_string(_epsilon_float_compare).c_str() );
	Parameters::updateParameter( PARAM_NUM_GENERATIONS, 		to_string(_num_generations).c_str() );
	Parameters::updateParameter( PARAM_SAMPLE_SIZE, 			to_string(_sample_size ).c_str() );
	
	Parameters::updateParameter( PARAM_SIM_ID, 					_uid.c_str() );
	
	Parameters::updateParameter( PARAM_NUM_ALLELES, 			to_string(_num_alleles).c_str() );
	
	Parameters::updateParameter( PARAM_BOTTLENECK_INTERVAL, 	to_string(_bottleneck_interval).c_str() );
	Parameters::updateParameter( PARAM_BOTTLENECK_SIZE, 		to_string(_bottleneck_size ).c_str() );
	
	Parameters::updateParameter( PARAM_LOGISTIC_GROWTH, 		to_string(_logistic_growth).c_str() );
	Parameters::updateParameter( PARAM_LOGISTIC_GROWTH_K, 		to_string(_logistic_growth_K).c_str() );
	Parameters::updateParameter( PARAM_LOGISTIC_GROWTH_r, 		to_string(_logistic_growth_r).c_str() );
	
	Parameters::updateParameter( PARAM_GENERATION_SHIFT, 		to_string(_generation_shift ).c_str() );
	Parameters::updateParameter( PARAM_FIRST_SIM_IN_BATCH, 		to_string(_first_sim_in_batch).c_str() );
	
	Parameters::updateParameter( PARAM_ALLELE_NO_INIT_FREQS, 	to_string(_no_init_freqs_as_parameters).c_str() );
	
	Parameters::updateParameter( PARAM_FIXED_UPPER_THRESHOLD, 	to_string(_fixation_upper_threshold ).c_str() );
	Parameters::updateParameter( PARAM_FIXED_LOWER_THRESHOLD, 	to_string(_fixation_lower_threshold ).c_str() );
	Parameters::updateParameter( PARAM_FIXED_CONSECUTIVE_GENERATIONS, to_string(_fixation_consecutive_generations).c_str() );
	
	
	// update only if previously existed
	if ( Parameters::getInt( PARAM_ALLELE_NO_INIT_FREQS, -1 ) > 0 ) {
 
 for ( int allele_count=0; allele_count<_num_alleles; allele_count++ ) {
 
 string param_name = PARAM_ALLELE_INIT_FREQ + to_string(allele_count);
 string param_val = to_string(_allele_init_freqs[allele_count]);
 
 Parameters::updateParameter( param_name, param_val.c_str() );
 }
	}
	
	
	if ( !_allele_max_fitness.empty() ) {
 for ( int allele_count=0; allele_count<_num_alleles; allele_count++ ) {
 
 string param_name = PARAM_ALLELE_MAX_FITNESS + to_string(allele_count);
 string param_val = to_string(_allele_max_fitness[allele_count]);
 
 Parameters::updateParameter( param_name, param_val.c_str() );
 }
	}
	
	if ( !_allele_min_fitness.empty() ) {
 for ( int allele_count=0; allele_count<_num_alleles; allele_count++ ) {
 
 string param_name = PARAM_ALLELE_MIN_FITNESS + to_string(allele_count);
 string param_val = to_string(_allele_min_fitness[allele_count]);
 
 Parameters::updateParameter( param_name, param_val.c_str() );
 }
	}
	
	if ( _allele_max_fitness.empty() && _allele_min_fitness.empty() && !_allele_fitness.empty() ) {
 for ( int allele_count=0; allele_count<_num_alleles; allele_count++ ) {
 
 string param_name = PARAM_ALLELE_FITNESS + to_string(allele_count);
 string param_val = to_string(_allele_fitness[allele_count]);
 
 Parameters::updateParameter( param_name, param_val.c_str() );
 }
	}
 
 
	for ( int i=0; i<_num_alleles; i++ ) {
 
 for ( int j=0; j<_num_alleles; j++ ) {
 
 string param_name = PARAM_MUTATION_RATE + to_string(i) + "_" + to_string(j);
 string param_val = to_string( _mutation_rates[i][j] );
 
 Parameters::updateParameter( param_name, param_val.c_str() );
 }
	}
	
	
 }
 */

