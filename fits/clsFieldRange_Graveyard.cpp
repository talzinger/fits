//
//  clsFieldRange_Graveyard.cpp
//  fits
/* supposedly unused functions */

/*
 int clsFieldRange::UseSubRange( vector<float> min_limits, vector<float> max_limits )
 {
	_subrange_min_vector = min_limits;
	_subrange_max_vector = max_limits;
	_use_subrange = true;
	
	int new_combinations_in_subrange = 1;
	for ( int i=0; i<min_limits.size(); i++ ) {
 new_combinations_in_subrange *= static_cast<int>((max_limits[i] - min_limits[i]) / _subrange_increment + 1.0f);
	}
 
	int old_combinations_in_subrange = 1;
 
	for ( int i=0; i<min_limits.size(); i++ ) {
 new_combinations_in_subrange *= static_cast<int>((max_limits[i] - min_limits[i]) / _increment + 1.0f);
	}
	
	_combinations_with_subrange = _num_of_combinations + new_combinations_in_subrange - old_combinations_in_subrange;
	
	return (new_combinations_in_subrange - old_combinations_in_subrange);
 }
 */

/*
 string clsFieldRange::GetDataAsString( string separator = "," )
 {
	string tmpstr = to_string( _data_vector[0] );
	
	for (int i=1; i<_data_vector.size(); i++) {
 tmpstr += separator + to_string( _data_vector[i] );
	}
	
	return tmpstr;
 }
 */

/*
 void clsFieldRange::PrintCurrent()
 {
	for ( int i=0; i<_data_vector.size(); i++ ) {
 cout << _data_vector[i] << "\t";
	}
	
	cout << endl;
 }
 */

/*
 bool clsFieldRange::MaxReached()
 {
	return _max_reached;
 }
 */

/*
 vector<float> clsFieldRange::GetData()
 {
	return _data_vector;
 }
 */
/*
 float clsFieldRange::FractionCompleted()
 {
	return (float)_num_of_steps_passed / (float)_num_of_combinations;
 }
 */
/*
 int clsFieldRange::PercentCompleted()
 {
	// to avoid rounding issues limiting to 99
	if ( MaxReached() ) {
 return 100;
	}
	
	return _num_of_steps_passed*100 / _num_of_combinations;
 }
 */

/*
 void clsFieldRange::Reset()
 {
	_data_vector = _min_vector;
	_num_of_steps_passed = 0;
	_is_data_reset = true;
	_max_reached = false;
 }
 */

/*
 bool clsFieldRange::Next()
 {
	if (_combination_iterator == _all_combinations.end()) {
 return false;
	}
 
	_data_vector = *_combination_iterator;
	++_combination_iterator;
	return true;
 }
 */


/*
std::vector<MATRIX_TYPE> clsFieldRange::GetRandomCombinationMatrices( std::size_t num )
{
    std::vector<MATRIX_TYPE> vector_of_matrices(0);
    
    MATRIX_TYPE mtrx_delta(_min_matrix.size1(), _min_matrix.size2());
    
    std::uniform_real_distribution<FLOAT_TYPE> dist_uniform(0, 1);
    
    
    for ( auto i = 0; i < mtrx_delta.size1(); ++i ) {
        
        for ( auto j = 0; j < mtrx_delta.size2(); ++j ) {
            
            mtrx_delta(i, j) = _max_matrix(i, j) - _min_matrix(i, j);
        }
    }
    
    // add matrices to vector
    for ( auto count_matrx = 0; count_matrx < num; ++count_matrx ) {
        
        MATRIX_TYPE mtrx_mutation_rates(_min_matrix.size1(), _min_matrix.size2());
        
        // generate the random rates, (i,i) will be overwritten
        for ( auto i = 0; i < mtrx_mutation_rates.size1(); ++i ) {
            
            for ( auto j = 0; j < mtrx_mutation_rates.size2(); ++j ) {
                
                mtrx_mutation_rates(i, j) = _min_matrix(i, j) + mtrx_delta(i, j) * dist_uniform(_rnd_engine);
            }
        }
        
        // sum all rates in a row and put their complement to 1 in the diagonal
        for ( auto i = 0; i < mtrx_mutation_rates.size1(); ++i ) {
            
            FLOAT_TYPE sum_row = 0.0f;
            
            for ( auto j = 0; j < mtrx_mutation_rates.size2(); ++j ) {
                
                if ( i != j ) {
                    sum_row += mtrx_mutation_rates(i, j);
                }
            }
            
            mtrx_mutation_rates(i, i) = 1.0f - sum_row;
            
            std::cout << " sum is " << sum_row << " and complement is " << mtrx_mutation_rates(i, i) << std::endl;
            
            if ( mtrx_mutation_rates(i, i) < 0.0f ) {
                std::cerr <<  "GetRandomCombinationMatrices: negative complement mutation rate: " << mtrx_delta(i, i) << std::endl;
                throw "GetRandomCombinationMatrices: negative complement mutation rate.";
            }
        }
        
        
        std::cout << "clsFieldRange: adding matrix of mutation rates:" << std::endl;
        std::cout << mtrx_mutation_rates << std::endl;
        
        vector_of_matrices.push_back(std::move(mtrx_mutation_rates));
    }
    
    return vector_of_matrices;
}
*/
