//
//  ActualDataFile.cpp
//  fits
//
//  Created by Tal Zinger on 16/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#include "ActualDataFile.hpp"


ActualDataFile::ActualDataFile() :
_actual_data(),
_is_initialized(false)
{}


ActualDataFile::ActualDataFile( const ActualDataFile& other ) :
_actual_data(other._actual_data),
_is_initialized(other._is_initialized)
{}


// based on code taken from the GUI
// 2016-04-12
void ActualDataFile::LoadActualData( std::string filename )
{
    _actual_data.clear();
    
    std::ifstream infile(filename);
    
    if (!infile.is_open()) {
        std::cerr << "error while opening actual data file: " << filename << std::endl;
        throw "error while opening actual data file";
    }
    
    std::string tmp_line;
    bool is_first_line = true;
    //std::vector<ActualDataEntry> tmp_actual_data_vec;
    
    
    // todo: check if newline needs to be normalized
    while (std::getline(infile, tmp_line)) {
        
        if (is_first_line) {
            is_first_line = false;
            continue;
        }
        
        
        std::vector<std::string> line_fields;
        try {
            boost::split(line_fields, tmp_line, boost::is_any_of("\t"));
            
        }
        catch (std::exception& e) {
            std::cerr << "Error while parsing actual data file to columns. Line: " << std::endl << tmp_line << std::endl;
            throw e;
        }
        catch (...) {
            std::cerr << "Unknown error while parsing actual data file to columns. Line: " << std::endl << tmp_line << std::endl;
            std::cerr << "NOT TERMINATING" << std::endl;
        }
        
        
        if (line_fields.size() <=1) {
            std::cout << "skipping line" << std::endl;
            continue;
        }
        
        // 2016-11-3 had enough of useleff information being read so not forcing all columns to be present
        // except generation, allele, frequency.
        if (line_fields.size() < ACTUAL_DATA_COLUMNS) {
            std::cerr << "incompatible actual data file with " << line_fields.size() << " columns, expected " << ACTUAL_DATA_COLUMNS << std::endl;
            throw "incompatible actual data file.";
        }
        
        /*
         if (line_fields.size() > ACTUAL_DATA_COLUMNS && warning_flag) {
         warning_flag = false;
         std::cout << "note: using only generation, frequency and allele from actual data file." << std::endl;
         }
         */
        
        ActualDataEntry tmp_data_entry;
        
        try {
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_GENERATION]);
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_ALLELE]);
            boost::trim(line_fields[ACTUAL_DATA_COLUMN_FREQ]);
        }
        catch (...) {
            std::cerr << "Error while parsing actual data file (trim)." << std::endl;
            throw "Error while pasring actual data file (trim).";
        }
        
        try {
            tmp_data_entry.gen = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_GENERATION]);
            tmp_data_entry.base = boost::lexical_cast<int>(line_fields[ACTUAL_DATA_COLUMN_ALLELE]);
            tmp_data_entry.freq = boost::lexical_cast<FLOAT_TYPE>(line_fields[ACTUAL_DATA_COLUMN_FREQ]);
            tmp_data_entry.pos = -1;
            tmp_data_entry.ref = -1;
            tmp_data_entry.read_count = -1;
        }
        catch (...) {
            std::cerr << "Error while parsing actual data file (cast)." << std::endl;
            throw "Error while pasring actual data file (cast).";
        }
        
        
        try {
            _actual_data.push_back(tmp_data_entry);
        }
        catch (...) {
            std::cerr << "Error while adding actual data entry" << std::endl;
            throw "Error while adding actual data entry";
        }
    }
    
    try {
        std::sort(_actual_data.begin(), _actual_data.end());
    }
    catch (...) {
        std::cerr << "Error while sorting actual data." << std::endl;
        throw "Error while sorting actual data.";
    }
    
    infile.close();
    
    std::cout << "Done." << std::endl;
    
    _is_initialized = true;
}


std::vector<int> ActualDataFile::GetActualGenerations()
{
    if (!_is_initialized) {
        std::cerr << "GetActualGenerations - object uninitialized" << std::endl;
    }
    
    std::vector<int> tmp_gen_vec;
    
    for (auto current_entry : _actual_data) {
        tmp_gen_vec.push_back(current_entry.gen);
    }
    
    // keep only unique values
    // sort because unique works on censecutive repeats
    // erase because no deletions occur but only running over existing data
    std::sort(tmp_gen_vec.begin(), tmp_gen_vec.end());
    auto last = std::unique(tmp_gen_vec.begin(), tmp_gen_vec.end());
    tmp_gen_vec.erase(last, tmp_gen_vec.end());
    
    return tmp_gen_vec;
}


// get vector of the frequencies of the alleles in the first generation
// assumes this is sorted
std::vector<FLOAT_TYPE> ActualDataFile::GetInitFreqs()
{
    if (!_is_initialized) {
        std::cerr << "GetInitFreqs - object uninitialized" << std::endl;
    }
    
    auto first_generation = _actual_data[0].gen;
    std::vector<FLOAT_TYPE> tmp_freqs(0);
    
    for (auto current_entry : _actual_data) {
        if (current_entry.gen == first_generation) {
            tmp_freqs.push_back(current_entry.freq);
        }
    }
    
    return tmp_freqs;
}

std::vector<FLOAT_TYPE> ActualDataFile::GetActualFrequencies()
{
    if (!_is_initialized) {
        std::cerr << "GetActualFrequencies - object uninitialized" << std::endl;
    }
    
    std::vector<FLOAT_TYPE> tmp_freq_vec;
    
    for (auto current_entry : _actual_data) {
        tmp_freq_vec.push_back(current_entry.freq);
    }
    
    return tmp_freq_vec;
}


int ActualDataFile::GetFirstGeneration()
{
    if (!_is_initialized) {
        std::cerr << "GetFirstGeneration - object uninitialized" << std::endl;
    }
    
    std::sort(_actual_data.begin(), _actual_data.end());
    auto ret_gen = _actual_data[0].gen;
    return ret_gen;
}


std::size_t ActualDataFile::GetNumberOfAlleles()
{
    std::size_t tmp_allele_count = 2;
    
    if ( _actual_data.empty() ) {
        throw "Get number of alleles - actual data vector is empty.";
    }
    
    if ( _actual_data.size() <= 1 ) {
        throw "Get number of alleles - actual data vector contains only 1 entry.";
    }
    
    std::sort(_actual_data.begin(), _actual_data.end());
    
    auto current_allele = _actual_data[0].base;
    
    while ( _actual_data[tmp_allele_count].base > _actual_data[tmp_allele_count].base ) {
        ++tmp_allele_count;
    }
    
    return tmp_allele_count;
}


std::vector<FLOAT_TYPE> ActualDataFile::GetSDPerAllele()
{
    auto num_alleles = GetNumberOfAlleles();
    
    std::vector<FLOAT_TYPE> result_vector(num_alleles);
    
    std::vector< boost::accumulators::accumulator_set<
    float,
    boost::accumulators::stats<
    boost::accumulators::tag::median,
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > acc_vec_sd( num_alleles );

    
    for ( auto current_entry : _actual_data ) {
        acc_vec_sd[ current_entry.base ]( current_entry.freq );
        
        std::cout << "allele " << current_entry.base << " freq " << current_entry.freq << std::endl;
    }
    
    for ( auto current_allele=0; current_allele<num_alleles; ++current_allele ) {
        result_vector[current_allele] = std::sqrt( boost::accumulators::variance(acc_vec_sd[current_allele]) );
    }
    
    return result_vector;
}
