/*
 void clsCMulatorABC::DivideEachAllele( std::size_t start_idx, std::size_t end_idx, std::vector<FLOAT_TYPE> value_vector )
 {
 
 for ( auto current_idx=start_idx; current_idx<end_idx; ++current_idx ) {
 
 //std::cout << "before: " << _simulation_result_vector[current_idx].sim_data_matrix << std::endl;
 
 for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
 
 boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_col( _simulation_result_vector[current_idx].sim_data_matrix, current_allele );
 
 current_col = current_col / value_vector[current_allele];
 }
 //std::cout << "after: " <<_simulation_result_vector[current_idx].sim_data_matrix << std::endl;
 }
 
 }
 */



/*
 OLD COMPOSITE DISTRIB
 
 // 2017-02-07 changed the parameters to improve the distribution visually
 // define distributions which compose the desired distribution
 boost::random::gamma_distribution<float> dist_gamma(1.0f, 0.1f); // dist_gamma(0.5f, 1.0f);
 boost::random::lognormal_distribution<float> dist_lognorm(0.0f, 2.0f); // (0.092, 6.0f)
 boost::random::normal_distribution<float> dist_norm(1.0f,0.2f);
 boost::random::discrete_distribution<int> dist_disc( { _fitness_lth_prob, _fitness_del_prob, _fitness_neu_prob, _fitness_adv_prob});
 
 
// 2017-02-07
// may cause some bias if the loop is outside so create small loops for each case
//do {
// first step - what category
auto chosen_category = dist_disc(_rnd_gen);

switch (chosen_category) {
        
    case fits_constants::FITNESS_ADV_IDX:
        do {
            tmp_val = FITNESS_NEU_VAL + dist_gamma(_rnd_gen);
        } while ( tmp_val > max || tmp_val < min );
        break;
        
    case fits_constants::FITNESSS_NEU_IDX:
        tmp_val = FITNESS_NEU_VAL;
        // tmp_fitness = 0.0 - std::fabs(dist_norm(_rnd_engine));
        break;
        
    case fits_constants::FITNESS_DEL_IDX:
        do {
            tmp_val = 1.0f - dist_lognorm(_rnd_gen);// * dist_norm(_rnd_gen);
        } while ( tmp_val > max || tmp_val < min );
        break;
        
    case fits_constants::FITNESS_LTH_IDX:
        tmp_val = 0.0f;
        break;
        
    default:
        break;
} // switch


//} while ( tmp_val > max || tmp_val < min );
*/


/*
 void SimulationResult::KeepOnlyActualGenerations()
 {
 if ( actual_generations.empty() ) {
 std::cerr << " Actual Generation list is empty" << std::endl;
 throw "Actual Generation list is empty";
 }
 
 auto num_alleles = fitness_values.size();
 auto rows_to_return = actual_generations.size();
 
 MATRIX_TYPE tmp_matrix(rows_to_return, num_alleles);
 
 for ( auto gen_idx=0; gen_idx<actual_generations.size(); ++gen_idx ) {
 
 auto shifted_gen = actual_generations[gen_idx] - generation_shift;
 
 boost::numeric::ublas::matrix_row<MATRIX_TYPE> current_src_row( sim_data_matrix, shifted_gen );
 boost::numeric::ublas::matrix_row<MATRIX_TYPE> current_dest_row( tmp_matrix, gen_idx );
 
 current_dest_row = current_src_row;
 }
 
 //std::cout << "original matrix: " << sim_data_matrix << std::endl;
 //std::cout << "new matrix: " << tmp_matrix << std::endl;
 
 sim_data_matrix = tmp_matrix;
 }
 */

/*
 std::string ReportStatsFitness(const std::vector<SimulationResult>& result_vector, std::string filename = "")
 {
 std::string output_str;
 
 if (result_vector.empty()) {
 std::cerr << "ReportStats: empty result vector." << std::endl;
 throw "ReportStats: empty result vector.";
 }
 
 std::ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);
 
 
 // open the file now, so if it'll fail we won't even start generating data
 if (!outfile.is_open() && !filename.empty() ) {
 std::cerr << "unable to open file for writing: " << filename << std::endl;
 throw "unable to open file for writing: " + filename;
 }
 
 output_str += "Fitness Report\n";
 output_str += "===============\n";
 
 //cout << "Processing " << result_vector.size() << " results." << endl;
 output_str += "Processing " + std::to_string(result_vector.size()) + " results.\n";
 
 // int _num_alleles = result_vector[0].fitness_values.size();
 
 ResultsStats results_stats;
 results_stats.CalculateStatsFitness(result_vector);
 
 for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
 output_str += "\tallele" + std::to_string(i);
 //cout << "\t" << "allele" << i;
 
 }
 output_str += "\n";
 //cout << endl;
 
 //cout << "median";
 
 // replacing boost's format so that we could put this inside a string
 //std::setprecision(3);
 
 output_str += "median";
 for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
 std::ostringstream tmp_out;
 tmp_out << std::fixed << std::setprecision(REPORT_FLOAT_PRECISION) << results_stats.allele_median_fitness[i];
 output_str += "\t" + tmp_out.str();
 
 // cout << "\t" << boost::format("%1.2f") % results_stats.allele_median_fitness[i];
 // output_str += "\t" + to_string(results_stats.allele_median_fitness[i]);
 }
 output_str += "\n";
 // cout << endl;
 
 output_str += "mean";
 // cout << "mean";
 for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
 std::ostringstream tmp_out;
 tmp_out << std::fixed << std::setprecision(REPORT_FLOAT_PRECISION) << results_stats.allele_mean_fitness[i];
 output_str += "\t" + tmp_out.str();
 
 // cout << "\t" << boost::format("%1.2f") % results_stats.allele_mean_fitness[i];
 // output_str += "\t" + to_string(results_stats.allele_mean_fitness[i]);
 }
 output_str += "\n";
 //cout << endl;
 
 //cout << "lbound";
 output_str += "lbound";
 for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
 std::ostringstream tmp_out;
 tmp_out << std::fixed << std::setprecision(REPORT_FLOAT_PRECISION) << results_stats.allele_min_fitness[i];
 output_str += "\t" + tmp_out.str();
 
 
 // output_str += "\t" + to_string(results_stats.allele_min_fitness[i]);
 // cout << "\t" << boost::format("%1.2f") % results_stats.allele_min_fitness[i];
 }
 //cout << endl;
 output_str += "\n";
 
 output_str += "ubound";
 //cout << "ubound";
 for (auto i = 0; i < result_vector[0].fitness_values.size(); i++) {
 std::ostringstream tmp_out;
 tmp_out << std::fixed << std::setprecision(REPORT_FLOAT_PRECISION) << results_stats.allele_max_fitness[i];
 output_str += "\t" + tmp_out.str();
 
 // output_str += "\t" + to_string(results_stats.allele_max_fitness[i]);
 // cout << "\t" << boost::format("%1.2f") % results_stats.allele_max_fitness[i];
 }
 //cout << endl;
 output_str += "\n";
 
 // distribution
 // todo: do this more generally with boost? only 4 bins so here it may be more work than benefit
 
 output_str += "\nDistribution:\n";
 // output_str += "\tlet.\tdel.\tneu.\tadv.\n";
 output_str += "\tdel.\tneu.\tadv.\n";
 
 for (int current_allele = 0; current_allele < results_stats._num_alleles; current_allele++) {
 
 output_str += "allele" + std::to_string(current_allele);
 //cout << "allele" << current_allele;
 
 // output_str += "\t" + to_string(results_stats.lethal_percent[current_allele]) + "%";
 output_str += "\t" + std::to_string(results_stats.deleterious_percent[current_allele]) + "%";
 output_str += "\t" + std::to_string(results_stats.neutral_percent[current_allele]) + "%";
 output_str += "\t" + std::to_string(results_stats.advantageous_percent[current_allele]) + "%";
 
 
 output_str += "\t=>";
 //cout << "\t=> ";
 
 switch ( results_stats.allele_category[current_allele] ) {
 
 case AlleleCategory::Undefined:
 output_str += "???";
 break;
 
 case AlleleCategory::Adventageous:
 output_str += "ADV";
 break;
 
 case AlleleCategory::Possible_advantageous:
 output_str += "?ADV";
 break;
 
 case AlleleCategory::Deleterious:
 output_str += "DEL";
 break;
 
 case AlleleCategory::Possible_deleterious:
 output_str += "?DEL";
 break;
 
 case AlleleCategory::Neutral:
 output_str += "NEU";
 break;
 case AlleleCategory::Possible_neutral:
 output_str += "?NEU";
 break;
 
 case AlleleCategory::Lethal:
 output_str += "LTH";
 break;
 
 case AlleleCategory::Possible_lethal:
 output_str += "?LTH";
 break;
 
 case AlleleCategory::WT:
 output_str += "WT";
 break;
 }
 
 output_str += "\n";
 }
 
 
 if ( outfile.is_open() ) {
 outfile << output_str;
 }
 
 std::cout << output_str;
 
 return output_str;
 }
 
 
 */


/*
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
 */
