/* CMulatorABC Graveyard
supposedly unused functions */
/*
// Run simulations in serial mode
vector<CMulatorToolbox::SimulationResult> clsCMulatorABC::RunSimulationsSerial(int num_repeats, bool use_threshold)
{
    // vector<CMulatorToolbox::SimulationResult> tmp_vec1;
    
    auto tmp_vec2 = RunFitnessInferenceBatch(num_repeats, use_threshold);
    
    return tmp_vec2;
}
*/

/*
CMulatorToolbox::SimulationResult clsCMulatorABC::RunSingleSimulation( CMulator my_sim )
{
    my_sim.Reset_Soft();
    my_sim.EvolveAllGenerations();
    
    CMulatorToolbox::SimulationResult sim_result(my_sim);
    
    // if generations not selected, remain with the default
    if (_selected_actual_generations.size() > 0) {
        auto raw_frequency_data = my_sim.GetRawFrequencyData(_selected_actual_generations);
        
        // calculate distance instead of default
        if (_actual_data_raw_freqs.size() > 0) {
            sim_result.distance_from_actual = CMulatorToolbox::GetDistanceSimActual(_actual_data_raw_freqs, raw_frequency_data);
        }
    }
    
    
    return sim_result;
}
*/

/*
CMulatorToolbox::SimulationResult clsCMulatorABC::RunSingleSimulation()
{
    _sim_object.Reset_Soft();
    _sim_object.EvolveAllGenerations();
    
    CMulatorToolbox::SimulationResult sim_result(_sim_object);
    
    // if generations not selected, remain with the default
    if (_selected_actual_generations.size() > 0) {
        auto raw_frequency_data = _sim_object.GetRawFrequencyData(_selected_actual_generations);
        
        // calculate distance instead of default
        if (_actual_data_raw_freqs.size() > 0) {
            sim_result.distance_from_actual = CMulatorToolbox::GetDistanceSimActual(_actual_data_raw_freqs, raw_frequency_data);
        }
    }
    
    
    return sim_result;
}
*/

/*
 clsCMulatorABC::clsCMulatorABC(CMulator pre_initialized_sim_object, clsFieldRange fitness_range, std::vector<int> selected_actual_generations, std::vector<FLOAT_TYPE> actual_raw_freqs)
	: _rejection_threshold(0.0f),
	_sim_object(pre_initialized_sim_object),
	_actual_data_raw_freqs(actual_raw_freqs),
	_selected_actual_generations(selected_actual_generations),
	_fitness_range(fitness_range),
	_repeats(1)
 {
 if ( _actual_data_raw_freqs.empty() ) {
 throw std::range_error("Actual data vector is empty");
 }
 
 if ( _selected_actual_generations.empty() ) {
 throw std::range_error("No selected actual generations");
 }
 }
 */
