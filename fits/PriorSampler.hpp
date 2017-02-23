//
//  PriorSampler.hpp
//  fits
//
//  Created by Tal Zinger on 14/11/2016.
//  Copyright © 2016 Stern Lab. All rights reserved.
//

#ifndef PriorSampler_hpp
#define PriorSampler_hpp

#include <iostream>
#include <vector>
#include <chrono>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>

#include <limits>

#include "CMulator.h"
#include "ZParams.h"


enum PriorDistributionType {
    UNIFORM,
    LOG_UNIFORM,
    FITNESS_COMPOSITE
};

template<class C>
class PriorSampler {
public:
    void SetManualCategoryProportions(ZParams zparams)
    {
        bool at_least_one_found = false;
        
        auto manual_lth = zparams.GetFloat( "_prior_fraction_lth", -1.0f );
        auto manual_del = zparams.GetFloat( "_prior_fraction_del", -1.0f );
        auto manual_neu = zparams.GetFloat( "_prior_fraction_neu", -1.0f );
        auto manual_adv = zparams.GetFloat( "_prior_fraction_adv", -1.0f );
        
        at_least_one_found = at_least_one_found || (manual_lth >= 0.0);
        at_least_one_found = at_least_one_found || (manual_del >= 0.0);
        at_least_one_found = at_least_one_found || (manual_neu >= 0.0);
        at_least_one_found = at_least_one_found || (manual_adv >= 0.0);
        
        if (at_least_one_found) {
            auto all_found = (manual_del >= 0.0) && (manual_neu >= 0.0) && (manual_adv >= 0.0);
            
            if ( !all_found) {
                std::cerr << "Prior fraction LTH : " << ( manual_lth >= 0 ? "FOUND" : "MISSING" ) << std::endl;
                std::cerr << "Prior fraction DEL : " << ( manual_del >= 0 ? "FOUND" : "MISSING" ) << std::endl;
                std::cerr << "Prior fraction NEU : " << ( manual_neu >= 0 ? "FOUND" : "MISSING" ) << std::endl;
                std::cerr << "Prior fraction ADV : " << ( manual_adv >= 0 ? "FOUND" : "MISSING" ) << std::endl;
                
                throw "Coomposite prior: At least one but not all manual fractions are given.";
            }
            
            auto tmp_sum = manual_lth + manual_del + manual_neu + manual_adv;
            
            if (tmp_sum - 1.0f > std::numeric_limits<float>::epsilon() ) {
                std::cerr << "Prior fractions LTH/DEL/NEU/ADV do not sum up to 1.0 (" << tmp_sum << ")" << std::endl;
                throw "Prior fractions LTH/DEL/NEU/ADV do not sum up to 1.0";
            }
            
            _fitness_lth_prob = manual_lth;
            _fitness_del_prob = manual_del;
            _fitness_neu_prob = manual_neu;
            _fitness_adv_prob = manual_adv;
        }
        
        // std::cout << "Prior updated." << std::endl;
    }
                                 
private:
    std::vector<C> _min_vector;
    std::vector<C> _max_vector;
    PriorDistributionType _distrib_type;
    
    boost::mt19937 _rnd_gen;
    unsigned int _rnd_seed;
    
    // define fitness category index, to be used with the discrete distribution
    static const int FITNESS_LTH_IDX = 0;
    static const int FITNESS_DEL_IDX = 1;
    static const int FITNESSS_NEU_IDX = 2;
    static const int FITNESS_ADV_IDX = 3;
    
    // define probability for each fitness category
    const float FITNESS_LTH_PROB = 0.1;
    const float FITNESS_DEL_PROB = 0.79;
    const float FITNESS_NEU_PROB = 0.1;
    const float FITNESS_ADV_PROB = 0.01;
    
    float _fitness_lth_prob;
    float _fitness_del_prob;
    float _fitness_neu_prob;
    float _fitness_adv_prob;
    
    const float FITNESS_NEU_VAL = 1.0;
    
    
public:
    PriorSampler() :
    _min_vector(0),
    _max_vector(0),
    _rnd_gen(0),
    _fitness_lth_prob(FITNESS_LTH_PROB),
    _fitness_del_prob(FITNESS_DEL_PROB),
    _fitness_neu_prob(FITNESS_NEU_PROB),
    _fitness_adv_prob(FITNESS_ADV_PROB),
    _distrib_type(PriorDistributionType::UNIFORM)
    {
        _rnd_seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        _rnd_gen.seed(_rnd_seed);
    }
    
    PriorSampler( std::vector<C> min_values, std::vector<C> max_values, PriorDistributionType distrib_type ) :
    _min_vector(min_values),
    _max_vector(max_values),
    _distrib_type(distrib_type),
    _rnd_gen(0),
    _fitness_lth_prob(FITNESS_LTH_PROB),
    _fitness_del_prob(FITNESS_DEL_PROB),
    _fitness_neu_prob(FITNESS_NEU_PROB),
    _fitness_adv_prob(FITNESS_ADV_PROB)
    {
        _rnd_seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        _rnd_gen.seed(_rnd_seed);
    }
    
    PriorSampler( PriorSampler<C> &original ) :
    _min_vector(original._min_vector),
    _max_vector(original._max_vector),
    _rnd_gen(0),
    _distrib_type(PriorDistributionType::UNIFORM),
    _fitness_lth_prob(original._fitness_lth_prob),
    _fitness_del_prob(original._fitness_del_prob),
    _fitness_neu_prob(original._fitness_neu_prob),
    _fitness_adv_prob(original._fitness_adv_prob)
    {
        _rnd_seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        _rnd_gen.seed(_rnd_seed);
    }
    
    std::vector<std::vector<C>> SamplePrior( std::size_t num_samples, std::size_t items_in_sample )
    {
        if ( items_in_sample != _min_vector.size() || items_in_sample != _min_vector.size()  ) {
            throw "Number of items to be sampled is not compatible with min/max vectors.";
        }
        
        std::vector<std::vector<C>> return_vector(0);
        
        for ( auto current_sample=0; current_sample<num_samples; ++current_sample ) {
            
            std::vector<C> tmp_vector(0);
            
            for ( auto current_item=0; current_item<items_in_sample; ++current_item ) {
                
                tmp_vector.push_back( SampleFromDistribution(_min_vector[current_item], _max_vector[current_item]) );
            }
            
            return_vector.push_back(tmp_vector);
        }
        
        return return_vector;
    }
    
private:
    C SampleFromDistribution( C min, C max )
    {
        C tmp_val = 0.0;
        
        if ( min == max ) {
            return min;
        }
        
        switch (_distrib_type) {
            case PriorDistributionType::UNIFORM: {
                boost::random::uniform_01<> dist;
                tmp_val = min + (max-min) * dist(_rnd_gen);
                break;
            }
            
            case PriorDistributionType::LOG_UNIFORM: {
                boost::random::uniform_01<> dist;
                auto power = min + (max-min) * dist(_rnd_gen);
                tmp_val = std::pow(10.0, power);
                break;
            }
                
            case PriorDistributionType::FITNESS_COMPOSITE: {
                
                
                
                // 2017-02-07 changed the parameters to improve the distribution visually
                // define distributions which compose the desired distribution
                boost::random::gamma_distribution<float> dist_gamma(1.0f, 0.1f); // dist_gamma(0.5f, 1.0f);
                boost::random::lognormal_distribution<float> dist_lognorm(0.0f, 2.0f); // (0.092, 6.0f)
                boost::random::normal_distribution<float> dist_norm(1.0f,0.2f);
                boost::random::discrete_distribution<int> dist_disc( { _fitness_lth_prob, _fitness_del_prob, _fitness_neu_prob, _fitness_adv_prob});
                
                /*
                std::cout << " prior distrib. probs: " << _fitness_lth_prob
                << "  " << _fitness_del_prob
                << "  " << _fitness_neu_prob
                << "  " << _fitness_adv_prob
                << "  " << std::endl;
                */
                
                // 2017-02-07
                // may cause some bias if the loop is outside so create small loops for each case
                //do {
                    // first step - what category
                    auto chosen_category = dist_disc(_rnd_gen);
                    
                    switch (chosen_category) {
                            
                        case FITNESS_ADV_IDX:
                            do {
                                tmp_val = FITNESS_NEU_VAL + dist_gamma(_rnd_gen);
                            } while ( tmp_val > max || tmp_val < min );
                            break;
                            
                        case FITNESSS_NEU_IDX:
                            tmp_val = FITNESS_NEU_VAL;
                            // tmp_fitness = 0.0 - std::fabs(dist_norm(_rnd_engine));
                            break;
                            
                        case FITNESS_DEL_IDX:
                            do {
                                tmp_val = 1.0f - dist_lognorm(_rnd_gen);// * dist_norm(_rnd_gen);
                            } while ( tmp_val > max || tmp_val < min );
                            break;
                            
                        case FITNESS_LTH_IDX:
                            tmp_val = 0.0f;
                            break;
                            
                        default:
                            break;
                    } // switch
                    
                    
                //} while ( tmp_val > max || tmp_val < min );

                break;
            }
        }
        
        return tmp_val;
    }
    
    
};

#endif /* PriorSampler_hpp */
