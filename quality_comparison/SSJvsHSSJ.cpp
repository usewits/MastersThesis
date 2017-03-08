#include <iostream>
#include <cmath>
#include <cstring>
#include <assert.h>
#include <algorithm>
#include <numeric>
#include <set>
#include <functional>
#include "sampleJoins.h"

using namespace std;

//Generate weights (n elements), with a selected skew ratio and number of discrete values.
vector<double> get_distribution(int n, double skew, double ratio = 0.0, double n_discrete = 0.0) {
    if(ratio == 0.0) {
        //Either the ratio or n_discrete has to be defined 
        assert(n_discrete != 0);
        ratio = n_discrete;
    }
    vector<double> w(n);
    for(int i=0; i<n; i++)
        w[i] = pow(mtwist_drand(mt), skew); //the weights are in [0,1[ with a (polynomial) skew
    double max_w = *max_element(w.begin(), w.end());
    double sum_w = 0;
    for(int i=0; i<n; i++) {
        w[i] = (w[i]/max_w)*(ratio-1.0)+1.0; //the weights are in [1, ratio[
        sum_w += w[i];
    }
    //normalise the weights
    for(int i=0; i<n; i++)
        w[i] /= sum_w;

    if(n_discrete > 0) {
        max_w = *max_element(w.begin(), w.end());
        for(int i=0; i<n; i++) {
            w[i] *= n_discrete/max_w;
            w[i] = round(w[i]);
        }
    }
    return w;    
}

//Generic function to estimate aggregates over joins
//It can be used to obtain SSJ, HSSJ, WS-Join, HWS-Join or US-Join estimates (both filtered and unfiltered)
//Runs in O(n_1 log n_1+n_2 log n_2) (not as fast as possible, in favour of shorter code)
//Output probability is h1*h2
double generic_sample_join(function<double(double,double)> h1, function<double(double)> h2, int m,
                                const vector<pdd>& R1, const vector<pdd>& R2,
                                function<vector<int>(int, vector<double>)> range_sampler,
                                            //sample(sample_size, weights)
                                function<double(double, double, double)> aggregation_f,
                                function<bool(double, double)> R1_filter,//Ri_filter are predicates; true => selected
                                function<bool(double, double)> R2_filter, bool filtered_estimator) {

    //Compute (filtered) stratum weights and cdfs
    Tstrat R2_stratified = stratify(R2);
    map<double, double> R2_stratum_weights;
    map<double, double> R2_filtered_stratum_weights;
    map<double, vector<double> > R2_stratum_cdfs;
    for(auto stratum : R2_stratified) {
        double key = stratum.first;
        double norm = 0.0;
        double filtered_norm = 0.0;
        vector<double> stratum_weights(stratum.second.size());
        for(int i=0; i<stratum.second.size(); i++) {
            pdd t2 = stratum.second[i];
            stratum_weights[i] = h2(t2.second);
            norm += stratum_weights[i];
            if(R2_filter(t2.first, t2.second))
                filtered_norm += stratum_weights[i]; 
        }
        R2_stratum_weights[key] = norm;
        R2_filtered_stratum_weights[key] = norm;
        R2_stratum_cdfs[key] = get_cdf(stratum_weights);
    }

    //Compute normalisation factors
    double normalisation = 0.0;          //Total weight of all elements in J
    double filtered_normalisation = 0.0; //Total weight of selection sigma(J)
    vector<double> R1_sample_weights(R1.size());                //Sampling weights in R1
    vector<double> R1_filtered_sample_weights(R1.size(), 0.0);  //Filtered sampling weights
    for(int i=0; i<R1.size(); i++) {
        pdd t1 = R1[i];
        R1_sample_weights[i] = h1(t1.first, t1.second) * R2_stratum_weights[t1.first];
        normalisation += R1_sample_weights[i];
        if(R1_filter(t1.first, t1.second)) {
            R1_filtered_sample_weights[i] = h1(t1.first, t1.second) * R2_filtered_stratum_weights[t1.first];
            filtered_normalisation += R1_filtered_sample_weights[i];
        }
    }
    
    vector<int> S_indices = range_sampler(m, R1_sample_weights);
    vector<pdd> S(m);
    vector<double> S_weights(m);
    for(int i=0; i<m; i++) {
        S[i] = R1[S_indices[i]];
        S_weights[i] = R1_sample_weights[S_indices[i]];
    }
    vector<tdd> sample = minijoin(S, R2_stratified);
    
    double estimate = 0.0;
    
    int filtered_sample_size = 0;
    for(int i=0; i<m; i++) {
        tdd t = sample[i];
        double p_t = h1(get<0>(t), get<1>(t))*h2(get<2>(t));
        if(R1_filter(get<0>(t), get<1>(t)) && R2_filter(get<0>(t), get<2>(t))) {
            estimate += aggregation_f(get<0>(t), get<1>(t), get<2>(t))/p_t;
            filtered_sample_size++;
        }
    }
    
    if(filtered_estimator) {
        estimate *= filtered_normalisation/(double) filtered_sample_size;
    } else {
        estimate *= normalisation/(double) m;
    }
    return estimate;
}


int main() {
    mt = mtwist_new();
    mtwist_seed(mt, 832982837UL);

    int m = 20;
    double k_factor = 1.0;
    double sigma = 0.99;
    double sigma_factor = 1.0/log(1/sigma);

    int n1          = 10000;
    double skew1    = 1.0;
    double ratio1   = 20.0;
    int n_discrete1 = 10.0;
    vector<double> R1A = get_distribution(n1,skew1,ratio1,n_discrete1);
    vector<double> R1B = get_distribution(n1,1.0,n1);
    vector<pdd> R1 = zipvec(R1A, R1B);
    auto stratR1 = stratify(R1);

    int n2          = 1000;
    double skew2    = 1.0;
    double ratio2   = 50.0;
    int n_discrete2 = 10.0;
    vector<double> R2A = get_distribution(n2,skew2,ratio2,n_discrete2);
    vector<double> R2C = get_distribution(n2,1.0,n2);
    vector<pdd> R2 = zipvec(R2A, R2C);
    Tstrat stratR2 = stratify(R2);
   
    auto aggregate_f = [] (double A, double B, double C) -> double {return C;};

    auto h1_unif =     []   (double A, double B) -> double {return 1.0;};
    auto h1_US = [&stratR2] (double A, double B) -> double {return 1.0/(double)(stratR2[A].size());};
    auto h1_weighted = []   (double A, double B) -> double {return 1.0;};

    auto h2_unif = [] (double C) -> double {return 1.0;};
    auto h2_weighted = [] (double C) -> double {return C;};

    auto exact_sampler = [] (int m, vector<double> w) -> vector<int> {
                                return weighted_sample_indices(w.size(), get_cdf(w), m);
                            };
    auto heuristic_sampler = [&sigma,&k_factor] (int m, vector<double> w) -> vector<int> {
                                double w_max = (*max_element(w.begin(), w.end()));
                                double w_min = (*min_element(w.begin(), w.end()));
                                double sigma_factor = 1.0/log(1.0/sigma);

                                int k = k_factor*sigma_factor* m*m * w_max/w_min;
    
                                vector<int> U = sample_indices(w.size(), k);
                                vector<double> U_w(k);
                                for(int i=0; i<k; i++) {
                                    U_w[i] = w[U[i]];
                                }
                                return weighted_sample(U, get_cdf(U_w), m);
                            };

    auto no_filter = [] (double X, double Y) -> bool {return true;};
    auto range_filter = [] (double X, double Y) -> bool {return Y > 1.5e-3;};

    auto R1_filter = no_filter;
    auto R2_filter = range_filter;
 
    auto J = join(stratR1, stratR2);
    
    //Compute the sampling weights required for SSJ
    vector<double> SSJ_prob(n1);
    for(int i=0; i<n1; i++) {
        double key = R1[i].first;
        SSJ_prob[i] = stratR2[key].size();
           //note stratR2[key].size() = m_2(t_1.A)
    }
    vector<double> SSJ_c_prob = get_cdf(SSJ_prob);
    
    set<int> sampling_methods_used = {0,1,2,3,4};
    string                          sample_types[] = {"SSJ     ","HSSJ    ","WS-Join ","HWS-Join",  "US-Join "};
    function<double(double,double)> h1_functions[] = { h1_unif,   h1_unif,  h1_weighted,h1_weighted, h1_US};
    function<double(double)>        h2_functions[] = { h2_unif,   h2_unif,  h2_weighted,h2_weighted, h2_unif};
    function<vector<int>(int, vector<double>)> samplers[] = 
                        { exact_sampler, heuristic_sampler, exact_sampler, heuristic_sampler, exact_sampler};
    
    set<int> filter_methods_used = {0,1,2};
    string filter_types[] = {"full      ","filtered  ","fltr.naive"};
    function<bool(double,double)> R1_filters[] = {no_filter, R1_filter, R1_filter};
    function<bool(double,double)> R2_filters[] = {no_filter, R2_filter, R2_filter};
    bool filtered_estimations[] = {false, true, false};

    vector<double> true_aggregates(3, 0.0);
    vector<double> selectivities(3, 0.0);
    for(int i_f : filter_methods_used) {
        for(auto j : J)
            if(R1_filters[i_f](get<0>(j), get<1>(j)) && R2_filters[i_f](get<0>(j), get<2>(j))) {
                true_aggregates[i_f] += aggregate_f(get<0>(j), get<1>(j), get<2>(j));
                selectivities[i_f] ++;
            }

        selectivities[i_f] /= (double) J.size();
        cout << "Exact aggregation (" << filter_types[i_f] << ") :" << true_aggregates[i_f] << " (selectivity " << selectivities[i_f]*100 << "%)" << endl;
    }
     
    int nruns = 100;
    
    map< pair<int, int>, vector<double> > relative_errors;
    
    for(int i_s : sampling_methods_used)
    for(int i_f : filter_methods_used) {
        relative_errors[make_pair(i_s, i_f)] = vector<double>(nruns, 0);
    }
    
    int k = k_factor * m*m * sigma_factor
            * (*max_element(SSJ_prob.begin(), SSJ_prob.end()))
            / (*min_element(SSJ_prob.begin(), SSJ_prob.end()));//HWS-heuristic
    cout << "k          :" << k << " (should be smaller than " << R1.size() << ")"<< endl;
    if(k > R1.size()) {
        cout << "WARNING: Skipping Heuristic methods!" << endl;
        sampling_methods_used.erase(1);
        sampling_methods_used.erase(3);
    }

    
    for(int run_i=0; run_i<nruns; run_i++) {
        if(round(10*run_i/(double)nruns) > round(10*(run_i-1)/(double)nruns)) {
            cout << round(100*run_i/(double)nruns) << "% done.." << endl;
        }
        for(int i_s : sampling_methods_used)
        for(int i_f : filter_methods_used) {
            double estimate = generic_sample_join(h1_functions[i_s], h2_functions[i_s], 
                                                  round(m/selectivities[i_f]), R1, R2, samplers[i_s], aggregate_f, 
                                                  R1_filters[i_f], R2_filters[i_f], filtered_estimations[i_f]);
            relative_errors[make_pair(i_s, i_f)][run_i] = abs(true_aggregates[i_f]-estimate)/true_aggregates[i_f];
        }
    }
    
    for(int i_f : filter_methods_used)
    for(int i_s : sampling_methods_used) {
        cout << sample_types[i_s] << "(" << filter_types[i_f] << ", sample size " << round(m/selectivities[i_f]) << "):" << endl;
        show_sigma_levels(relative_errors[make_pair(i_s, i_f)]);
    }
    
    return 0;
}
