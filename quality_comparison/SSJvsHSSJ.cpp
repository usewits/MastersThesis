#include <iostream>
#include <cmath>
#include <assert.h>
#include <algorithm>
#include <numeric>
#include "sampleJoins.h"

using namespace std;

vector<double> get_distribution(int n, double skew, double ratio = 0.0, double n_discrete = 0.0) {
    if(ratio == 0.0) {
        assert(n_discrete != 0);
        ratio = n_discrete;
    }
    vector<double> w(n);
    for(int i=0; i<n; i++)
        w[i] = pow(mtwist_drand(mt), skew);
    double max_w = *max_element(w.begin(), w.end());
    double sum_w = 0;
    for(int i=0; i<n; i++) {
        w[i] = (w[i]/max_w)*(ratio-1.0)+1.0;
        sum_w += w[i];
    }
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

double aggregate_f(double A, double B, double C) {
    return A+B*C;
}


int main() {
    mt = mtwist_new();
    mtwist_seed(mt, 832982837UL);

    int m = 100;
	double k_factor = 1.0;
	double sigma = 0.99;
	double sigma_factor = 1.0/log(1/sigma);

    int n1          = 100;
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
    auto stratR2 = stratify(R2);
    auto R2stratcounts = stratify_counts(R2);
    
    auto J = join(stratR1, stratR2);
    double true_aggregate = 0.0;
    for(auto j : J)
        true_aggregate += aggregate_f(get<0>(j), get<1>(j), get<2>(j));
    cout << "true  :" << true_aggregate << endl;

    
    //Compute the sampling weights required for SSJ
    vector<double> SSJ_prob(n1);
    for(int i=0; i<n1; i++) {
        double key = R1[i].first;
        SSJ_prob[i] = R2stratcounts[key];
           //note R2stratcounts[key] = m_2(t_1.A)
    }
    vector<double> SSJ_c_prob = get_cdf(SSJ_prob);
    
    int nruns = 1;
    for(int run_i=0; run_i<nruns; run_i++) {
	//SSJ
        vector<pdd> SSJ_S = weighted_sample(R1, SSJ_c_prob, m);
		vector<tdd> SSJ_samp = minijoin(SSJ_S, stratR2);
		double SSJ_estimate = 0.0;
		for(auto t : SSJ_samp)
			SSJ_estimate += aggregate_f(get<0>(t), get<1>(t), get<2>(t));
		SSJ_estimate *= J.size()/(double)SSJ_samp.size();
		cout << "SSJ   :" << SSJ_estimate << endl;

	//HSSJ
		int k = k_factor * m*m * sigma_factor
				* (*max_element(SSJ_prob.begin(), SSJ_prob.end()))
				/ (*min_element(SSJ_prob.begin(), SSJ_prob.end()));//HWS-heuristic
        vector<int> HSSJ_U_indices = sample_indices(R1.size(), k);
		vector<double> HSSJ_U_prob(k);
		vector<pdd> HSSJ_U(k);
		for(int i=0; i<k; i++) {
			HSSJ_U_prob[i] = SSJ_prob[HSSJ_U_indices[i]];
			HSSJ_U[i] = R1[HSSJ_U_indices[i]];
		}
        vector<pdd> HSSJ_S = weighted_sample(HSSJ_U, get_cdf(HSSJ_U_prob), m);
		vector<tdd> HSSJ_samp = minijoin(HSSJ_S, stratR2);
		double HSSJ_estimate = 0.0;
		for(auto t : HSSJ_samp)
			HSSJ_estimate += aggregate_f(get<0>(t), get<1>(t), get<2>(t));
		HSSJ_estimate *= J.size()/(double)HSSJ_samp.size();
		cout << "HSSJ  :" << HSSJ_estimate << " (k = " << k << ")" << endl;
    }
    
    return 0;
}
