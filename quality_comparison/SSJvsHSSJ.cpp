#include <iostream>
#include <cmath>
#include <assert.h>
#include <algorithm>
#include <numeric>
#include <functional>
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

double get_WS_join_estimate(function<double(pdd)> *h1, function<double(double)> *h2, vector<double>* c_prob, int m,
			map<double,double> *R2_norm, const vector<pdd>& R1, const vector<pdd>& R2, const Tstrat& stratR2, const vector<tdd>& J) {
	assert(h1 || !h2);//if h2 is defined, than h1 should be defined as well
										 //(h2 => h1) <=> (!(!h1 && h2)) <=> (h1 || !h2)
	bool is_SSJ = !h1 && !h2;//If no output weights are defined, assume uniform output-probability (SSJ)
	bool c_prob_is_local = (c_prob == NULL);
	bool compute_R2_normalization = (h2 && !R2_norm);
	if(c_prob_is_local) {
		//Compute the sampling weights required for WS-join
		vector<double> prob( R1.size() );
		for(int i=0; i<R1.size(); i++) {
			double key = R1[i].first;
			if(is_SSJ) {
				prob[i] = stratR2.at(key).size();//= m_2(t_1.A)
			} else if(!h2) {
				prob[i] = stratR2.at(key).size() * (*h1)(R1[i]);
			} else {
				prob[i] = 0.0;
				for(pdd t2 : stratR2.at(key)) {
					prob[i] += (*h2)(t2.second);
				}
				prob[i] *= (*h1)(R1[i]);
			}
		}
		c_prob = new vector<double>(get_cdf(prob));
	}
	if(compute_R2_normalization) {
		R2_norm = new map<double,double>();
		for(pair<double, vector<pdd> > R2stratum : stratR2) {
			double key = R2stratum.first;
			double normalisation = 0.0;
			for(pdd t2 : R2stratum.second) {
				normalisation += (*h2)(t2.second);
			}
			(*R2_norm)[key] = normalisation;
		}
	}
    vector<int> S_indices = weighted_sample_indices(R1.size(), *c_prob, m);
	vector<pdd> S(m);
	for(int i=0; i<m; i++)
		S[i] = R1[S_indices[i]];
	vector<tdd> samp = minijoin(S, stratR2);

	//Use Hansen-Hurwitz to estimate the result. This reduces to the uniform estimator in the SSJ case.
	double estimate = 0.0;
	for(int i=0; i<m; i++) {
		tdd t = samp[i];
		int index = S_indices[i];
		double tA = get<0>(t);
		double tB = get<1>(t);
		double tC = get<2>(t);
		if(h2) {
			//double p2 = (*h2)(tC)/((*R2_norm)[tC]);
			estimate += aggregate_f(tA,tB,tC)/((*h1)(make_pair(tA,tB))*(*h2)(tC));
		} else {
			if(h1) {
				estimate += aggregate_f(tA,tB,tC)/(*h1)(make_pair(tA, tB));
			} else {
				double p1 = (*c_prob)[index];
				if(index != 0)
					p1 -= (*c_prob)[index-1];//p1 is the selection probability in R1
				p1 *= stratR2.at(tA).size()/R2.size();
				estimate += aggregate_f(tA,tB,tC)/p1;
			}
		}
	}
	estimate /= (double)samp.size();

	if(c_prob_is_local)
		delete c_prob;
	if(compute_R2_normalization)
		delete R2_norm;

	return estimate;
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
		double SSJ_estimate = get_WS_join_estimate(NULL, NULL, &SSJ_c_prob, m, NULL, R1, R2, stratR2, J);
		cout << "SSJ   :" << SSJ_estimate << endl;

	//WS-join
		function<double(pdd)> h1 = [] (pdd t1) -> double {return t1.first+t1.second;};
		function<double(double)> h2 = [] (double t2C) -> double {return t2C;};
		double WS_join_estimate = get_WS_join_estimate(&h1, &h2, &SSJ_c_prob, m, NULL, R1, R2, stratR2, J);
		cout << "WSjoin:" << WS_join_estimate << endl;

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
