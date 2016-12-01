Rcpp::NumericVector R1A(R_R1A);
Rcpp::NumericVector R1B(R_R1B);
Rcpp::NumericVector R2A(R_R2A);
Rcpp::NumericVector R2C(R_R2C);

Rcpp::NumericVector Rprob(R_prob);//The weights of R1 
					//(Rprob ``emulates a function'', so should depend only on attribute values and/or index)

int n1 = R1A.size();
int n2 = R2A.size();

int m = as<int>(R_m);
double delta = as<double>(R_delta);
int n_replications = as<int>(R_n_replications);

int C_weight_contribution = as<int>(R_C_weight_contribution);
//enum C_WEIGHT_CONTRIBUTION {C_UNIFORM=1, C_LINEAR=2, C_INVERSE=3};

int aggregation_index = as<int>(R_aggregation_index);
//enum AGGREGATION_COLUMN {COL_A=1, COL_B=2, COL_C=3, COL_X=4};

int aggregation_approach = as<int>(R_aggregation_approach);
//enum AGGREGATION_APPROACH {UNIFORM_AGGREGATION=1, KEY_STRATA_AGGREGATION=2, RANGE_STRATA_AGGREGATION=3, WEIGHTED_AGGREGATION=4 };

int sampling_approach = as<int>(R_sampling_approach);
//enum SAMPLE_JOIN_APPROACH {U_SAMPLE_JOIN=1, WS_SAMPLE_JOIN=2, AWS_SAMPLE_JOIN=3};

int use_stratum_correction = as<int>(R_use_stratum_correction);


mt = mtwist_new();
mtwist_seed(mt, 832982837UL);

//Convert columns to combined vectors
vector<pdd> R1(n1);
vector<pdd> R2(n2);
vector<int> R1_indices(n1);

for(int i=0; i<n1; i++) {
	R1[i]=make_pair(R1A[i],R1B[i]);
	R1_indices[i]=i;
}

for(int i=0; i<n2; i++) {
	R2[i]=make_pair(R2A[i],R2C[i]);
}

map<double, int> R1stratcounts = stratify_counts(R1);
map<double, int> R2stratcounts = stratify_counts(R2);
map<double, int> Jstratcounts = strat_count_multiply(R1stratcounts,R2stratcounts);

int join_size = 0;
for(map<double, int>::iterator Jstratcount_it=Jstratcounts.begin(); Jstratcount_it!=Jstratcounts.end(); ++Jstratcount_it) {
	join_size += Jstratcount_it->second;
}

Tstrat R2strat = stratify(R2);

map<double, vector<double> > R2weights;
map<double, vector<double> > R2c_weights;


if(C_weight_contribution != C_UNIFORM && aggregation_approach == UNIFORM_AGGREGATION) {
	cout << "WARNING: OVERWRITING C_weight_contribution SO SAMPLE-JOIN OUTPUT IS (APPROXIMATELY) UNIFORM" << endl;
	C_weight_contribution=C_UNIFORM;
}

//Make sure all keys in R2 can be found in R1stratcounts
for(Tstrat::iterator R2it = R2strat.begin(); R2it != R2strat.end(); ++R2it) {
	double key = R2it->first;
	if(R1stratcounts.find(key) == R1stratcounts.end()) {
		R1stratcounts[key] = 0;
	}
	vector<pdd> &stratum = R2it->second;

	//Create distributions within strata based on C
	double norm = 0;
	R2weights.insert(make_pair(key,vector<double>(stratum.size())));
	for(int i=0; i<stratum.size(); i++) {
		double C = stratum[i].second;
		if(C_weight_contribution == C_UNIFORM) {
			R2weights[key][i] = 1;
		} else if(C_weight_contribution == C_LINEAR) {
			R2weights[key][i] = C;
		} else if(C_weight_contribution == C_INVERSE) {
			R2weights[key][i] = 1.0/C;
		}

		norm+=R2weights[key][i];
	}
	for(int i=0; i<stratum.size(); i++) {
		R2weights[key][i]/=norm;
	}
	R2c_weights.insert(make_pair(key,get_cdf(R2weights[key])));
}

if(Rprob.size() != n1) {
	cout << "ERROR: prob IS NOT THE SAME LENGTH AS R1" << endl;
	return NULL;
}

vector<double> prob(n1);

if(aggregation_approach == UNIFORM_AGGREGATION) {
	//In this case we require the sample to be (approximately) uniform	
	if(sampling_approach == U_SAMPLE_JOIN) {
		cout << "ERROR: CANNOT PERFORM UNIFORM AGGREGATION OVER (non-uniform) U-JOIN" << endl;
		return NULL;
	}
	cout << "WARNING: OVERWRITING prob SO SAMPLE-JOIN OUTPUT IS (APPROXIMATELY) UNIFORM" << endl;
	use_stratum_correction=true;
	for(int i=0; i<n1; i++) {
		Rprob[i]=1;
	}
}

//For now, do not remove elements that do not join to make the effect of this visible

////Remove elements that will not join
//for(int i=0; i<n1; i++) {
//	double key = R1[i].first;
//	if(R2stratcounts[key] == 0) {
//		Rprob[i]=0;
//	}
//}

//Correct for strata if requisted by user
if(use_stratum_correction) {
	for(int i=0; i<n1; i++) {
		double key = R1[i].first;
		prob[i] = R2stratcounts[key]*Rprob[i];
				//R2stratcounts[key] = m_2(t_1.A)
	}
}

double normalisation = 0;
//Normalize prob
for(int i=0; i<n1; i++) {
//	double key = R1[i].first;
//	if(R2stratcounts[key] == 0) continue;
	normalisation+=prob[i];
}

double prob_max = 0;
double prob_min = 1;
for(int i=0; i<n1; i++) {
	prob[i]/=normalisation;

	double key = R1[i].first;
	if(R2stratcounts[key] == 0) continue;

	if(prob[i] > prob_max)
		prob_max = prob[i];
	if(prob[i] < prob_min)
		prob_min = prob[i];
}

double w_ratio = prob_max / prob_min;

if(w_ratio*m*m/(double)n1 > 1.0) {
	cout << "ERROR: sampling fraction is too high for AWS" << endl;
	vector<int> result(n_replications);
	for(int i=0; i<n_replications; i++)
		result[i]=0;
	return Rcpp::wrap(result);
}

if(sampling_approach == U_SAMPLE_JOIN) {
	cout << "WARNING: OVERWRITING prob TO REPRESENT U_SAMPLE_JOIN" << endl;
	normalisation = 0;
	for(int i=0; i<n1; i++) {
		double key = R1[i].first;
		prob[i] = 1;
		normalisation+=prob[i];
	}
	for(int i=0; i<n1; i++) {
		prob[i]/=normalisation;
	}
}

//Obtain cumilative probabilities
vector<double> c_prob = get_cdf(prob);


double output_probability_normalization=-1;//Normalization of the output probabilities (norm of h)
               //Not to be confused with normalization of prob, the selection probabilies in R1 (norm of w)
               //In the case where we use a uniform minijoin, they are the same thing
if(aggregation_approach == WEIGHTED_AGGREGATION) {
    //We obtain this in linear time here, in practice we assume this to be precomputed
    output_probability_normalization=0;
    for(int i=0; i<n1; i++) {
        double key = R1[i].first;
		if(R2stratcounts[key] == 0)
			continue;
		vector<double> &R2_stratum_weights = R2weights[key];
		for(int j=0; j< R2_stratum_weights.size(); j++) {
        	output_probability_normalization+=prob[i]*R2_stratum_weights[j];
		}
            //for uniform minijoin; divided by R2stratcounts[key], times R2stratcounts[key] for size of stratum
    }
}

Rcpp::NumericVector results(n_replications);

double average_m = 0;

//Repeatedly sample and estimate aggregation
for(int iteration=0; iteration<n_replications; iteration++) {
	double result = 0;
	
	//obtain sample to join
	vector<pdd> S(m);
	vector<int> S_indices;

	if(sampling_approach == U_SAMPLE_JOIN) {
		//naive weighted approach:
		S_indices = sample(R1_indices,m);
	} else if(sampling_approach == WS_SAMPLE_JOIN) {
		//custom weighted approach (exact):
		S_indices = weighted_sample(R1_indices, c_prob, m);
	} else if(sampling_approach == AWS_SAMPLE_JOIN) {
		//custom weighted approach (approximate):
		S_indices = approximate_weighted_sample(R1_indices, c_prob, m, w_ratio, delta);
	}

	//Obtain sample values from indices
	for(int i=0; i<m; i++) {
		S[i] = R1[S_indices[i]];
	}
	
	//Correct for number of elements in S that can actually be joined
	//Alternatives:
		//-Set prob[i] to zero for these values
		//-Avoid empty strata
		//-Remove values from S
	int effective_m = 0;
	for(int i=0; i<m; i++) {
		double key = R1[S_indices[i]].first;
		if(R2stratcounts[key] > 0) {
			effective_m++;
		}
	}
	average_m += effective_m/(double)n_replications;

	if(aggregation_approach == WEIGHTED_AGGREGATION) {
		for(int i = 0; i < m; i++) {
			int S_index = S_indices[i];

			double key = R1[S_index].first;
            if(R2strat.find(key) == R2strat.end()) {
                continue;//The element in S does not join!
            }
			//Obtain uniform random element in R2 that joins.
			int R2element_index = *weighted_sample_indices(R2strat[key].size(),R2c_weights[key], 1).begin();
			pdd R2element = R2strat[key][R2element_index];

			double X_term = X_function(R1[S_index].first, R1[S_index].second, R2element.second);
			
			double prob_i=(prob[S_index]*R2weights[key][R2element_index])/output_probability_normalization;
					//this is h(S_i) (not to be confused with w(S_i)=prob[S_index]
	
            if(aggregation_index == COL_A) {
			    result += key/prob_i;
            }
            if(aggregation_index == COL_B) {
			    result += R1[S_index].second/prob_i;
            }
            if(aggregation_index == COL_C) {
			    result += R2element.second/prob_i;
            }
            if(aggregation_index == COL_X) {
			    result += X_term/prob_i;
            }
			//result = 1/|S| * sum_i (s_i/p_i)
			//NOTE: we use that prob is normalized here!
		}
        result /= (double)(effective_m);
	} else {//For UNIFORM_AGGREGATION, KEY_STRATA_AGGREGATION and RANGE_STRATA_AGGREGATION we use the stratification of S
		Tstrat Sstrat = stratify(S);

		//Loop over strata (stratified by A, the join attribute)
		for(Tstrat::iterator Sstrat_it = Sstrat.begin(); Sstrat_it != Sstrat.end(); ++Sstrat_it) {
			double key = Sstrat_it->first;
			vector<pdd>& Sstrat = Sstrat_it->second;

			//Stratum does not persist after join
			if(Jstratcounts.find(key) == Jstratcounts.end()) {
				continue;
			}

			if(aggregation_approach == UNIFORM_AGGREGATION) {
				if(aggregation_index == COL_A) {
					result += key * Sstrat.size() * join_size/(double)(effective_m);
					//Note that key = avg(Sstrat, 0)
					//avg * Sstrat.size is the contribution to the sum over S. join_size/m is |J|/|S|
				}
				if(aggregation_index == COL_B) {
					result += avg(Sstrat, 1) * Sstrat.size() * join_size/(double)(effective_m);
					//avg * Sstrat.size is the contribution to the sum over S. join_size/m is |J|/|S|
				}
				if(aggregation_index == COL_C) {
					int Sstrat_size = Sstrat.size();
					for(int i=0; i<Sstrat_size; i++) {
						pdd R2element = *weighted_sample(R2strat[key], R2c_weights[key], 1).begin();
						result += R2element.second * join_size/(double)(effective_m);
						//R2element.second is the C column of one joined element. join_size/m is |J|/|S|
					}
				}
				if(aggregation_index == COL_X) { //``Generic'' aggregation
					for(vector<pdd>::iterator element_it = Sstrat.begin();
											  element_it != Sstrat.end();
											++element_it) {
						//Obtain uniform random element in R2 that joins.
						pdd R2element = *weighted_sample(R2strat[key], R2c_weights[key], 1).begin();

						double X_term = X_function(element_it->first, element_it->second, R2element.second);

						result += X_term * join_size/(double)(effective_m);
					}
				}
			}


			if(aggregation_approach == KEY_STRATA_AGGREGATION) {
				if(aggregation_index == COL_A) {
					result += key * Jstratcounts[key];
					//Note that key = avg(Sstrat, 0)
				}
				if(aggregation_index == COL_B) {
					result += avg(Sstrat, 1) * Jstratcounts[key];
				}
				if(aggregation_index == COL_C) {
					double est_sum_C = 0;
					int Sstrat_size = Sstrat.size();
					//We use the Hansen Hurwitz estimator within a stratum, since C_weight_contribution may be weighted
					//Notice that the Hansen Hurwitz estimator reduces to the average in the uniform case
					for(int i=0; i<Sstrat_size; i++) {
						int R2element_index = *weighted_sample_indices(R2strat[key].size(), R2c_weights[key], 1).begin();
						pdd R2element = R2strat[key][R2element_index];
						est_sum_C += R2element.second/(R2weights[key].size() * R2weights[key][R2element_index]);
					}
					result += est_sum_C/(double)Sstrat_size * Jstratcounts[key];
				}
				if(aggregation_index == COL_X) { //``Generic'' aggregation
					double est_sum_X = 0;
					int Sstrat_size = Sstrat.size();
					//We use the Hansen Hurwitz estimator within a stratum, since C_weight_contribution may be weighted
					//Notice that the Hansen Hurwitz estimator reduces to the average in the uniform case
					for(vector<pdd>::iterator element_it = Sstrat.begin();
											  element_it != Sstrat.end();
											++element_it) {

						int R2element_index = *weighted_sample_indices(R2strat[key].size(), R2c_weights[key], 1).begin();
						pdd R2element = R2strat[key][R2element_index];

						double X_term = X_function(element_it->first, element_it->second, R2element.second);

						est_sum_X += X_term/(R2weights[key].size() * R2weights[key][R2element_index]);
					}
					result += est_sum_X/(double)Sstrat_size * Jstratcounts[key];
				}
			}

            if(aggregation_approach == RANGE_STRATA_AGGREGATION) {
                //TODO

            }
		}
	}

	results[iteration]=result;
}

if(average_m < 0.8*m)
	printf("WARNING: ");
if(average_m < 0.999*m)
	printf("On average, joined %.1f out of %d tuples\n", average_m,m);



return Rcpp::wrap(results);

