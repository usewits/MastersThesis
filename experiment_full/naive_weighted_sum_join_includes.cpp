#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <stdlib.h>

#define pdd pair<double, double>
#define Tstrat map<double, vector<pdd> >
using namespace std;

#define MTWIST_N 624
#define MTWIST_M 397

/* Mersenne Twister library state */
struct mtwist_s {
    /* MT buffer holding N 32 bit unsigned integers */
    unsigned int state[MTWIST_N];

    /* Pointer into above - next int to use */
    unsigned int* next;

    /* Number of remaining integers in state before an update is needed */
    unsigned int remaining;

    /* 1 if a seed was given */
    unsigned int seeded : 1;
};

/* Mersenne Twister state */
typedef struct mtwist_s mtwist;



#define MTWIST_UPPER_MASK 0x80000000UL
#define MTWIST_LOWER_MASK 0x7FFFFFFFUL
#define MTWIST_FULL_MASK 0xFFFFFFFFUL
#define MTWIST_MATRIX_A 0x9908B0DFUL

#define MTWIST_MIXBITS(u, v) (((u)&MTWIST_UPPER_MASK) | ((v)&MTWIST_LOWER_MASK))
#define MTWIST_TWIST(u, v)       \
    ((MTWIST_MIXBITS(u, v) >> 1) ^ \
     ((v)&1UL ? MTWIST_MATRIX_A : 0UL))

/**
 * mtwist_new:
 *
 * Construct a Mersenne Twister object
 *
 * Return value: new MT object or NULL on failure
 */
mtwist* mtwist_new(void) {
    mtwist* mt;

    mt = (mtwist*)calloc(1, sizeof(*mt));
    if (!mt) return NULL;

    mt->remaining = 0;
    mt->next = NULL;
    mt->seeded = 0;

    return mt;
}

/**
 * mtwist_free:
 * @mt: mt object
 *
 * Destroy a Mersenne Twister object
 */
void mtwist_free(mtwist* mt) {
    if (mt) free(mt);
}

/**
 * mtwist_seed:
 * @mt: mt object
 * @seed: seed (lower 32 bits used)
 *
 * Initialise a Mersenne Twister with an unsigned 32 bit int seed
 */
void mtwist_seed(mtwist* mt, unsigned int seed) {
    int i;

    if (!mt) return;

    mt->state[0] = (unsigned int)(seed & MTWIST_FULL_MASK);
    for (i = 1; i < MTWIST_N; i++) {
        mt->state[i] =
            (1812433253UL * (mt->state[i - 1] ^ (mt->state[i - 1] >> 30)) +
             i);
        mt->state[i] &= MTWIST_FULL_MASK;
    }

    mt->remaining = 0;
    mt->next = NULL;

    mt->seeded = 1;
}

static void mtwist_update_state(mtwist* mt) {
    int count;
    unsigned int* p = mt->state;

    for (count = (MTWIST_N - MTWIST_M + 1); --count; p++)
        *p = p[MTWIST_M] ^ MTWIST_TWIST(p[0], p[1]);

    for (count = MTWIST_M; --count; p++)
        *p = p[MTWIST_M - MTWIST_N] ^ MTWIST_TWIST(p[0], p[1]);

    *p = p[MTWIST_M - MTWIST_N] ^ MTWIST_TWIST(p[0], mt->state[0]);

    mt->remaining = MTWIST_N;
    mt->next = mt->state;
}

/**
 * mtwist_u32rand:
 * @mt: mt object
 *
 * Get a random unsigned 32 bit integer from the random number generator
 *
 * Return value: unsigned int with 32 valid bits
 */
unsigned int mtwist_u32rand(mtwist* mt) {
    unsigned int r;

    if (!mt) return 0UL;

    if (!mt->seeded) mtwist_seed(mt, 0);

    if (!mt->remaining) mtwist_update_state(mt);

    r = *mt->next++;
    mt->remaining--;

    /* Tempering */
    r ^= (r >> 11);
    r ^= (r << 7) & 0x9D2C5680UL;
    r ^= (r << 15) & 0xEFC60000UL;
    r ^= (r >> 18);

    r &= MTWIST_FULL_MASK;

    return r;
}

/**
 * mtwist_drand:
 * @mt: mt object
 *
 * Get a random double from the random number generator
 *
 * Return value: random double in the range 0.0 inclusive to 1.0 exclusive;
 *[0.0, 1.0) */
double mtwist_drand(mtwist* mt) {
    unsigned int r;
    double d;

    if (!mt) return 0.0;

    r = mtwist_u32rand(mt);

    d = r / 4294967296.0; /* 2^32 */

    return d;
}


/**
 * mtwist_uniform_int:
 * @a, b; two integers such that a<=b
 *
 * Get an int in an interval uniform randomly from the
 * random number generator.
 *
 * Return value: random interval in range a inclusive to b inclusive;
 * [a,b] 
 */
int mtwist_uniform_int(mtwist* mt, int a, int b) {
    if(b < a) {//invalid range!
        return 0;
    }
    unsigned int range = b-a+1;
    unsigned int scale = 4294967295UL/range;
        //4294967295UL=2^32-1=RAND_MAX for this Mersenne Twister
    unsigned int max_x = range*scale;
    //x will be uniform in [0, max_x[
    //Since past%range=0, x%range will be uniform in [0,range[
    unsigned int x; 
    do {
        x = mtwist_u32rand(mt);
    } while(x >= max_x);

    return a+(x/scale);
    //x is uniform in [0,max_x[ = [0,range*scale[
    //hence x/scale is uniform in [0,range[=[0,b-a+1[
    //thus a+(x/scale) is uniform in [a,b]
    
    //alternative: return a+(x%range); 
    //x is uniform in [0,max_x[ = [0,range*scale[
    //hence (x%range) is uniform in [0,range[=[0,b-a+1[
    //thus a+(x%range) is uniform in [a,b]
}

//The only global variable
mtwist* mt;

enum C_WEIGHT_CONTRIBUTION {C_UNIFORM=1, C_LINEAR=2, C_INVERSE=3};
enum AGGREGATION_APPROACH {UNIFORM_AGGREGATION=1, KEY_STRATA_AGGREGATION=2, RANGE_STRATA_AGGREGATION=3, WEIGHTED_AGGREGATION=4 };
enum SAMPLE_JOIN_APPROACH {U_SAMPLE_JOIN=1, WS_SAMPLE_JOIN=2, AWS_SAMPLE_JOIN=3};
enum AGGREGATION_COLUMN {COL_A=1, COL_B=2, COL_C=3, COL_X=4};


vector<double> get_cdf(const vector<double>& w) {
    vector<double> result(w.size());
	result[0]=w[0];
    for(int i=1; i<w.size(); i++) //set result to cumsum of w
		result[i] = result[i-1]+w[i];
    for(int i=0; i<w.size(); i++) //normalize result (last element should be 1)
        result[i]/=result[w.size()-1];
    return result;
}

Tstrat stratify(const vector<pdd>& R) {//Stratifies R by first argument
	Tstrat result;
	for(vector<pdd>::const_iterator Rit = R.begin(); Rit!=R.end(); ++Rit) {
		if(result.find(Rit->first) == result.end()) {
			result.insert(make_pair(Rit->first,vector<pdd>()));//for first element in strata, insert empty vector
		}
		result[Rit->first].push_back(*Rit);
	}
	return result;
}

map<double, int> stratify_counts(const vector<pdd>& R) {//Stratifies R by first argument
	Tstrat Rstrat = stratify(R);
	
	map<double, int> result;
	for(Tstrat::const_iterator Rit = Rstrat.begin(); Rit!=Rstrat.end(); ++Rit) {
		result[Rit->first]=(Rit->second).size();
	}
	return result;
}

map<double, int> strat_count_multiply(const map<double, int>& a, const map<double, int>& b) {
	map<double, int> result;
	for(map<double, int>::const_iterator ait = a.begin(); ait != a.end(); ++ait) {
		map<double, int>::const_iterator bit = b.find(ait->first);
		if(bit != b.end()) {
			result.insert(make_pair(ait->first, ait->second*bit->second));
		}
	}
	return result;
}

template <typename T> vector<T> sample(const vector<T>& R, int k) {//Obtain sample with replacement of size k
	vector<T> result(k);
	int n = R.size();
	for(int i=0; i<k; i++) {
		result[i]=R[mtwist_uniform_int(mt,0,n-1)];
	}
	return result;
}

vector<int> sample_indices(int n, int k) {//Obtain sample with replacement of size k from a list of indices 0:n-1
	vector<int> result(k);
	for(int i=0; i<k; i++) {
		result[i]=mtwist_uniform_int(mt,0,n-1);
	}
	return result;
}

template <typename T> vector<T> weighted_sample(const vector<T>& R, const vector<double>& c_p, int m) {
	//Obtain weighted sample with replacement of size m
	//NOTE: we assume c_p is cumilative probability distribution from 0 to 1!
	vector<T> result(m);
	int n = R.size();
	for(int i=0; i<m; i++) {
		double random_variate = mtwist_drand(mt);
		vector<double>::const_iterator c_p_it = upper_bound(c_p.begin(), c_p.end(), random_variate);
		int index = c_p_it-c_p.begin();
		result[i]=R[index];
	}
	return result;
}

vector<int> weighted_sample_indices(int n, const vector<double>& c_p, int m) {
	vector<int> result(m);
	for(int i=0; i<m; i++) {
		double random_variate = mtwist_drand(mt);
		vector<double>::const_iterator c_p_it = upper_bound(c_p.begin(), c_p.end(), random_variate);
		int index = c_p_it-c_p.begin();
		result[i]=index;
	}
	return result;
}


template <typename T> vector<T> approximate_weighted_sample(const vector<T>& R, const vector<double>& c_w, int m, double w_ratio, double delta) {
	//NOTE: we assume c_w is cumilative weight distribution from 0 to 1!
	//delta is the desired maximum normalised absolute difference between the projected total weight from U and the total weight in R
	int min_k = ceil(w_ratio*m*m);

	if(min_k > R.size()) {//NOTE: in practice we may want to switch from AWS to WS at much lower sampling fraction, depending on runtime cutoff.
						  //This switch is just here to avoid run-away runtimes, and will not influence determination of crossover since it is much beyond the crossover point
		printf("WARNING: oversampling in AWS! Early switch to regular sampling.\n");
		printf("AWS effective sampling fraction 1.0\n");
		return(weighted_sample(R, c_w, m));
	}

	double memory_factor = 2.0;
	double U_weight = 0;

	vector<double> w_U;
	vector<int> U_indices;
	w_U.reserve(ceil(memory_factor*min_k));
	U_indices.reserve(ceil(memory_factor*min_k));

	int n=R.size();
	int k;
	for(k=0; k<min_k || abs(U_weight*n/(double)k-1.0) > delta; k++) {
						//Note; U_weight*n/(double)k is the projected total weight from U
		//int j = *sample_indices(n,1).begin();
		int j=mtwist_uniform_int(mt,0,n-1);
		U_indices.push_back(j);

		if(j == 0)
			w_U.push_back(c_w[j]);
		else
			w_U.push_back(c_w[j] - c_w[j-1]);

		U_weight += w_U[k];
		if(k > R.size()) {
			printf("WARNING: oversampling in AWS! Late switch to regular sampling.\n");
			printf("AWS effective sampling fraction 1.0\n");
			return(weighted_sample(R, c_w, m));
		}
	}
	if(k > min_k) {//Only print effective sampling fraction if it is not min_k/n1
		cout << "AWS effective sampling fraction " << k/(double)R.size() << " (default sampling fraction is " << min_k/(double)R.size() << ")" << endl;
	}

	vector<double> c_w_U = get_cdf(w_U);
	
	vector<int> S_indices = weighted_sample(U_indices, c_w_U, m);
	vector<T> S(m);
	for(int i=0; i<m; i++) {
		S[i] = R[S_indices[i]];
	}
	return S;
}
	


double avg(const vector<pdd>& R, int column) {
	double result =0;
	if(column == 0) {
		for(vector<pdd>::const_iterator Rit = R.begin(); Rit != R.end(); ++Rit) {
			result += Rit->first;
		}
	}
	if(column == 1) {
		for(vector<pdd>::const_iterator Rit = R.begin(); Rit != R.end(); ++Rit) {
			result += Rit->second;
		}
	}
	return result/R.size();
}

double X_function(double el_A, double el_B, double el_C) {
	return el_C*el_B+el_A; //Insert custom aggregation function here
}


