#include <cmath>
#include <map>

using namespace std;

double log_prod_fac(map<int,int>& m) {
	double r=0;
	for(map<int,int>::iterator mi = m.begin(); mi != m.end(); ++mi) {
		r+=lgamma((*mi).second+1);//lgamma(x+1)=log(x!)
	}
	return r;
}

