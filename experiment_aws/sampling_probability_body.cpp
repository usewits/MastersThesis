Rcpp::IntegerVector Ubar(R_Ubar);
Rcpp::IntegerVector S(R_S);
Rcpp::NumericVector prob(R_prob);
int n = prob.length();
int m = S.length();
int k = Ubar.length()+m;

double norm_U=0;
for(int i=0; i<m; i++) {
	norm_U += prob[S[i]-1];
}
for(int i=0; i<k-m; i++) {
	norm_U += prob[Ubar[i]-1];
}

map<int,int> counts;
for(int i=0; i<k-m; i++) {
	counts[Ubar[i]-1]++;//obtain counts in Ubar only
}
double log_n_perms_Ubar = lgamma(k-m+1)-log_prod_fac(counts);

for(int i=0; i<m; i++) {//add counts in S to get total counts S+Ubar=U
	counts[S[i]-1]++;
}
double log_n_perms_U = lgamma(k+1)-log_prod_fac(counts);

double log_p_aws_w_given_u = 0;
for(int i=0; i<m; i++) {
	int s = S[i]-1;
int s_count=0;
map<int,int>::iterator s_count_it = counts.find(s);
if(s_count_it != counts.end()) {
	s_count = (*s_count_it).second;
}
	log_p_aws_w_given_u += log(prob[s]) + log(s_count); //here counts[s] is total count in U
}

counts.clear();//new counts to obtain counts in S only
for(int i=0; i<m; i++) {
	counts[S[i]-1]++;
}
double log_n_perms_S = lgamma(m+1)-log_prod_fac(counts);

log_p_aws_w_given_u += log_n_perms_S-log(norm_U)*m;

double log_prob = log_p_aws_w_given_u-k*log(n)+log_n_perms_U-log_n_perms_Ubar;

return Rcpp::wrap(log_prob);
