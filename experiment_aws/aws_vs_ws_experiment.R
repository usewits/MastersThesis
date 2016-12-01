source("solvers.R")
source("sampling_utils.R")
source("distribution_generation.R")
source("sampling_probability.R")

library(doParallel)
library(memoise)

get_relative_errors_given_distribution = function(n, k, m, skew, max_weight_ratio, n_outer, n_inner, n_cores=1) {
    if(skew < 0) {
        print("ERROR: invalid parameters");
        print("invalid skew!");
        return(NULL);
    }
    if(n<k || k<m) {
        print("ERROR: invalid parameters");
        print("we need n>=k>=m!");
        return(NULL);
    }

    #The weights
	w=runif(n)^skew;
    w=w/max(w);                 #w[i] is in [0,1]
    w=w*(max_weight_ratio-1)+1; #w[i] is in [1,max_weight_ratio]
    w=w/sum(w);                 #w[i] is normalised (but still has max(w)/min(w)=max_weight_ratio)
	#Note that skew and max_weight_ratio are not independent
	#Note than n influences max(w) at the start, so will also influence stuff
	#They could be used to find independent variables however
    
	#The data
	R=1:n;

    #The main calculation!
    if(n_cores == 1){
        relative_errors = sort(replicate(n_outer,aws_w(w,k,sample(R,m,prob=w,replace=FALSE),n_inner)$"relative_error"));
    } else {
        registerDoParallel(cores=n_cores);
		relative_errors = sort(unlist(
			foreach(i=1:n_outer) %dopar% aws_w(w,k,sample(R,m,prob=w,replace=FALSE),n_inner)$"relative_error"
			));
    }
    
}

#Obtain parameters
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!= 8) {
    print("ERROR: Not the right number of arguments!");
    quit();
}
n 					= as.numeric(args[1]);
k 					= as.numeric(args[2]);
m 					= as.numeric(args[3]);
skew 				= as.numeric(args[4]);
max_weight_ratio 	= as.numeric(args[5]);
n_outer 			= as.numeric(args[6]);
n_inner 			= as.numeric(args[7]);
n_cores 			= as.numeric(args[8]);

print(args);

#Call function
relative_errors = get_relative_errors_given_distribution(n, k, m, skew, max_weight_ratio, n_outer, n_inner, n_cores);

#Write full data to stdout
print(relative_errors);

show_sigma_levels(relative_errors);
