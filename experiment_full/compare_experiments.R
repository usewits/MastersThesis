#Single Cell
source("distribution_generation.R");
source("sampling_utils.R");
library(doParallel)
library(inline)
library(memoise)
library(dplyr)

src_filename = "custom_weighted_sum_join_body.cpp";
include_filename = "naive_weighted_sum_join_includes.cpp";
src = readChar(src_filename, file.info(src_filename)$size);
include = readChar(include_filename, file.info(include_filename)$size);

cat("Compiling C++...\n")
get_approx_join_aggr = cxxfunction(
                        signature(R_R1A="NumericVector",R_R1B="NumericVector",
                                  R_R2A="NumericVector",R_R2C="NumericVector",
                                  R_prob="NumericVector",
                                  R_m="integer", R_use_stratum_correction="numeric", R_delta="numeric",
                                  R_n_replications="integer",
								  R_C_weight_contribution="integer", R_aggregation_index="integer",
								  R_aggregation_approach="integer", R_sampling_approach="integer"),
                                         body = src, includes = include, plugin="Rcpp");

C_UNIFORM=1
C_LINEAR=2
C_INVERSE=3

C_weight_contribution_names = c( "C_UNIFORM",
								 "C_LINEAR",
								 "C_INVERSE")

UNIFORM_AGGREGATION=1
KEY_STRATA_AGGREGATION=2
RANGE_STRATA_AGGREGATION=3
WEIGHTED_AGGREGATION=4

aggregation_approach_names = c( "UNIFORM_AGGREGATION",
                                "KEY_STRATA_AGGREGATION",
                                "RANGE_STRATA_AGGREGATION",
                                "WEIGHTED_AGGREGATION");

U_SAMPLE_JOIN=1
WS_SAMPLE_JOIN=2
AWS_SAMPLE_JOIN=3

sampling_approach_names = c("U_SAMPLE_JOIN",
                            "WS_SAMPLE_JOIN",
                            "AWS_SAMPLE_JOIN");

COL_A=1
COL_B=2
COL_C=3
COL_X=4

aggregation_index_names = c("COL_A",
                            "COL_B",
                            "COL_C",
                            "COL_X");


get_distribution = function(n, skew, ratio=0, n_discrete=0) {
	if(ratio == 0) {#In this case ratio is chosen to get integers in the range 1,n_discrete
		if(n_discrete == 0) {
			cat("ERROR: must set at least either ratio or n_discrete\n");
			return();
		}
		ratio = n_discrete;	
	}
    w=runif(n)^skew;
    w=w/max(w);                 #w[i] is in [0,1]
    w=w*(ratio-1)+1; #w[i] is in [1,max_weight_ratio]
    w=w/sum(w);                 #w[i] is normalised (but still has max(w)/min(w)=max_weight_ratio)
    #Note that skew and max_weight_ratio are not independent
    #Note than n influences max(w) at the start, so will also influence stuff
    #They could be used to find independent variables however

	if(n_discrete > 0) {
		w=w*n_discrete/max(w);#w[i] is in [1, n_discrete]
		w=round(w);#w is discrete
	}
    return(w);
}


get_aggregations = function(R1, R2, J_full, w_index, m, n_inner,
                            aggregation_index = COL_X, delta = 0.01, output_str = "") {

    n1 = length(R1$A);
	
	n1ones = replicate(n1,1);
	
	w_types = list(NULL,R1$A, R1$B,R1$A, R1$B,R1$A, R1$B, n1ones,n1ones,n1ones);
	w_index_names = c("0-ERROR weights","R1A (C_UNIFORM)","R1B (C_UNIFORM)","R1A (C_LINEAR)",
					"R1B (C_LINEAR)","R1A (C_INVERSE)","R1B (C_INVERSE)","UNIFORM (C_UNIFORM)",
					"UNIFORM (C_LINEAR)","UNIFORM (C_INVERSE)");
	use_stratum_correction = FALSE;
	if(w_index <= 3 || w_index == 8) {
		C_weight_contribution = C_UNIFORM; 
	}
	if(w_index == 8 || w_index == 9 || w_index == 10) {
		use_stratum_correction = TRUE;
	}
	if(w_index == 4 || w_index == 5 || w_index == 9) {
		use_stratum_correction = TRUE;
		C_weight_contribution = C_LINEAR; 
	}
	if(w_index == 6 || w_index == 7 || w_index == 10) {
		use_stratum_correction = TRUE;
		C_weight_contribution = C_INVERSE; 
	}

	if(w_index == 1) {
		#####ZERO ERROR WEIGHTS (COL_B)
			##key/((prob[S_index]/R2stratcounts[key])/output_probability_normalization)
			##Assuming all strata join:
			##key/((prob[S_index]/R2stratcounts[key]))
			##B/((prob[S_index]/R2stratcounts[key]))
			##Let prob[S_index] = B*R2stratcounts[key]/sum(B)
			##B/((B*R2stratcounts[key]/sum(B))/R2stratcounts[key])
			##1/((R2stratcounts[key]/sum(B))/R2stratcounts[key])
			##1/(1/sum(B))
			##sum(B)
		if(aggregation_index != COL_B) {
			cat("ERROR: 0-error aggregation not implemented for aggregation indices other than COL_B!\n");
			return();
		}
		w = 1:n1
		cat(paste("obtaining 0-error estimates (",round(n1*n2/1e6,1),"M ops)\n",sep=""));
		sum_J_full_B = sum(J_full$B);
		for(i in c(1:n1)) {
			stratum=R2A[R2A == R1$A[i]];
			w[i]=R1B[i]*sum(R2A == R1$A[i])/sum_J_full_B;
			#Here sum(R2A == R1A[i]) = R2stratcounts[key]
		}
	} else {
		w = w_types[w_index][[1]];
	}

    #Weights on R1 for the (A)WS-step (relevant when using sampling_approach 2 or 3),
    prob = w/sum(w);#Note that our C++ code assumes this is normalized
    w_ratio = max(prob)/min(prob);

	AWS_samp_fraction = ceiling(m*m*w_ratio)/n1;
    cat(paste("Default AWS-sampling fraction is ", AWS_samp_fraction,"\n"));

	result = "";

    #TODO: RANGE_STRATA_AGGREGATION
	#TODO: think of heuristics to select weight function, run experiments with filter, etc
	aggregation_column = c("A","B","C","X")[aggregation_index];
	sum_true = sum(J_full[,aggregation_column]);

	for(aggregation_approach in c(WEIGHTED_AGGREGATION, KEY_STRATA_AGGREGATION, UNIFORM_AGGREGATION)) {
	for(sampling_approach in c(U_SAMPLE_JOIN, WS_SAMPLE_JOIN, AWS_SAMPLE_JOIN)) {
		if(aggregation_approach == UNIFORM_AGGREGATION && sampling_approach == U_SAMPLE_JOIN) {
			cat("\n---Skipped UNIFORM_AGGREGATION over U_SAMPLE_JOIN---\n");
			next;
		}
		if(sampling_approach == AWS_SAMPLE_JOIN && AWS_samp_fraction >= 0.4) {
			cat("\n---Skipped AWS_SAMPLE_JOIN since sampling rate >= 0.4---\n");
			next;
		}
		cat("\n\n------------------------------------------------------\n")
		cat(paste(aggregation_approach_names[aggregation_approach],
				"with", sampling_approach_names[sampling_approach],
				"over",aggregation_index_names[aggregation_index],
				"( and weights w =",w_index_names[w_index],")\n"));
		cat("------------------------------------------------------\n")
		sum_ests = sort(get_approx_join_aggr(R1$A, R1$B, R2$A, R2$C, prob, m, delta, n_inner,
											 R_C_weight_contribution=C_weight_contribution,
											 R_aggregation_index=aggregation_index,
											 R_aggregation_approach=aggregation_approach,
											 R_sampling_approach=sampling_approach,
											 R_use_stratum_correction=use_stratum_correction));
		cat(paste("average estimate",mean(sum_ests),"\n"));
		cat(paste("      true value",sum_true,"\n"));
		relative_errors = sort(abs(sum_ests - sum_true)/sum_true);
		show_sigma_levels(relative_errors);
		cat("\n");

		#cat(paste("CRES ", output_str, " ", w_index," ",aggregation_index," ",aggregation_approach," ",sampling_approach," ",sep="")); comp_friendly_sigma_levels(relative_errors); cat("\n");
		result = paste(result,paste("CRES ", output_str, " ", w_index," ",aggregation_index," ",aggregation_approach," ",sampling_approach," ",sep="")); 
		result = paste(result,get_friendly_sigma_levels(relative_errors),sep="");
	}}
	return(result);
}

n_cores = 8;
#			      1  2  3  4   5   6   7   8
R2C_skewlist =  c(1, 2, 4, 8, 16, 32, 64,128);
R2C_ratiolist = c(1024,1024,1024,1024,1024,1024,1024,1024);

registerDoParallel(cores=n_cores);

foreach(job_id=1:n_cores) %dopar% {

	R2C_skew = R2C_skewlist[as.numeric(job_id)];
	R2C_ratio = R2C_ratiolist[as.numeric(job_id)];

	m = 1000;
	n_inner = 100;

	n1=1e3;
	n2=1000;

	delta = 0.01;

	cat("Generate distributions...\n");

	R1B = get_distribution(n1, skew=8, ratio=1.1);
	R1A = get_distribution(n1, skew=1, n_discrete=100);
	R2A = get_distribution(n2, skew=1, n_discrete=100);
	R2C = get_distribution(n2, skew=R2C_skew, ratio=R2C_ratio);

	R1A = sort(R1A);#without loss of generality we can sort R1A.
	#From this point on, sorting the other columns induces correlations
	#R2A = sort(R2A, decreasing=TRUE);

	#cat(paste("CRES exp_A: uniform 100-discrete join-key. R2C_skew =",R2C_skew,"R2C_ratio =",R2C_ratio,"\n"));

	R1 = as.data.frame(cbind(R1A,R1B));
	colnames(R1) = c("A","B");
	R2 = as.data.frame(cbind(R2A,R2C));
	colnames(R2) = c("A","C");

	est_join_size = max(table(R1A))*max(table(R2A))*max(length(table(R1A)),length(table(R2A)));
	cat(paste("Joining R1 and R2 (estimated size = ",round(est_join_size/1e6,1),"M)...\n",sep=""));
	#cat(paste("Joining R1 and R2 (estimated size = ",est_join_size,")...\n",sep=""));
	J_full = inner_join(R1, R2, by="A");
	cat(paste("(true size = ",round(length(J_full$A)/1e6,1),"M)...\n",sep=""));
	#cat(paste("(true size = ",length(J_full$A),")...\n",sep=""));
	J_full$X=J_full$B*J_full$C+J_full$A; #Add custom ``aggregation function'' here (for attribute X)
	#cat("CRES J_full$B*J_full$C+J_full$A\n")
	#cat(paste("CRES",n1,n2,m,delta,"\n"));

	#Call function
	result = "";
	for(aggregation_index in c(COL_C)) {#c(COL_A, COL_B, COL_C, COL_X)) {
	for(w_index in c(8,9)) {
		partial_res=get_aggregations(R1=R1, R2=R2, J_full=J_full, w_index=w_index, m=m, n_inner=n_inner,
											aggregation_index=aggregation_index, delta=delta,
											output_str = paste(job_id,R2C_skew, R2C_ratio, n1, n2, m,delta,sep=" "));
		result = paste(result,partial_res,sep="");
	}}
	return(result);
}
