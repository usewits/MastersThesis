#FORMAT OF FUNCTION NAMES:
#
# [[log_]p_](ws|aws)_(w|wo)[_given_u]
#    ^   ^     ^       ^       ^
#    |   |     |       |       +----------- conditional probability given U
#    |   |     |       |
#    |   |     |       +------------------- w=with replacement, wo=without replacement
#    |   |     |
#    |   |     +--------------------------- ws=weighted sampling, aws=approximate weighted sampling
#    |   |
#    |   +--------------------------------- p_ is prepended if the function only returns a probability
#    |
#    +------------------------------------- log_ is prepended if the returned probability is logarithmic

#  unif. wght.
# R---->U---->S
# n     k     m

library("inline");

src_filename = "sampling_probability_body.cpp";
include_filename = "sampling_probability_includes.cpp";
src = readChar(src_filename, file.info(src_filename)$size);
include = readChar(include_filename, file.info(include_filename)$size);

get_log_prob = cxxfunction(signature(R_Ubar="IntegerVector", R_S="IntegerVector", R_prob="NumericVector"),
                body = src, includes = include, plugin="Rcpp");

count = function(S) {
	return(as.vector(table(S)));
}

n_perms = function(S) {
    return(factorial(length(S))/prod(factorial(count(S))));
}

log_n_perms = function(S) {
    #only works if max(count(S)) is not too big! Will give warning otherwise.
    return(log_factorial(length(S))-sum(log(factorial(count(S)))));
}

p_ws_w = function(S, prob) {
    base_probability = prod(prob[S]/sum(prob));
    return(base_probability*n_perms(S));
}

log_p_ws_w = function(S, prob) {
    log_base_probability = sum(log(prob[S])-log(sum(prob)));
    return(log_base_probability+log_n_perms(S));
}

simulate_p_ws_w = function(S, prob, n_replications) {
    S=sort(S);
    m=length(S);
	n=length(prob);
	
	R=1:n;

    samples = sapply(1:n_replications,
    FUN = function(id) {
        St = sort(sample(R,m,prob=prob,replace=TRUE));
        if(sum(abs(St-S))==0) {
            return(1);
        } else {
            return(0);
        }
    })
    return(mean(samples));
}

p_aws_w_given_u = function(S, Ubar, prob) {
    #return(factorial(length(S))*prod(prob[S]*(count_occurences(Ubar,S)+1)/sum(prob[Ubar],prob[S])));
    return(prod(prob[S]*(count_occurences(c(Ubar,S),S))/sum(prob[Ubar],prob[S]))*n_perms(S))
}

simulate_p_aws_w_given_u = function(S, Ubar, prob, n_replications) {
    S=sort(S);
    m=length(S);
    U=c(S,Ubar);
    Uprob = prob[U];

    samples = sapply(1:n_replications,
    FUN = function(id) {
        St = sort(sample(U,m,prob=Uprob,replace=TRUE));
        if(sum(abs(St-S))==0) {
            return(1);
        } else {
            return(0);
        }
    })
    return(mean(samples));
}

log_p_aws_w_given_u = function(S, Ubar, prob) {
     log_norm_U = log(sum(prob[Ubar],prob[S]));
     return(sum(log(prob[S])+log(count_occurences(c(Ubar,S),S)))-log_norm_U*length(S)+log_n_perms(S));
}

simulate_aws_w_log_factor = function(prob,k,S,n_replications) {
    n=length(prob);
    m=length(S);
    R=1:n;
    Us = replicate(n_replications,sample(R,k,replace=TRUE),simplify=FALSE);
    contains=unlist(lapply(Us,
        FUN=function(U) {
            if(length(setdiff(S,U))==0) {
                return(1);
            } else {
                return(0);
            }
        }));
    exact_n_nonzero = n^(k-m);
    exact_total = n^k;
    print(exact_n_nonzero)
    print(exact_total)
    print(exact_n_nonzero/exact_total)
    return(mean(contains));
}


aws_w = function(prob,k,S,n_replications,cpp=TRUE) {#NOTE: only works if S does not contain duplicates!
    n=length(prob);
    m=length(S);
    log_factor = (k-m)*log(n);#num ordered Ubar's is n^(k-m, log(n^(k-m))=(k-m)*log(n)
    R=1:n;
    
    Ubars = replicate(n_replications,sample(R,k-m,replace=TRUE),simplify=FALSE);#sample from all sets (different order)
    #log_p_ws_w(c(Ubar,S),replicate(n,1)) is the uniform selection probability of U=c(Ubar,S) with replacement (unordered)
    #-k*log(n)+log_n_perms(c(S,Ubar)) is the uniform selection probability of U=c(Ubar,S) with replacement (ordered)
    #need to subtract log_n_perms(Ubar) to correct for the fact that we draw unordered Ubar (order matters)

    if(cpp == TRUE) {
        log_probs = sapply(Ubars, 
                        FUN = function(Ubar) get_log_prob(Ubar,S,prob)
                    );
    } else {
        log_probs = sapply(Ubars, 
                        FUN = function(Ubar) log_p_aws_w_given_u(S,Ubar,prob)-k*log(n)+log_n_perms(c(S,Ubar))-log_n_perms(Ubar)
                    );
    }
    
    log_p_aws_w = log_factor + log_mean(log_probs);
    log_p_ws_w = log_p_ws_w(S,prob);
    log_absolute_error = log_abs_diff(log_p_ws_w,log_p_aws_w);
    log_relative_error = log_absolute_error - log_p_ws_w;
    relative_error = exp(log_relative_error);

    return(list(log_p_aws_w=log_p_aws_w,relative_error=relative_error,log_relative_error=log_relative_error));
}

exact_aws_w = function(prob,k,S) {
     n=length(prob);
    m=length(S);

    R=1:n;
    
    Ubars = generate_ordered_sets(k-m,n);#all sets up to permutation

    return(list(p_aws_w=sum(sapply(Ubars, FUN = function(Ubar) p_aws_w_given_u(S,Ubar,prob)*p_ws_w(c(Ubar,S),replicate(n,1))))));
}

simulation_aws_w = function(prob,k,S,n_replications) {
    S=sort(S);
    n=length(prob);
    m=length(S);
    R = 1:n;
    samples = sapply(1:n_replications,
        FUN = function(id) {
            Ut = sample(R,k,replace=TRUE);
            St = sort(sample(Ut,m,prob=prob[Ut],replace=TRUE));
            if(sum(abs(St-S))==0) {
                return(1);
            } else {
                return(0);
            }
        })
    return(mean(samples));
}


#Note that this is incorrect! Does not take all orders into account!
#p_ws_wo = m! * product(w(s_i)/(w_n-w_[1 ... i-1]), {i,1,m)
p_ws_wo = function(S, prob, norm=sum(prob)) {
	stop("called p_ws_wo!")
    result = 1;
    for(i in c(1:length(prob[S]))) {
        s = prob[S[i]];
        result=result*i*s/norm;#note that i contributes to a factor factorial(length(S))
        norm=norm-s;
    }
    return(result);
}

log_p_ws_wo = function(S, prob, norm=sum(prob)) {
    result = 0;
    for(i in c(1:length(prob[S]))) {
        s = prob[S[i]];
        result=result+log(i*s)-log(norm);#note that i contributes to a factor factorial(length(S))
        norm=norm-s;
    }
    return(result);
}

p_ws_wo_given_u = function(S, Ubar, prob) {
    norm = sum(prob[Ubar])+sum(prob[S]);
    return(p_ws_wo(S,prob,norm));
}

log_p_ws_wo_given_u = function(S, Ubar, prob) {
    norm = sum(prob[Ubar])+sum(prob[S]);
    return(log_p_ws_wo(S,prob,norm));
}

aws_wo = function(prob,k,S,n_replications) {
    #only implemented log version of this function, because typically very small probabilities are involved
    n=length(prob);
    m=length(S);

    if(k>n) {
        print("WARNING: U bigger than R!");
    }
    base_weight = sum(prob[S]);
    indices_not_S = setdiff(c(1:n),S);
    
    #constant_factor = choose(n-m,k-m)/choose(n,k);
    log_constant_factor = sum(log((k-c(0:(m-1)))/(n-c(0:(m-1)))));
   
    log_p_aws_results = replicate(n_replications,
                            (function() {
                                Ubar=sample(indices_not_S,k-m,replace=FALSE);#INNER LOOP
                                return(log_p_ws_wo_given_u(S,Ubar,prob));     #The sampling takes roughly half the time
                            })()
                           );
    
    log_p_aws_result = log_mean(log_p_aws_results)+log_constant_factor;
    
    log_p_ws_wo = log_p_ws_wo(S,prob);
    log_absolute_error = log_abs_diff(log_p_ws_wo,log_p_aws_result);
    log_relative_error = log_absolute_error - log_p_ws_wo;
    relative_error = exp(log_relative_error);

    return(list(log_p_aws_wo=log_p_aws_result,relative_error=relative_error,log_relative_error=log_relative_error));
}

exact_aws_wo = function(prob,k,S) {
    #only implemented log version of this function, because typically very small probabilities are involved
    n=length(prob);
    m=length(S);

    if(k>n) {
        print("WARNING: U bigger than R!");
    }
    
    indices_not_S = setdiff(c(1:n),S);
    #Loop over all size k-m sets in indices_not_S
    all_u_bars = combn(indices_not_S,k-m,simplify=FALSE);
    all_u_bar_weights = lapply(all_u_bars, FUN = function(ubar) {
                                    factor=1;
                                    norm = sum(prob[S])+sum(prob[ubar]);
                                    for(i in 1:m) {
                                        factor=factor*1/(norm);
                                        norm = norm - prob[S[i]];
                                    }
                                    return(prod(prob[S])*factor);
                                }
                            );
    result = factorial(m)/choose(n,k)*sum(unlist(all_u_bar_weights));
    return(result);
}



#for every key in keys, counts number of occurences in data
#output: vector of results
count_occurences = function(data, keys) {
    data_dict = as.list(table(data));
    results = lapply(keys, FUN = function(n) data_dict[as.character(n)][[1]]);
    results[sapply(results,is.null)]=0;
    return(unlist(results));
}

#Generates all ordered sets of length len with maximum value max_val
generate_ordered_sets = function(len, max_val, stub=c()) {
    if(length(stub) == 0) {
        min_val = 1;
    } else {
        min_val = stub[length(stub)];
    }
    if(len > 0) {
        result_list = unlist(lapply(min_val:max_val, FUN = function(val) generate_ordered_sets(len-1,max_val,c(stub,val))),recursive=FALSE);
    } else {
        return(list(stub));
    }
    return(result_list);
}



log_mean = function(log_list) {
    #returns log(mean(exp(log_list))), but allows for high precision calculation if log_list[i] is very negative
    #We use the identity log(exp(a)+exp(b)) = a + log(1 + exp(b - a)) repeatedly
    #The sorting imporves the quality
    #We have to substract log(length(log_list)) to get from log(sum(exp(log_list))) to the desired result
    return(Reduce(function(a,b) {a + log(1 + exp(b - a))},sort(log_list,decreasing = TRUE))-log(length(log_list)));
}

log_abs_diff = function(log_a, log_b) {
    x = min(log_a, log_b);
    y = max(log_a, log_b);
    return(x+log(exp(y-x)-1));
}


log_factorial_raw = function(n) {
    n = as.integer(n);
    if(n > 1) {
        return(sum(log(1:n)));
    }
    if(n == 1) {
        return(0);
    }
    print("ERROR: negative log_factorial!")
}

log_factorial_table = sapply(1:3000, log_factorial_raw);

log_factorial = function(n) {
    n = as.integer(n);
    if(n > 1) {
        if(n <= 3000) {
            return(log_factorial_table[n]);
        } else {
            x=n+1;
            return((x-0.5)*log(x)-x+0.5*log(2*pi)+1/(12*x))-1/(360*x*x*x);
        }
    }
    if(n == 1) {
        return(0);
    }
    print("ERROR: negative log_factorial!")
}


d_prob=c(1,2,3,4,1,2,3,4);
d_S = c(1,2);
d_Ubar = c(3,4);
