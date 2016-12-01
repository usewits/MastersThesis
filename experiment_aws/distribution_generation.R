
calc_sparsity = function(l) {
    return(length(unique(sort(l)))/length(l));
}

calc_skew = function(l) {
    # variance of a uniform sample is S^2/n(1-f)
    # S^2 = sum (y_i - Y_avg)^2/(N-1)
    l_avg = mean(l);
    #l_s_sq = sum((l-l_avg)^2)/(length(l)-1);
    #return(l_s_sq/length(l));
    return(sum(((l-l_avg)/l_avg)^2)/length(l));#note that n, N and f are constant
}

linear_to_v_B = function(linear_B) {
    return(3^linear_B);
}

linear_to_v_alpha = function(linear_alpha) {
    return(2^linear_alpha);
}

get_sparsity_geom = function(vars) {
    v_alpha = vars[1]; #exponential distribution
    v_B = vars[2];     #rounding factor
    v_n = vars[3];     #amount sampled
    
    #create a list of measured sparsity and skew for different samples
    
    ###results = replicate(n_replications,(function(x) c(calc_sparsity(x),calc_skew(x)))(round(rexp(v_n, v_alpha)*v_B)/v_B));
    results = replicate(n_replications,(function(x) c(calc_sparsity(x),calc_skew(x)))(round(runif(v_n, 0, 1)^v_alpha*v_B)/v_B));
    
    return(c(v_alpha,v_B,v_n,mean(results[1,]),mean(results[2,])));
}

#This function should give 0 once we are close enough to some target sparsity and skew
#Could ``linearize'' B and alpha using parametrizations I found before (see above cell)
root_finding_function = function(linear_B, linear_alpha, v_n, target_sparsity, target_skew, n_replications) {
    v_B = linear_to_v_B(linear_B);
    v_alpha = linear_to_v_alpha(linear_alpha);
    #print(paste("params; ", linear_B, ", ", linear_alpha));
    ###results = replicate(n_replications,(function(x) c(calc_sparsity(x),calc_skew(x)))(round(rexp(v_n, v_alpha)*v_B)/v_B));
    results = replicate(n_replications,(function(x) c(calc_sparsity(x),calc_skew(x)))(round(runif(v_n, 0, 1)^v_alpha*v_B)/v_B));
                                        
    round(runif(v_n, 0, 1)^v_alpha*v_B)/v_B
    mean_sparsity = mean(results[1,]);
    mean_skew = mean(results[2,]);
    #print(c(mean_sparsity,mean_skew));
    #distance = (mean_sparsity - target_sparsity)^2 + (mean_skew - target_skew)^2;
    #to make this converge, we could set distance to 0 if it is smaller than some epsilon here.
    #Epsilon should be slightly bigger than the uncertainty of our measurement.
    #print(c((mean_sparsity - target_sparsity)^2, (mean_skew - target_skew)^2));
    return(((mean_sparsity - target_sparsity)/target_sparsity)^2+((mean_skew - target_skew)/target_skew)^2);
}



get_sample = function(n, sparsity, skew, time) {#This is difficult because of the sparcity!
    if(skew>50*(1.5/50.0)^sparsity) {
        print("WARNING: difficult sparsity/skew relation! Try lowering sparsity or skew.");
    }
    v_n = n;
    target_sparsity = sparsity;
    target_skew = skew;
    n_replications = 1;
    
    n_iterations_available = round(2000*time/(n/1000));
    n_iterations_sim_annealling = round(0.5 * n_iterations_available);
    n_iterations_brute_force = round(0.5 * n_iterations_available);
    
    tmp = function(linear_B_and_alpha) root_finding_function(linear_B_and_alpha[1], linear_B_and_alpha[2], v_n, target_sparsity, target_skew, n_replications);
    solution_linear_B_and_alpha = sim_annealling_find_root(c(0,0),c(12,6),c(1,3),tmp,n_iterations_sim_annealling);
    v_B = linear_to_v_B(solution_linear_B_and_alpha[1]);
    v_alpha = linear_to_v_alpha(solution_linear_B_and_alpha[2]);
    
    best_result = brute_force_find_root(
            c(v_B, v_alpha, v_n),
            ###function(params) round(rexp(params[3], params[2])*params[1])/params[1],
            function(params) round(runif(params[3], 0, 1)^params[2]*params[1])/params[1],
            function(x) ((calc_sparsity(x)-target_sparsity)/target_sparsity)^2+((calc_skew(x)-target_skew)/target_skew)^2,
            n_iterations_brute_force
        );
    
    print(paste("v_B = ", v_B));
    print(paste("v_alpha = ", v_alpha));
    found_sparsity = calc_sparsity(best_result);
    found_skew = calc_skew(best_result);
       
    score = (target_sparsity-found_sparsity)^2+(target_skew-found_skew)^2;
    
    if(score > 0.01) {
        print("WARNING: result did not converge! Try raising the timelimit.");
    }
    
    print(paste("score = ", score, ". Results; sparsity=",found_sparsity,", skew=",found_skew));

    
    return(best_result);
}

sample_assym_sqrt = function(s) {
    sampTab = as.list(table(s));
    left_freqs = unlist(lapply(sampTab, function(x) divisor_closest_to_sqrt(x)), use.names=FALSE);
    right_freqs = unlist(lapply(sampTab, function(x) x/divisor_closest_to_sqrt(x)), use.names=FALSE);
    values = sort(unique(s));
    left_sqrt = unlist(mapply(function(left_freq, value) replicate(left_freq, value), left_freqs, values));
    right_sqrt = unlist(mapply(function(right_freq, value) replicate(right_freq, value), right_freqs, values));
    #Notice that full_join(left_sqrt, right_sqrt) equals sort(sampA)
    return(list(l=left_sqrt,r=right_sqrt));
}
 

