minijoin = function(X, Y, joinIndexX=1, joinIndexY=1) {
    lapply(X, function(Xelement) {
        Yjoinable = list();
        for(Yelement in Y) {
            if(Yelement[joinIndexY] == Xelement[joinIndexX]) {
                Yjoinable = c(Yjoinable,list(Yelement));
            }
        }
        if(length(Yjoinable) == 0) {
            return(NULL);
        } else {
            Yelement = sample(Yjoinable,1)[[1]];
            return(c(Xelement,Yelement[-joinIndexY]));
        }
    });
}

stratify_list = function(data, stratum_keys, stratum_attribute_index=1) {
    result = list();
    for(index in stratum_keys) {
        result[as.character(index)][[1]]=list();
    }
    for(element in data) {
        index = as.character(element[stratum_attribute_index]);
        result[index][[1]] = c(result[index][[1]],list(element));
    }
    return(result);
}

#This is the ``naive weighting'' approach
estimate_sum_over_join = function(R1,R2,k,stratum_keys, aggregation_index, join_index) {
    #2) create S \subset R1 (uniform)
    S = sample(R1, k, replace=TRUE);

    #3) let J = S <mini-join> R2
    J_sample = minijoin(S,R2);
    J_full = join(R1, R2, stratum_keys);

    #4) estimate sum from J (assuming strata based on R2, with weights based on counts)
    J_sample_strat = stratify_list(J_sample, stratum_keys, join_index);
    J_strat = stratify_list(J_full, stratum_keys, join_index); #NOTE: we assume strata are in same order!!
    stratum_sizes = sapply(J_strat, length);

    stratum_means = sapply(J_sample_strat, 
                            function(stratum) {
                                if(length(stratum) == 0) {
                                    return(0);
                                } else {
                                    aggregation_column = sapply(stratum, function(element) element[aggregation_index]);
                                    return(mean(aggregation_column));
                                }
                            }
                          );
    sum_est = sum(stratum_means*stratum_sizes);
    return(sum_est);
}
