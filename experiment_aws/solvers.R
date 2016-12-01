sim_annealling_find_root = function(x_min, x_max, x_start, #vectors
									score_function, niters, 
									initial_step_size_factor = 0.1,
									step_size_reduction_factor = 0.9999,
									initial_temperature = 1,
									temperature_reduction_factor = 0.9) {
    
    step_x = (x_max - x_min)*initial_step_size_factor;
    
    x_curr = x_start;
    score_curr = score_function(x_curr);
    score_best = score_curr;
    x_best = x_curr;
    
    for(iter in 1:niters) {
        T = initial_temperature*temperature_reduction_factor^iter;
        x_try = x_curr + runif(2,-1,1)*step_x;   #update vector
        x_try = mapply(min, x_try, x_max);  #keep it within [x_min, x_max]
        x_try = mapply(max, x_try, x_min);
        
        score_try = score_function(x_try);
        
        if(score_try < score_curr || exp(-(score_try-score_curr)/T) > runif(1,0,1) ) {
            x_curr = x_try;
            score_curr = score_try;
            if(score_curr < score_best) {
                x_best = x_curr;
                score_best = score_curr;
            }
        }
        step_x = step_size_reduction_factor*step_x;
    }
    return(x_best);
}

#Function for optimizing randomized value by trying multiple times
brute_force_find_root = function(parameters, randomized_generator, score_function, niters) {
    #assumes we want to minimize score_function
    best_result = -1;
    best_score = -1;
    for(iter in 1:niters) {
        try_result = randomized_generator(parameters);
        try_score = score_function(try_result);
        if(try_score < best_score || best_score < 0) {
            best_result = try_result;
            best_score = try_score;
        }
    }
    return(best_result);
}


