###Dump of weighed linear sampling code. Stopped developement because;
##-Did not work yet (approximations are nowhere near to correct values)
##-Paper that describes these methods is very vague and theoretical
##-The code below is pretty slow. Optimisation should be possible, but is not trivial
##	> already sped up compute_birthday_period by approximation
##	> calls to stream_sample for retrival of single elements is the current bottleneck
##		- could be fixed by cheating; store full join and get random elements from this
##			> in R this was not enough due to overhead of function call and/or sampling implementation
##				- could be fixed by rewriting this to ``consume'' a list of results

##depedency needed for gamma function related stuff
#install.packages("pracma",repos="http://cran.r-project.org")
#library("pracma")


compute_birthday_period = function(B) {
    print(paste("bday period",B));
    if(B<2000) {
        b = (B-(1:B))/B;
        a = unlist(lapply(1:B, function(i) prod(b[1:i-1])));
        return(sum(a));
    } else {
        return(1.25*sqrt(B));#rough approximation, relative precision ~0.001 in first 1'000'000 values (I think, check with Mathematica)
    }
    #This function can also be expressed as; -1 + B * E^B * ExpIntegralE[-B, B]
    #See Wolfram World for more info
    #This implementation is O(B^2)
    #Can be precomputed, or expressed as (Taylor) series
}

#exp_integral_e_1 = function(x, epsilon) {
#    if(x < 30) {
#        incomplete_gamma = gammainc(x,0);
#    } else {#this term will be < 1e-14
#        incomplete_gamma = 0;
#    }
#    EulerMacheroni = 0.577215664901533;
#    return(-EulerMascheroni-log(x)-incomplete_gamma)
#}
#exp_integral_e = function(n, x) {
#    
#}

#returns TRUE if <= b buckets
#returns FALSE if >= b(1+epsilon) buckets
bucket_number = function(b, epsilon, delta, binary_linear_weighted_sampler) {
    r = compute_birthday_period(b);
    
    #TODO: constants
    #c is a constant such that; E(r((1+epsilon)B)) > (1+c*epsilon)E(r(B))
    #What is c??
    c = 0.5*1/epsilon*(compute_birthday_period((1+epsilon)*b)/r-1);
    #What are c1, c2??
    c1 = 1;#??
    c2 = 1;#??
    #What is maxl? approx sqrt(sum)?
    maxl = 20;
    
    k1 = round(c1 * log(1/delta));
    k2 = round(c2 / (epsilon^2));
    print("ks")
    print(k1)
    print(k2)
    s = c();
    for(i in 1:k1) {
        vec_r = c();
        for(j in 1:k2) {
            #sample until repeated bucket, let rj be the # of samples
            ids_found = c();
            for(l in 1:maxl) {
                #see call in LWSE for binary_linear_weighted_sampler definition
                ids_found = append(ids_found, binary_linear_weighted_sampler());
                
                #terminate if we repeat ``bucket''
                if(length(unique(ids_found)) != length(ids_found)) {
                    vec_r[j] = l;#should this be l or l-1?
                    break;
                }
            }
        }
        rhs = (1+c*epsilon/2)*r;
        lhs = mean(vec_r);
        s[i] = (lhs <= rhs);
    }
    #return TRUE iff more than half of si is TRUE
    return(mean(s) > 0.5);
}

LWSE = function(n, epsilon,delta, linear_weighted_sampler) {
    #TODO: implement linear_weighted_sampler with sample size counter
    #TODO: test functions (correctness, determine constants)
    #TODO: optimize for better performance
    
    #TODO:What is maxk??
    maxk=n;#??
    epsilon_one = epsilon/3;
    L = linear_weighted_sampler()$values;
    for(k in 0:maxk) {
        print(k);
        flush.console();
        a = L*(1+epsilon_one)^k/n;
        b = n*(1+epsilon_one)/epsilon_one;
        buck_val = bucket_number(b,epsilon,delta, 
                        function() {
                            while(TRUE) {
                                lwsample = linear_weighted_sampler();#list(id, value)
                                s = lwsample$values;
                                zeta = epsilon_one*a;
                                m = floor(s/zeta);
                                r = s-zeta*m;
                                #reject with probability r/s
                                if(runif(1,0,1)>r/s) {
                                    #l = sample(1:m,1);
                                    l = ceil(runif(1,0,m-1e-300));
                                    return(c(lwsample$ids,l));
                                }
                            }
                        }()
                      );
        if(buck_val) {
            return(a*n);
        }
    }
    print("LWSE did not terminate!")
    return(-1);
}

#Dump of experiments;
##LWSE(length(AB),epsilon = 0.5, delta = 0.5, function()stream_sample(sampA,sampB,1,weight_modifier = linear_weight_function))
#LWSE(length(AB),epsilon = 0.5, delta = 0.5, function() {
#    return(
#        (function(x)(c(x,AB[x])))(sample(1:length(AB),1,prob=AB)));
#    }
#    );

###
#kaasa= function() return((function(x)(c(x,AB[x])))(sample(1:length(AB),1,prob=AB)));
#kaasa();
