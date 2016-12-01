show_sigma_levels = function(rel_errors, sigmas=c(0.8,0.9,0.95,0.99,0.999)  ) {
    invisible(sapply(sigmas,
           FUN=function(sigma) {
			   epsilon = rel_errors[round(sigma*length(rel_errors))];
			   if(epsilon >= 0.01) {
               	   cat(paste("Found a (sigma = ",sigma,", epsilon = ",round(epsilon*100,1),"%) approximation\n",sep=""));
			   } else if(epsilon >= 1e-6) {
               	   cat(paste("Found a (sigma = ",sigma,", epsilon = ",round(epsilon*100,3),"%) approximation\n",sep=""));
			   } else {
               	   cat(paste("Found a (sigma = ",sigma,", epsilon = ",round(epsilon*1e6,3),"e-6) approximation\n",sep=""));
			   }
           }
           ));
}

comp_friendly_sigma_levels = function(rel_errors, sigmas=c(0.8,0.9,0.95,0.99,0.999)  ) {
    cat(paste(sigmas, sep=" "));
	cat(" ");
    invisible(sapply(sigmas,
           FUN=function(sigma) {
			   epsilon = rel_errors[round(sigma*length(rel_errors))];
			   cat(paste(epsilon," ", sep=""));
           }
           ));
}

get_friendly_sigma_levels = function(rel_errors, sigmas=c(0.8,0.9,0.95,0.99,0.999)  ) {
	result="";
    result=paste(result,(paste(sigmas, sep=" ")));
	result=paste(result,(" "));
   	for(sigma in sigmas) { 
	   epsilon = rel_errors[round(sigma*length(rel_errors))];
	   result=paste(result,(paste(epsilon," ", sep="")));
    }
	return(result);
}


get_sigma_levels = function(rel_errors, sigmas=c(0.8,0.9,0.95,0.99,0.999) ) {
    rel_errors = sort(rel_errors);
	sig_levels = sapply(sigmas,
           FUN=function(sigma) {
               return( c(sigma,rel_errors[round(sigma*length(rel_errors))]) );
           }
           );
	return(sig_levels);
}

get_alpha_levels = function(rel_errors, alphas=c(0.2,0.1,0.05,0.01,0.001) ) {
    rel_errors = sort(rel_errors);
	alpha_levels = sapply(alphas,
           FUN=function(alpha) {
               return( c(alpha,rel_errors[round((1-alpha)*length(rel_errors))]) );
           }
           );
	return(alpha_levels);
}
