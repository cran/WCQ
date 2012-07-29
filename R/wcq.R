wcq <-
function(marker_data, trait_data, alleles) {
#	library(msm);
#	source("twosamplemeantest.R");
#	source("onesamplemeantest.R");
#	source("twosamplemeantest_qn.R");
#	source("onesamplemeantest_qn.R");
	num_samples <- dim(marker_data)[2];
	num_markers <- dim(marker_data)[1];
	num_traits <- dim(trait_data)[1];
	pval_matrix <- array(0, dim=c(num_markers, num_traits));

	#allele_sample_prob <- array(0, dim=c(length(alleles),num_samples));
	#for (i in 1:length(alleles)) {
	#	for (j in 1:num_samples) {
	#		allele_sample_prob[i,j] <- sum(marker_data[,j]==alleles[i]);
	#	}
	#}	
	#allele_sample_prob <- allele_sample_prob / num_markers;
		
	for (i in 1:num_markers) {		
		for (j in 1:num_traits) {
			current_trait <- trait_data[j,];
			mean_current_trait <- mean(current_trait);
			stdev_current_trait <- sqrt(var(current_trait));
			X_indices_1 <- which(marker_data[i,]==alleles[1], arr.ind=T);
			#current_trait[X_indices_1] <- current_trait[X_indices_1] * (allele_sample_prob[1,X_indices_1] * length(X_indices_1) / num_samples);
			X_indices_2 <- which(marker_data[i,]==alleles[3], arr.ind=T);
			#current_trait[X_indices_2] <- 0.75 * current_trait[X_indices_2];
			X_indices_tot <- c(X_indices_1, X_indices_2);
			#X_indices_tot <- c(X_indices_1);
			X <- as.vector(current_trait)[X_indices_tot];
			X <- X / pnorm(X, mean=mean_current_trait, sd=stdev_current_trait);
			#X <- X / pnorm(X);
			
			Y_indices <- which(marker_data[i,]==alleles[2], arr.ind=T);
			#Y_indices <- c(Y_indices, X_indices_2);
			#current_trait[Y_indices] <- current_trait[Y_indices] * (allele_sample_prob[2,Y_indices] * length(Y_indices) / num_samples);
			Y <- as.vector(current_trait)[Y_indices];
			Y <- Y / pnorm(Y, mean=mean_current_trait, sd=stdev_current_trait);
			#Y <- Y / pnorm(Y);
			
			pval <- 0;			
			X_arr <- array(X, dim=c(1,length(X_indices_tot)));
  			Y_arr <- array(Y, dim=c(1,length(Y_indices)));
			if ((length(X) > 1) && (length(Y) > 1)) {
				pval <- twosamplemeantest(X_arr,Y_arr);
				if (is.nan(pval)) { #Chen-Qin test "failed"
					if (length(X) > length(Y)) {
						X_arr_normalized <- X_arr - array(mean(X_arr), c(1, length(X_indices_tot)));
						pval <- onesamplemeantest(X_arr_normalized);
					} else if (length(Y) > length(X)){
						Y_arr_normalized <- Y_arr - array(mean(Y_arr), c(1, length(Y_indices)));
						pval <- onesamplemeantest(Y_arr_normalized);
					} else {
						XY_diff <- X_arr - Y_arr;
						pval <- onesamplemeantest(XY_diff);
					}
				}
			} else { #only one group has at least sample size = 2
				if (length(X) > 1) {
					X_arr_normalized <- X_arr - array(mean(X_arr), c(1, length(X_indices_tot)));
					pval <- onesamplemeantest(X_arr_normalized);
				} else if (length(Y) > 1){
					Y_arr_normalized <- Y_arr - array(mean(Y_arr), c(1, length(Y_indices)));
					pval <- onesamplemeantest(Y_arr_normalized);
				}
			}
			if (is.nan(pval)) { #if all possble tests fail
				pval <- 0;
			}
			# transform and express the result in terms of log(p)
			#pval <- 1 - pval;
			#if (pval != 0) {
			#	pval = -log10(pval);
			#}
			pval_matrix[i,j] <- pval;
		}
	}
	return(pval_matrix);
}
