onesamplemeantest <-
function(X1) {
	#library(msm);
	n1 <- dim(X1)[2];
	p <- dim(X1)[1];
	
	#write("Chen-Qin Two Sample Mean Test:", "");
	#write(paste("X1_n =", n1), "");
	#write(paste("X2_n =", n2), "");
	#write(paste("p =", p), "");
	
	exp1 <- 0;
	T_n <- 0;
	F <- 0;
	pval <- 0;

	for (i in 1:n1) {	
		for (j in 1:n1) {
			if (i != j) {
				term1 <- t(X1[,i])%*%X1[,j];	
				exp1 <- exp1 + term1[1];
			}
		}
	}
	
	#write(paste("expression1 =", exp1), "");

	
	#write(paste("expression3 =", exp3), "");
	
	F_n <- (exp1/(n1*(n1-1)));
	#write(paste("T_n =", T_n), "");

	trS1 <- 0;
	sigma_n1 <- 0;
	
	for (j in 1:n1)	{
		for (k in 1:n1) {
			if (j != k) {
				tempMean <- (rowSums(X1) - X1[,j] - X1[,k])/(n1 - 2);
				term1 <- (X1[,j] - tempMean);
				term1 <- term1 %*% t(X1[,j]);
				term1 <- term1 %*% (X1[,k] - tempMean);
				term1 <- term1 %*% t(X1[,k]);
				trS1 <- trS1 + term1;
			}
		}
	}
	trS1 <- sum(diag(trS1));
	trS1 <- trS1/(n1 * (n1 - 1));

	
	sigma_n1 <- ((2 / (n1 * (n1 - 1))) * trS1);
	sigma_n1 <- sqrt(sigma_n1);
	
	
	Q_n <- F_n/sigma_n1;
	#pval <- ptnorm(Q_n, lower=0, upper=2);
	pval <- 1 - pnorm(Q_n);
	#pval <- pt(Q_n, df=10000);
	return(pval);
}
