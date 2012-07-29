twosamplemeantest <-
function(X1,X2) {
	#library(msm);
	n1 <- dim(X1)[2];
	n2 <- dim(X2)[2];
	p <- dim(X1)[1];
	
	#write("Chen-Qin Two Sample Mean Test:", "");
	#write(paste("X1_n =", n1), "");
	#write(paste("X2_n =", n2), "");
	#write(paste("p =", p), "");
	
	exp1 <- 0;
	exp2 <- 0;
	exp3 <- 0;
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

	for (i in 1:n2) {
		for (j in 1:n2) {
			if (i != j) {
				term2 <- t(X2[,i])%*%X2[,j];	
				exp2 <- exp2 + term2[1];
			}
		}
	}
	
	#write(paste("expression2 =", exp2), "");
	
	for (i in 1:n1) {
		for (j in 1:n2) {
			term3 <- t(X1[,i])%*%X2[,j];
			exp3 <- exp3 + term3[1];
		}
	}
	
	#write(paste("expression3 =", exp3), "");
	
	T_n <- (exp1/(n1*(n1-1))) + (exp2/(n2*(n2-1))) - (2*(exp3/(n1*n2)));
	#write(paste("T_n =", T_n), "");

	trS1 <- 0;
	trS2 <- 0;
	trS1S2 <- 0;
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
	#write(paste("trS1 =", trS1), "");

	for (j in 1:n2)	{
		for (k in 1:n2) {
			if (j != k) {
				tempMean <- (rowSums(X2) - X2[,j] - X2[,k])/(n2 - 2);
				term2 <- (X2[,j] - tempMean);
				term2 <- term2 %*% t(X2[,j]);
				term2 <- term2 %*% (X2[,k] - tempMean);
				term2 <- term2 %*% t(X2[,k]);
				trS2 <- trS2 + term2;
			}
		}
	}
	trS2 <- sum(diag(trS2));
	trS2 <- trS2/(n2 * (n2 - 1));
	#write(paste("trS2 =", trS2), "");

	for (j in 1:n1)	{
		for (k in 1:n2) {
			tempMean1 <- (rowSums(X1) - X1[,j])/(n1 - 1);
			tempMean2 <- (rowSums(X2) - X2[,k])/(n2 - 1);
			term3 <- (X1[,j] - tempMean1);
			term3 <- term3 %*% t(X1[,j]);
			term3 <- term3 %*% (X2[,k] - tempMean2);
			term3 <- term3 %*% t(X2[,k]);
			trS1S2 <- trS1S2 + term3;
		}
	}
	trS1S2 <- sum(diag(trS1S2));
	trS1S2 <- trS1S2/(n1 * n2);
	#write(paste("trS1S2 =", trS1S2), "");
	
	sigma_n1 <- ((2 / (n1 * (n1 - 1))) * trS1) + 
				((2 / (n2 * (n2 - 1))) * trS2) + 
				((4 / (n1 * n2)) * trS1S2);
	sigma_n1 <- sqrt(sigma_n1);
	#write(paste("sigma_n1 =", sigma_n1), "");
	
	#psi_alpha <- qnorm(0.95);
	#write(paste("psi_alpha (at alpha = 0.05) =", psi_alpha), "");
	
	Q_n <- T_n/sigma_n1;
	#write(paste("Q_n =", Q_n), "");
	#pval <- ptnorm(Q_n, lower=0, upper=2);
	pval <- 1 - pnorm(Q_n);	
	#pval <- pt(Q_n, df=10000);
	return(pval);
}
