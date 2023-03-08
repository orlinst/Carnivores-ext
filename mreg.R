mreg<-function(formula,data,binomial=FALSE,iter=10,summary=TRUE)	{
	form <- formula(formula)
	if (! missing(data))	{
		y <- model.frame(form,data)[,1]
		x <- as.matrix(model.frame(form,data)[,-1])
	} else	{
		y <- model.frame(form)[,1]
		x <- as.matrix(model.frame(form)[,-1])
	}
	if (! is.matrix(x) || ncol(x) < 2)
		stop('there must be at least two predictors')
	if (binomial == FALSE)
		y <- scale(y)
	x <- scale(data.frame(x))
	n <- ncol(x)
	m <- x
	m2 <- x
	for (i in 1:iter)	{
		for (j in 1:n)
			m2[,j] <- scale(resid(lm(m[,j] ~ m[,-j]))) + m[,j]
		m <- scale(m2)
	}
	if (summary == TRUE)	{
		cat('\nRows and columns in matrix:',nrow(x),'and',ncol(x),'\n')
		cat('Correlation between raw and residual data matrices:',sprintf('%.4f',cor(as.vector(unlist(x)),as.vector(unlist(m)))),'\n')
		r <- cor(m)
		cat('Mean R2 between residual variables:',mean(r[lower.tri(r)]^2),'\n')
	}
	if (binomial == FALSE)	{
		l <- lm(y ~ .,data=data.frame(m))
		cf <- summary(l)$coefficients
		cf <- cbind(cf[,1],cf[,1]^2,cf[,-1])
		colnames(cf)[1:2] <- c('Beta coef.','Beta^2')
	} else	{
		l <- glm(y ~ .,data=data.frame(m),family='binomial')
		cf <- summary(l)$coefficients
	}
	if (summary == TRUE)	{
		cat('\nCoefficients:\n\n')
		b2 <- sum(as.numeric(cf[,2]))
		for (i in 1:ncol(cf))
			cf[,i] <- sprintf('% #.6g',as.numeric(cf[,i]))
		colnames(cf) <- paste0(' ',colnames(cf))
		print(cf,quote=F)
		if (binomial == FALSE)	{
			cat('\nMultiple R-squared:',sprintf('%.4f',summary(l)$r.squared),'\n')
			cat('Sum of squared betas (usually close to R-squared):',sprintf('%.4f',b2),'\n\n')
		} else	{
			cat('\nPseudo-R-squared:',sprintf('%.4f',1 - l$deviance / l$null.deviance),'\n\n')
		}
	} else	{
		return(list(coefficients = cf, residuals = l$residuals, fitted.values = l$fitted.values, residual.matrix = m))
	}
}




d2 = data.matrix(dataX1[,c(5:10,13,14,17:27)])
p = matrix(NA,ncol(d2)+56,20)
s = matrix(NA,ncol(d2)+56,20)
d = list()

for (j in 1:20) {
  
  
    n = data.matrix(get(paste0("dataX", j))[,c(11)])
    loctype <- table(1:128, n)
    m1 = matrix(unlist(loctype), nrow(loctype), ncol(loctype), dimnames=list(1:128, colnames(loctype)))
    
    n = data.matrix(get(paste0("dataX", j))[,c(12)])
    diet.cat <- table(1:128, n)
    m2 = matrix(unlist(diet.cat), nrow(diet.cat), ncol(diet.cat), dimnames=list(1:128, colnames(diet.cat)))
    
    n = data.matrix(get(paste0("dataX", j))[,c(15)])
    activity <- table(1:128, n)
    m3 = matrix(unlist(activity), nrow(activity), ncol(activity), dimnames=list(1:128, colnames(activity)))
    
    n = data.matrix(get(paste0("dataX", j))[,c(16)])
    soc <- table(1:128, n)
    m4 = matrix(unlist(soc), nrow(soc), ncol(soc), dimnames=list(1:128, colnames(soc)))
    
    d[[j]] = data.matrix(get(paste0("dataX", j))[,c(5:10,13,14,17:27)])
    tempx = cbind(d[[j]], m1, m2, m3, m4)
    d0 = data.matrix(dataX1)
    m = lsolm(d0[,3] ~ ., data = data.frame(tempx),summary=F)$rotated.matrix
    rownames(m) <- data_final$species_name
    cat(colnames(d)[i],'\n')
  for (i in 1:ncol(m)) {
    
    p[i,j] <- summary(glm(as.factor(d0[,3]) ~ m[,i],family='binomial'))$coefficients[2,4]
    s[i,j] <- summary(glm(as.factor(d0[,3]) ~ m[,i],family='binomial'))$coefficients[2,1]
    
    w = as.array(m[,i])
    dimnames(w) <- list(data_final$species_name)
    phylop[i,j] <- summary(phyloglm(data_final$extant ~ m,method="logistic_MPLE", phy = tree[[1]]))$coefficients[2,4]
  }
}



cbind(rowMedians(p), colnames(m))