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




d2 = data.matrix(dataX1[,c(3,4,7:13,15:27)])
p = matrix(NA,ncol(d2),20)


for (i in 1:ncol(d2)) {
  for (j in 1:20) {
    d[j] = data.matrix(get(paste0("dataX", j))[,c(3,4,7:13,15:27)])
    d0 = data.matrix(dataX1)
    m = mreg(d0[,5] ~ d[[j]],summary=F)$residual.matrix
    cat(colnames(d)[i],'\n')
    p[i,j] <- summary(glm(as.factor(d0[,5]) ~ m[,i],family='binomial'))$coefficients[2,4]
  }
}

rowMeans(p)
