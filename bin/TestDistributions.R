library(fitdistrplus) # evaluate distribution

# Define function to be used to test, get the log lik and aic
tryDistrib <- function(x, distrib){
# deals with fitdistr error:
fit <- tryCatch(MASS::fitdistr(x, distrib), error=function(err) "fit failed")
return(list(fit = fit,
		loglik = tryCatch(fit$loglik, error=function(err) "no loglik computed"),
		AIC = tryCatch(fit$aic, error=function(err) "no aic computed")))
}

findGoodDist <- function(x, distribs, distribs2){
l =lapply(distribs, function(i) tryDistrib(x, i))
names(l) <- distribs
print(l)
listDistr <- lapply(distribs2, function(i){
			       if (i %in% "t"){
			       fitdistrplus::fitdist(x, i, start = list(df =2))
			       } else {
			       fitdistrplus::fitdist(x,i)
			       }}
			       )
par(mfrow=c(2,2))
denscomp(listDistr, legendtext=distribs2)
cdfcomp(listDistr, legendtext=distribs2)
qqcomp(listDistr, legendtext=distribs2)
ppcomp(listDistr, legendtext=distribs2)
par(mfrow=c(1,1))
}

# The density plot and the CDF plot may be considered as the basic classical
# goodness-of-fits plots. The Q-Q plot emphasizes the lack-of-fit at the
# distribution tails while the P-P plot emphasizes the lack-of-fit at the
# distribution center. The nbinom distribution could be prefered for its better
# description of the tail of the empirical distribution.
#x <- pinwormsdata_bal$Aspiculuris_Syphacia
#x <- x[x>0]
#hist(x, breaks = 100)
#descdist(x)
#findGoodDist(x, distribs = c("normal", "negative binomial", "poisson"),
#             distribs2 = c("norm", "nbinom", "pois"))
