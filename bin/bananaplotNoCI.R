bananaplotNoCI <- function (mod, data, response, hybridIndex = seq(0, 1, 0.05),
				 cols = wesanderson::wes_palette("IsleofDogs1")[c(1, 3)],
				 group, islog10 = F)
{
data$response = data[[response]]
data$group = data[[group]]
getBananaDF <- function(mod, hybridIndex) {
fittedCoef <- bbmle::coef(mod)
if ("L2" %in% names(fittedCoef) == FALSE) {
fittedCoef <- c(fittedCoef, fittedCoef[names(fittedCoef) %in%
					    "L1"])
names(fittedCoef)[length(names(fittedCoef))] <- "L2"
}
myConfInt <- bbmle::confint(mod, method = "quad")
if ("L2" %in% rownames(myConfInt) == FALSE) {
myConfInt <- rbind(myConfInt[rownames(myConfInt) %in%
				     "L1"], myConfInt)
rownames(myConfInt)[1] <- "L2"
}
getInf <- function(paramname) {
myConfInt[rownames(myConfInt) == paramname][1]
}
getSup <- function(paramname) {
myConfInt[rownames(myConfInt) == paramname][2]
}
expectedResponse <- function(v, hybridIndex) {
L1 = v[1]
L2 = v[2]
alpha = v[3]
heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
expectedResponse <- (L1 + (L2 - L1) * hybridIndex) *
(1 - alpha * heterozygoty)
return(expectedResponse)
}
bananaDF = data.frame(HI = hybridIndex)
bananaDF$fit <- expectedResponse(c(fittedCoef[["L1"]],
					     fittedCoef[["L2"]], fittedCoef[["alpha"]]), hybridIndex = hybridIndex)
bananaDF2 = data.frame(HI = numeric(), min = numeric(),
			  max = numeric())
for (i in hybridIndex) {
maxLoad <- stats::optim(par = c(L1 = getSup("L1") -
				   getInf("L1"), L2 = getSup("L2") - getInf("L2"),
				   alpha = getSup("alpha") - getInf("alpha")), fn = expectedResponse,
				   lower = c(L1 = getInf("L1"), L2 = getInf("L2"),
						alpha = getInf("alpha")), upper = c(L1 = getSup("L1"),
										       L2 = getSup("L2"), alpha = getSup("alpha")),
						method = "L-BFGS-B", control = list(fnscale = -1),
						hybridIndex = i)
minLoad <- stats::optim(par = c(L1 = getSup("L1") -
				   getInf("L1"), L2 = getSup("L2") - getInf("L2"),
				   alpha = getSup("alpha") - getInf("alpha")), fn = expectedResponse,
				   lower = c(L1 = getInf("L1"), L2 = getInf("L2"),
						alpha = getInf("alpha")), upper = c(L1 = getSup("L1"),
										       L2 = getSup("L2"), alpha = getSup("alpha")),
						method = "L-BFGS-B", hybridIndex = i)
bananaDF2 = rbind(bananaDF2, data.frame(HI = i, min = minLoad$value,
					   max = maxLoad$value))
}
bananaDF$minAll <- bananaDF2$min
bananaDF$maxAll <- bananaDF2$max
bananaDF2 = data.frame(HI = numeric(), min = numeric(),
			  max = numeric())
for (i in hybridIndex) {
maxLoad <- stats::optim(par = getSup("alpha") - getInf("alpha"),
			    fn = parasiteLoad::MeanLoad, lower = getInf("alpha"),
			    upper = getSup("alpha"), L1 = fittedCoef[["L1"]],
			    L2 = fittedCoef[["L2"]], method = "L-BFGS-B",
			    control = list(fnscale = -1), hybridIndex = i)
minLoad <- stats::optim(par = getSup("alpha") - getInf("alpha"),
			    fn = parasiteLoad::MeanLoad, lower = getInf("alpha"),
			    upper = getSup("alpha"), L1 = fittedCoef[["L1"]],
			    L2 = fittedCoef[["L2"]], method = "L-BFGS-B",
			    hybridIndex = i)
bananaDF2 = rbind(bananaDF2, data.frame(HI = i, min = minLoad$value,
					   max = maxLoad$value))
}
bananaDF$minAlpha <- bananaDF2$min
bananaDF$maxAlpha <- bananaDF2$max
return(bananaDF)
}
if (is.list(mod) == FALSE) {
bananaDFtoplot <- getBananaDF(mod = mod, hybridIndex = hybridIndex)
bananaDFtoplot$group <- "all"
}
else {
bananaList <- lapply(mod, FUN = getBananaDF, hybridIndex = hybridIndex)
bananaDFA <- bananaList$groupA
bananaDFB <- bananaList$groupB
bananaDFA$group <- levels(data$group)[1]
bananaDFB$group <- levels(data$group)[2]
bananaDFtoplot <- rbind(bananaDFA, bananaDFB)
}
p <- ggplot2::ggplot() + ggplot2::geom_ribbon(data = bananaDFtoplot,
						   ggplot2::aes_string(x = "HI", ymin = "minAll", ymax = "maxAll",
									 group = "group"), fill = "grey", alpha = 0.4) + ggplot2::geom_line(data = bananaDFtoplot,
																		 size = 2, ggplot2::aes_string(x = "HI", y = "fit", col = "group")) +
ggplot2::geom_point(data = data,
			 ggplot2::aes_string(x = "HI", y = "response", fill = "group"),
			 pch = 21, size = 3, alpha = 0.5) +
ggplot2::scale_fill_manual(values = cols) +
ggplot2::theme_classic(base_size = 20) + {
if (islog10 == TRUE)
ggplot2::scale_y_log10()
} + ggplot2::ylab(label = response)
if (is.list(mod) == TRUE) {
p <- p + ggplot2::scale_color_manual(values = cols)
}
else {
p <- p + ggplot2::scale_color_manual(values = "black")
}
return(p)
}
