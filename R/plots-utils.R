.qqplot.pmdadj <- function(object, cutoff=0.99, ...){
	#require(ggplot2)
	xlab <- "expected"
	ylab <- "observed"
	main="QQ plot"
	p <- object$p
	pmd.adj <- object$pmd.adj
	id <- 1:length(pmd.adj)

	## identify any rows with completely missing data
	## or with infinite adj. pmd
	id.na <- which( is.na(pmd.adj))
	pmd.na <- pmd.adj[ id.na]
	##id.inf <- which( is.infinite(pmd.adj) )
	pmd.adj <- pmd.adj[ setdiff(id, id.na) ] 
	id <- id[ setdiff(id, id.na) ] 

	## plot
	n <- length(pmd.adj)
	id.order <- order(pmd.adj)
	pmd.adj <- pmd.adj[id.order]
	id <- id[id.order]
	pmd.exp <- qchisq((1:n - 0.5)/n, df=p)
	threshold <- qchisq( cutoff, df=p)
	outliers <- factor(ifelse(pmd.adj > threshold, 2, 1), levels=c(1,2), labels=c("N","Y"))

	## determine the min and max value of x and y
	## for plotting the thresholds
	min.x <- max( min(pmd.exp, na.rm=T), 1e-9 ); max.x <- max(pmd.exp, na.rm=T)
	min.y <- max( min(pmd.adj, na.rm=T), 1e-9 ); max.y <- max(pmd.adj, na.rm=T)
	min.xy <- min( min.x, min.y)
	max.xy <- min( max.x, max.y) 

	
	## main plot
	myqqplot <- qplot( x=pmd.exp, y=pmd.adj, xlab=xlab, ylab=ylab, main=main, color=outliers)
	## add threshold
	myqqplot <- myqqplot + geom_line(aes(x=c(min.xy, max.xy), y=c(threshold, threshold)), col="red", lty=2)
	## add identity line
	myqqplot <- myqqplot + geom_line(aes(x=c(min.xy, max.xy), y=c(min.xy, max.xy)), col="red", lty=1)
	## add coord_trans
	print( myqqplot + scale_colour_grey(end=0.7) + coord_trans(...) + theme(legend.position="top", 
                legend.title = element_text(size = 10, hjust = 3, vjust = 7)))	
}

.distplot.pmdadj <- function(object, cutoff=0.99, ...){
	#require(ggplot2)
	xlab <- "case number"
	ylab <- "observed adjusted squared distances"
	main="Index plot"
	p <- object$p
	pmd.adj <- object$pmd.adj
	id <- 1:length(pmd.adj)

	## identify any rows with completely missing data
	## or with infinite adj. pmd
	id.na <- which( is.na(pmd.adj))
	pmd.na <- pmd.adj[ id.na]
	##id.inf <- which( is.infinite(pmd.adj) )
	pmd.adj <- pmd.adj[ setdiff(id, id.na) ] 
	id <- id[ setdiff(id, id.na) ] 

	## plot
	threshold <- qchisq( cutoff, df=p)
	outliers <- factor(ifelse(pmd.adj > threshold, 2, 1), levels=c(1,2), labels=c("N","Y"))
	mydistplot <- qplot( x=id, y=pmd.adj, xlab=xlab, ylab=ylab, main=main, color=outliers)
	print(mydistplot + geom_abline(intercept=threshold, 
		slope=0, col="red", lty=2) + scale_colour_grey(end=0.7) + coord_trans(...) + theme(legend.position="top", 
                legend.title = element_text(size = 10, hjust = 3, vjust = 7)))
}

.ddplot.pmdadj <- function(object, cutoff=0.99, ...){
	#require(ggplot2)
	xlab <- "classical distances"
	ylab <- "robust distances"
	main="Distance-distance plot"
	p <- object$p
	pmd.adj <- object$pmd.adj
	id <- 1:length(pmd.adj)

	## obtain classical distances
	pmd.adj.EM <- CovEM(object$x, print.step=0)
	pmd.adj.EM <- pmd.adj.EM@pmd.adj
	
	## identify any rows with completely missing data
	## or with infinite adj. pmd
	id.na <- which( is.na(pmd.adj))
	pmd.na <- pmd.adj[ id.na]
	##id.inf <- which( is.infinite(pmd.adj) )
	pmd.adj <- pmd.adj[ setdiff(id, id.na) ] 
	pmd.adj.EM <- pmd.adj.EM[ setdiff(id, id.na)]
	id <- id[ setdiff(id, id.na) ] 

	## determine the min and max value of x and y
	## for plotting the thresholds
	min.x <- max( min(pmd.adj.EM, na.rm=T), 1e-9 ); max.x <- max(pmd.adj.EM, na.rm=T)
	min.y <- max( min(pmd.adj, na.rm=T), 1e-9 ); max.y <- max(pmd.adj, na.rm=T)
	min.xy <- min( min.x, min.y)
	max.xy <- min( max.x, max.y) 
	
	## plot
	threshold <- qchisq( cutoff, df=p)
	outliers <- factor(ifelse(pmd.adj > threshold | pmd.adj.EM > threshold, 2, 1), levels=c(1,2), labels=c("N","Y"))
	
	## main plot
	myddplot <- qplot( x=pmd.adj.EM, y=pmd.adj, xlab=xlab, ylab=ylab, main=main, color=outliers)
	## add threshold
	myddplot <- myddplot + geom_line(aes(x=c(min.xy, max.xy), y=c(threshold, threshold)), col="red", lty=2)
	myddplot <- myddplot + geom_line(aes(x=c(threshold, threshold), y=c(min.xy, max.y)), col="red", lty=2)	
	## add identity line
	myddplot <- myddplot + geom_line(aes(x=c(min.xy, max.xy), y=c(min.xy, max.xy)), col="red", lty=1)
	## add coord_trans
	print( myddplot + scale_colour_grey(end=0.7) + coord_trans(...) + theme(legend.position="top", 
                legend.title = element_text(size = 10, hjust = 3, vjust = 7)) )
}
