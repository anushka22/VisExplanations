rm(list=ls())
library(scagnostics)
library(reshape2)
library(plyr)
library(ggplot2)
library(gridExtra)
library(alphahull)
library(entropy)
library(ape) #---for MST
library(vegan) #---for spantree
set.seed(1357)

#--------
getMSRMeasure <- function(d, variables){
	m = mean(apply(d[,variables], 2, mean, na.rm=TRUE))
	s = mean(apply(d[,variables], 2, sd, na.rm=TRUE))
	model = lm(as.formula(paste(variables[1]," ~ .",sep="")), data=d)
	r = summary(model)$r.squared
	return(c(mean=m, stdev=s, r2=r))
}
#--------
getEntropyMeasure <- function(d, variables, r1, r2, N){
	y2d = discretize2d(d[[variables[1]]], d[[variables[2]]], numBins1=5, numBins2=5, r1=r1, r2=r2) 
	freqs = y2d/N
	H = -sum(ifelse(freqs > 0, freqs * log(freqs), 0))
	return(c(jent=H))
}
#--------
getScagnosticsMeasure <- function(d, variables){
	scagLabels <- c("Outlying"  ,"Skewed"  ,  "Clumpy" ,   "Sparse"  ,  "Striated" , "Convex" ,   "Skinny","Stringy" ,  "Monotonic","AlphaArea","MstLength", "HullPerimeter")
	result <- tryCatch({
		data <- d[,variables]
		data <- data[complete.cases(data),]
		res <- scagnostics(data[,variables])
		rr <- NULL
		for (si in 1:length(res))
			rr <- cbind(rr, res[si])
		rr <- as.data.frame(rr)
		colnames(rr) <- scagLabels
		return(as.data.frame(rr))
	}, error=function(err){
		rr <- t(as.data.frame(rep(0,12)))
		colnames(rr) <- scagLabels
		return(as.data.frame(rr))
	})
}
#--------
getPartitionScores <- function(d, pIndices, variables, r1, r2, N, scoreName){
	groupS <- NULL
	for (i in 1:dim(pIndices)[1]){
		sset = d[pIndices[i,1]:pIndices[i,2],]
		switch(scoreName,
			jent={
				groupS <- rbind(groupS, getEntropyMeasure(sset,variables,r1,r2, N))
			},
			scag={
				groupS <- rbind(groupS, getScagnosticsMeasure(sset, variables))
			},
			msr={
				groupS <- rbind(groupS, getMSRMeasure(sset,variables))
			},
			stop("Enter something that switches me!")
		)
	}
	return(groupS)
}
#--------
savePartitionData <- function(d, partitioner, variables, scoreName,rankValue, dscores, oscores, pIndices){
	nplots <- dim(pIndices)[1]
	ht <- 2*length(nplots)
	pplots = list()
	for (i in 1:nplots){
		sset = d[pIndices[i,1]:pIndices[i,2],]
		name = paste("V",i,sep="")
		pplots[[name]] <- ggplot(sset, aes_string(x=variables[1], y=variables[2])) + geom_point(shape=1) + theme(aspect.ratio = 1)
	}
#	subPlots <- do.call("grid.arrange", c(pplots,  ncol=1, nrow=nplots))

	a = data.frame(cbind(dscores, partitioner=letters[1:nplots]))
	mmin = NULL
	mmax = NULL
	for (i in 1:dim(dscores)[1]){
		mmin = rbind(mmin, min(oscores[i], min(dscores[i,])))
		mmax = rbind(mmax, max(oscores[i], max(dscores[i,])))
	}
	a= melt(a, id="partitioner")
	names(a)[which(names(a)=="partitioner")] = partitioner
	a$value = as.numeric(a$value)
	oscores = data.frame(cbind(oscores, partitioner=letters[1:nplots]))
	names(oscores) = c("values", partitioner)

	#---TODO ---this isn't plotting the correct vlines :S
	formula =paste(partitioner, "~ .",sep=" ")
	allHists <- ggplot(a, aes(x=value)) + geom_histogram() + ggtitle(paste("Partitioned by: ",partitioner,sep="")) 
	subHists <- allHists + facet_grid(formula, scales="free")
	subHists
	subHists <- subHists + geom_vline(data=oscores, aes(xintercept=as.numeric(values)), linetype=3, col="red") + geom_text(data=oscores, aes(x=as.numeric(values),y=0,label=as.character(values)), col="red", size=3) 
	
	#---TODO haven't figured out how to plot these together..
	fpath=paste(path,"simsMaxAbs/",scoreName,"/", rankValue, "-", partitioner, "-", variables[1], "-", variables[2], ".pdf", sep="")
	print(fpath)
	pdf(fpath)
	grid.newpage() # Open a new page on grid device
	pushViewport(viewport(layout = grid.layout(nplots, 1)))
	for (i in 1:nplots){
		print(pplots[i], vp = viewport(layout.pos.row = i, layout.pos.col = 1))
	}
	print(subHists, vp = viewport(layout.pos.row = 1:nplots, layout.pos.col = 2:2))
	dev.off()
	#g <-arrangeGrob(subPlots, subHists, ncol=2, widths=c(4,4), heights=c(ht,ht))
	#ggsave(, g, width=8,height=ht)
}

#----------------
assignToRandomKGroupsBySize <- function(d, ksizes){
	nrows <- dim(d)[1]
	td<-d[sample(1:nrows,nrows,replace=FALSE),]
	gnum <- NULL
	for (ks in ksizes){
		gnum <- c(gnum, rep(ks, ks))
	}
	return(cbind(td, as.factor(gnum)))
}

#----------------
computeKDE <- function(npts, nsplits, x, ptInterest){
	h = (1/npts)^nsplits
	Ku <- NULL
	for (i in 1:npts){
		u <- sum((x[,i]-ptInterest)^2)/(2*h*h)
		Ku <- cbind(Ku, (1/sqrt(2*pi) * exp(-u)))
	}
}
#----------------
computeShingles <- function(d, partitioner){
	t=d[complete.cases(d),]
	z = equal.count(d[[partitioner]])
	plot(z)
	t=levels(z)
	s <- sort(d[[partitioner]])
	inds <- NULL
	for (i in 1:length(t)){
		ss = which(s > t[[i]][1] & s < t[[i]][2])
		inds <- rbind(inds, c(ss[1], ss[length(ss)]))
	}
	return(inds)
}

#----------------
computeAllScores <- function(dall, d, partitioner, variables){
	scNames = c("scag")#, "jent", "msr" )
	NumSamples <- 100
	
	for (scoreName in scNames){
		cat("Computing score: ", scoreName, "\n")
		#---entropy related prep
		r1 <- range(dall[[variables[1]]], na.rm=TRUE)
		r2<- range(dall[[variables[2]]], na.rm=TRUE)
		N <- dim(dall)[1]

		#---split by known partitioner and save scores
		pIndices <- computeShingles(d, partitioner)
		allScores <- getPartitionScores(d, pIndices, variables, r1, r2, N, scoreName)
		#---Get known partitioner split sizes
		sizes <- pIndices[,2]-pIndices[,1] + 1
		cat("SIZES: ",sizes, "\n")
		
		if (scoreName == "msr"){
			vs = c( "mean", "stdev", "r2") 
		}else if (scoreName == "scag"){
			vs = c("Outlying"  ,"Skewed"  ,  "Clumpy" ,   "Sparse"  ,  "Striated" , "Convex" ,   "Skinny","Stringy" ,  "Monotonic","AlphaArea","MstLength", "HullPerimeter")
		}else{
			vs = c("jent")
		} 
		
		#---get scores from simulations
		scores <-list()		#---each row is one split's worth of sim scores
		for (i in 1:NumSamples){
			#---randomly permute
			nrows <- dim(d)[1]
			eg<-d[sample(1:nrows,nrows,replace=FALSE),]
			partScores <- getPartitionScores(eg, pIndices, variables, r1, r2, N, scoreName)
			for (v in vs){
				scores[[v]] <- cbind(scores[[v]], partScores[[v]])
			}
		}
		cat("Saving histograms...\n")
		zscores = list()
		#---loop through the score histograms
		for (v in vs){
			mu = apply(scores[[v]],1, mean, na.rm=TRUE)
			sdev = apply(scores[[v]],1, sd, na.rm=TRUE)
			zscores[[v]] <- max(abs(allScores[[v]] - mu)/sdev, na.rm=TRUE)
			cat("VAR: ",partitioner,"\nv: " , v, " mean: " ,mu, " sdev ", sdev, " zscore: " , zscores[[v]],"\n")
			#---save partition + ranking
			savePartitionData(d, partitioner, variables, v, zscores[[v]], scores[[v]], allScores[[v]], pIndices)
		}
	}	#---for scoreName
}


#######################################
#######################################
#---where to store data
path <- "/Users/aanand/Documents/RESEARCH Notes/RobustScagnostics/"
datafiles <- c("ourworld.csv", "IPEDS_data.csv", "genClustered.csv", "genDonut.csv", "genMonotonic.csv", "genQuadratic.csv", "genSerial.csv", "genStriatedChopped.csv", "genStriatedLines.csv", "adult.csv", "baseball2004.csv", "BOSTON.csv", "Sample - Superstore Sales.csv")

datafiles <- c("baseball2004.csv")
df = datafiles[1]
#for (df in datafiles){
	fname <- paste(path, "data/", df, sep="")
	dall <- read.csv(fname, header=TRUE)

	#---make binned fields factors
	bs=ifelse(grepl("(bin)",names(dall)), TRUE, ifelse(grepl("(cluster)",names(dall)), TRUE, FALSE))
	if (length(bs[which(bs==TRUE)]) > 0)
		dall[[names(dall)[bs]]] <- as.factor(dall[[names(dall)[bs]]])
	#---remove factors
	dfactors <- sapply(dall, is.factor)
	ps <- names(dall[dfactors])
	allVars <- names(dall)[!(names(dall) %in% ps)]
	
	#---for all pairs of vars that are Qs
	for (i in 1:length(allVars)){
		for (j in (i+1):length(allVars)){
			if (j > length(allVars))
				break
			for (k in 1:length(allVars)){
				if (k == i | k== j)
					next
				partitioner=allVars[k]
				variables = c(allVars[i], allVars[j])
				cat("Processing: ", variables, " and partitioner: ", partitioner, "\n")
				d <- dall[,c(variables,partitioner)]
				computeAllScores(dall, d, partitioner, variables)
			}	#---for k
		}	#---for j
	}	#---for i
#}

