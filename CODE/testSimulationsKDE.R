rm(list=ls())
library(scagnostics)
library(reshape2)
library(plyr)
library(ggplot2)
library(gridExtra)

set.seed(1357)

#--------
getMSRMeasure <- function(d, variables, partitioner){
	tmp = d[!names(d)==partitioner]
	m = mean(apply(tmp[,variables], 2, mean, na.rm=TRUE))
	s = mean(apply(tmp[,variables], 2, sd, na.rm=TRUE))
	model = NA
	model = tryCatch({
		 lm(as.formula(paste0(variables[1]," ~ .")), data=tmp,na.action=na.omit)
	}, error = function(e){})
	if (class(model) == "lm"){
		r = summary(model)$r.squared
		sl = summary(model)$coefficients[2]
		f <- summary(model)$fstatistic
		if (is.null(f))
			p=0
		else
			p = 1-unname(pf(f[1],f[2],f[3],lower.tail=F))
	}else{
		r = 0
		f = 0
		p = 0
		sl = 0
	}
	return(c(mean=m, stdev=s, r2=r, pval=p, slope=sl))
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
getPartitionScores <- function(d, partitioner, variables, r1, r2, N, scoreName){
	switch(scoreName,
		jent={
			groupS <- ddply(d, partitioner, .fun=function(x,col){ e=getEntropyMeasure(x,col,r1,r2, N)}, variables)
		},
		scag={
			groupS <- ddply(d, partitioner, .fun=function(x, col){s= getScagnosticsMeasure(x, col)}, variables)
		},
		msr={
			groupS <- ddply(d, partitioner, .fun=function(x,col,part){ e=getMSRMeasure(x,col, part)}, variables, partitioner)
		},
		stop("Enter something that switches me!")
	)
	return(groupS)
}
#--------
savePartitionData <- function(d, partitioner, variables, scoreName, groupS, rankValue, dscores, oscores, varpath){
	formula <- paste(partitioner, "~ .",sep=" ")
	allPlots <- ggplot(d, aes_string(x=variables[1], y=variables[2])) + geom_point(shape=1) + theme(aspect.ratio = 1) + theme(axis.title=element_text(size=5)) + theme(axis.text=element_text(size=4)) + geom_smooth(method='lm')
#	ggsave(paste0(path,folder,"/", varpath, ".pdf"), allPlots, width=2,height=2)	

	subPlots <- allPlots + facet_grid(formula)

	a = data.frame(cbind(dscores, partitioner=as.character(levels(d[[partitioner]]))))
	mmin = NULL
	mmax = NULL
	for (i in 1:dim(dscores)[1]){
		mmin = rbind(mmin, min(oscores[i,scoreName], min(dscores[i,])))
		mmax = rbind(mmax, max(oscores[i,scoreName], max(dscores[i,])))
	}
	a= melt(a, id="partitioner")
	names(a)[which(names(a)=="partitioner")] = partitioner
	a$value = as.numeric(a$value)

	oscores[is.na(oscores[,scoreName]),scoreName] = 0
	oscores[is.null(oscores[,scoreName]),scoreName] = 0
	allHists <- ggplot(a, aes(x=value)) + geom_histogram() + xlab(paste0(scoreName," scores")) + geom_vline(data=oscores, aes_string(xintercept=scoreName), linetype=3, col="red") + geom_text(data=oscores, aes_string(x=scoreName,y=0,label=paste0("round(",scoreName,",3)")), col="red", size=3) #+ ggtitle(paste0("Partitioned by: ",partitioner)) 
	subHists <- allHists + facet_grid(formula)
	subHists
	
	ht <- 2*length(levels(d[[partitioner]]))
	g <-arrangeGrob(subPlots, subHists, ncol=2, widths=c(4,4), heights=c(ht,ht))
	ggsave(paste0(path,folder,"/",scoreName,"/", varpath, "/", rankValue, "-", partitioner, ".pdf"), g, width=8,height=ht)	
	# "-", variables[1], "-", variables[2]
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
computeAllScores <- function(dall, d, partitioner, variables){
	scNames = c("msr")#, "jent", "msr" )
	NumSamples <- 1000
	varpath = paste(variables[1],variables[2],sep="-")
	for (scoreName in scNames){
		cat("Computing score: ", scoreName, "\n")
		#---entropy related prep
		r1 <- range(dall[[variables[1]]], na.rm=TRUE)
		r2<- range(dall[[variables[2]]], na.rm=TRUE)
		N <- dim(dall)[1]
		#---Get known partitioner split sizes
		tab <- table(d[[partitioner]], useNA="ifany")
		sizes <- sapply(tab, function(x) x[1])
		cat("SIZES: ",sizes, "\n")
		
		#---split by known partitioner and save scores
		d[[partitioner]] <-factor(d[[partitioner]],exclude=NULL)
		allScores <- getPartitionScores(d, partitioner, variables, r1, r2, N, scoreName)
		if (scoreName == "msr"){
			vs = c( "mean", "stdev", "r2","pval","slope") 
		}else if (scoreName == "scag"){
			vs = c("Outlying"  ,"Skewed"  ,  "Clumpy" ,   "Sparse"  ,  "Striated" , "Convex" ,   "Skinny","Stringy" ,  "Monotonic","AlphaArea","MstLength", "HullPerimeter")
		}else{
			vs = c("jent")
		} 
		
		#---get scores from simulations
		scores <-list()		#---each row is one split's worth of sim scores
		#---get scores aligned
		allScores = allScores[match(names(sizes), allScores[[partitioner]]),]
		for (i in 1:NumSamples){
			eg <- assignToRandomKGroupsBySize(d, sizes)
			eg <- eg[,-which(names(eg)== partitioner)]
			names(eg)[which(names(eg)=="as.factor(gnum)")] <- "randCluster"
			partScores <- getPartitionScores(d=eg, partitioner="randCluster", variables, r1, r2, N, scoreName)
			for (v in vs){
				scores[[v]] <- cbind(scores[[v]], partScores[[v]])
			}
		}
		cat("Saving histograms...\n")
		zscores = list()
		#---loop through the score histograms
		for (v in vs){
			if (! file.exists(paste(path,folder,v,sep="/")))
				dir.create(file.path(path, folder, v), showWarnings = FALSE)
			if (! file.exists(paste(path, folder, v, varpath,sep="/")))
				dir.create(file.path(path, folder, v, varpath), showWarnings = FALSE)
			mu = apply(scores[[v]],1, mean, na.rm=TRUE)
			sdev = apply(scores[[v]],1, sd, na.rm=TRUE)
			if (v == "pval")
				zscores[[v]] <- mean(abs(allScores[[v]] - mu)/sdev, na.rm=TRUE)
			else
				zscores[[v]] <- max(abs(allScores[[v]] - mu)/sdev, na.rm=TRUE)
			cat("VAR: ",partitioner,"\nv: " , v, " mean: " ,mu, " sdev ", sdev, " zscore: " , zscores[[v]], "  ", allScores[[v]],"\n")
			#---save partition + ranking
			pdata <- allScores[,c(partitioner,v)]
			pdata <- cbind(pdata, variable = rep(v, dim(allScores)[1]))
			pdata$variable <- as.factor(pdata$variable)
			names(pdata)[which(names(pdata) == v)] <- "value"
			savePartitionData(d, partitioner, variables, v, pdata, zscores[[v]], scores[[v]], allScores[c(partitioner,v)], varpath)
		}
	}	#---for scoreName
}


#######################################
#######################################
#---where to store data
path <- "/Users/aanand/Documents/RESEARCH Notes/RobustScagnostics/"
datafiles <- c("ourworld.csv")	# "IPEDS_data.csv")
#,"ourworld.csv"
#, "genClustered.csv", "genDonut.csv", "genMonotonic.csv", "genQuadratic.csv", "genSerial.csv", "genStriatedChopped.csv", "genStriatedLines.csv", "adult.csv", "baseball2004.csv", "BOSTON.csv", "Sample - Superstore Sales.csv")

#library(MASS)
#data(birthwt)
#dall=birthwt
#df="birthwt"

for (df in datafiles){
	fname <- paste0(path, "data/", df)
	dall <- read.csv(fname, header=TRUE)

	folder <- strsplit(x=df,"\\.")[[1]][1]
	dir.create(file.path(path, folder), showWarnings = FALSE)
	setwd(file.path(path, folder))

	#---make binned fields factors
	bs=ifelse(grepl("(bin)",names(dall)), TRUE, ifelse(grepl("(cluster)",names(dall)), TRUE, FALSE))
	if (length(bs[which(bs==TRUE)]) > 0)
		dall[[names(dall)[bs]]] <- as.factor(dall[[names(dall)[bs]]])
	#---get list of partitioners
	dfactors <- sapply(dall, is.factor)
	#---for birthwt
	#dfactors <- sapply(dall, function(x) ifelse(length(unique(x)) <10, TRUE, FALSE))
	ps <- names(dall[dfactors])
	partitioners <- vector()
	for (i in 1:length(ps)){
		dall[[ps[i]]] = as.factor(dall[[ps[i]]])
		splits <- split(dall,dall[[ps[i]]])
		avgSplitSize <- mean(sapply(splits, function(x) dim(x)[1])) 
		prop <- dim(dall)[1]/5
		#---for factors that have small cardinality
		#---and split the dataset into relatively big chunks
		if (length(levels(dall[,ps[i]])) < dim(dall)[1]/2 & avgSplitSize >= prop)
			partitioners <- c(partitioners, ps[i])
	}
	allVars <- names(dall)[!(names(dall) %in% ps)]
	
	#---for all pairs of vars and partitioners
	#allVars <- c("DEATH_RT", "BIRTH_RT")
	#allVars <- allVars[-c(1,2,3,4,5,6,8,10,13,14,17,19,21,24,30,31)]
	#allVars <- c("Percent.admitted...total" ,"Graduation.rate...Bachelor.degree.within.6.years..total")
	#allVars <- c("PutOutRate", "AssistRate")
	for (i in 1:length(allVars)){
		for (j in (i+1):length(allVars)){
			if (j > length(allVars))
				break
			for (partitioner in partitioners){
				if (length(levels(dall[[partitioner]])) <= 1){
					next
				}
				variables = c(allVars[i], allVars[j])
				cat("Processing: ", variables, " and partitioner: ", partitioner, "\n")
				d <- dall[,c(variables,partitioner)]
				computeAllScores(dall, d, partitioner, variables)
			}
		}
	}
}
