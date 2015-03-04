rm(list=ls())
library(scagnostics)
library(reshape)
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
getPartitionScores <- function(d, partitioner, variables, r1, r2, N, scoreName){
	switch(scoreName,
		jent={
			groupS <- ddply(d, partitioner, .fun=function(x,col){ e=getEntropyMeasure(x,col,r1,r2, N)}, variables)
		},
		scag={
			groupS <- ddply(d, partitioner, .fun=function(x, col){s= getScagnosticsMeasure(x, col)}, variables)
		},
		msr={
			groupS <- ddply(d, partitioner, .fun=function(x,col){ e=getMSRMeasure(x,col)}, variables)
		},
		stop("Enter something that switches me!")
	)
	return(groupS)
}
#--------
savePartitionData <- function(d, partitioner, variables, scoreName, groupS, rankValue){
	formula <- paste(partitioner, "~ .",sep=" ")
	allPlots <- ggplot(d, aes_string(x=variables[1], y=variables[2])) + geom_point(shape=1) + theme(aspect.ratio = 1)
	subPlots <- allPlots + facet_grid(formula)

	allHists <- ggplot(groupS, aes(x=variable, y=value)) + geom_bar(stat="identity") + geom_text(aes(y=value, ymax=value, label=value), position= position_dodge(width=0.9), vjust=-.5, color="black", size=2)  +ggtitle(paste("Partitioned by: ",partitioner,sep=""))+ theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))  
	#+ scale_y_discrete(limits = c(-1,5))
	subHists <- allHists + facet_grid(formula)

	ht <- 2
	g <-arrangeGrob(allPlots, allHists, ncol=2, widths=c(4,4), heights=c(ht,ht))
	ggsave(paste(path,"simsMaxAbs/",scoreName,"/ALL-", partitioner, "-", variables[1], "-", variables[2], ".pdf", sep=""), g, width=8,height=ht)
	
	ht <- 2*length(levels(d[[partitioner]]))
	g <-arrangeGrob(subPlots, subHists, ncol=2, widths=c(4,4), heights=c(ht,ht))
	ggsave(paste(path,"simsMaxAbs/",scoreName,"/", rankValue, "-", partitioner, "-", variables[1], "-", variables[2], ".pdf", sep=""), g, width=8,height=ht)
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
saveScoreHistograms <- function(v, scores, allScores, partitioner, variables, tab){
	for (i in 1:dim(scores[[v]])[1]){
		#---save histograms for scores from simulations
		pdf(paste(path,"simsMaxAbs/",v,"/", partitioner, "-", variables[1], "-", variables[2],"split-",tab[i],".pdf",sep=""))
		mmin = min(allScores[[v]][i], min(scores[[v]][i,]))
		mmax = max(allScores[[v]][i], max(scores[[v]][i,]))
		hist(scores[[v]][i,], xlim=c(mmin, mmax))
		abline(v=allScores[[v]][i],col="red")
		text(x=c(allScores[[v]][i]),y=c(0),labels=c(round(allScores[[v]][i], digits=3)),col="red")
		dev.off()
	}
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
	scNames = c("scag")#, "jent", "msr" )
	NumSamples <- 100
	
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
		d[[partitioner]] <- as.factor(d[[partitioner]])
		allScores <- getPartitionScores(d, partitioner, variables, r1, r2, N, scoreName)
		if (scoreName == "msr"){
			vs = c( "mean", "stdev", "r2") 
		}else if (scoreName == "scag"){
			vs = c("Outlying"  ,"Skewed"  ,  "Clumpy" ,   "Sparse"  ,  "Striated" , "Convex" ,   "Skinny","Stringy" ,  "Monotonic","AlphaArea","MstLength", "HullPerimeter")
		}else{
			vs = c("jent")
		} 
		
		#---get scores from simulations
		scores <-list()		#---each row is one split's worth of sim scores
		#---get scores aligned
		#splitNames=NULL
		on=names(sort(sizes))
		allScores = allScores[match(on, allScores[[partitioner]]),]
		for (i in 1:NumSamples){
			eg <- assignToRandomKGroupsBySize(d, sizes)
			eg <- eg[,-which(names(eg)== partitioner)]
			names(eg)[which(names(eg)=="as.factor(gnum)")] <- "randCluster"
			partScores <- getPartitionScores(eg, "randCluster", variables, r1, r2, N, scoreName)
			#plot(eg[which(eg$randCluster==6),1:2])
			#partScores
			#splitNames = partScores$randCluster
			for (v in vs){
				scores[[v]] <- cbind(scores[[v]], partScores[[v]])
			}
		}
		cat("Saving histograms...\n")
		zscores = list()
		#---loop through the score histograms
		for (v in vs){
			saveScoreHistograms(v, scores, allScores, partitioner, variables, allScores[[partitioner]])
			mu = apply(scores[[v]],1, mean, na.rm=TRUE)
			sdev = apply(scores[[v]],1, sd, na.rm=TRUE)
			zscores[[v]] <- max(abs(allScores[[v]] - mu)/sdev, na.rm=TRUE)
			cat("VAR: ",partitioner,"\nv: " , v, " mean: " ,mu, " sdev ", sdev, " zscore: " , zscores[[v]],"\n")
			#---save partition + ranking
			pdata <- allScores[,c(partitioner,v)]
			pdata <- cbind(pdata, variable = rep(v, dim(allScores)[1]))
			pdata$variable <- as.factor(pdata$variable)
			names(pdata)[which(names(pdata) == v)] <- "value"
			savePartitionData(d, partitioner, variables, v, pdata, zscores[[v]])
			
			x <- scores[[v]]
			kde <- computeKDE(dim(x)[2], dim(x)[1], x, allScores[[v]])
		}
	}	#---for scoreName
}


#######################################
#######################################
#---where to store data
path <- "/Users/aanand/Documents/RESEARCH Notes/RobustScagnostics/"
datafiles <- c("ourworld.csv", "IPEDS_data.csv", "genClustered.csv", "genDonut.csv", "genMonotonic.csv", "genQuadratic.csv", "genSerial.csv", "genStriatedChopped.csv", "genStriatedLines.csv", "adult.csv")

datafiles <- c("IPEDS_data.csv")
for (df in datafiles){
	fname <- paste(path, "data/", df, sep="")
	dall <- read.csv(fname, header=TRUE)

	#---make binned fields factors
	bs=ifelse(grepl("(bin)",names(dall)), TRUE, ifelse(grepl("(cluster)",names(dall)), TRUE, FALSE))
	if (length(bs[which(bs==TRUE)]) > 0)
		dall[[names(dall)[bs]]] <- as.factor(dall[[names(dall)[bs]]])
	#---get list of partitioners
	dfactors <- sapply(dall, is.factor)
	ps <- names(dall[dfactors])
	partitioners <- vector()
	for (i in 1:length(ps)){
		splits <- split(dall,dall[[ps[i]]])
		avgSplitSize <- mean(sapply(splits, function(x) dim(x)[1])) 
		prop <- dim(dall)[1]/5
		#---for factors that have small cardinality
		#---and split the dataset into relatively big chunks
		if (length(levels(dall[,ps[i]])) < dim(dall)[1]/2 & avgSplitSize >= prop)
			partitioners <- c(partitioners, ps[i])
	}
	allVars <- names(dall)[!(names(dall) %in% partitioners)]
	
	#---for all pairs of vars and partitioners
	#allVars <- c("DEATH_RT", "BIRTH_RT")
	allVars <- c("Percent.admitted...total" ,"Graduation.rate...Bachelor.degree.within.6.years..total")
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
