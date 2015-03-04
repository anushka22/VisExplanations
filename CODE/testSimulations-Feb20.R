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
getAlphaMSTMeasure <- function(d, variables){
	#---only handles non-duplicate data rows
	data <- unique(d[,variables])
	#---doesn't handle missing data or NAs
	data <- data[complete.cases(data),]
	if (dim(data)[1] <= 1){
		return(c( alen=0, ml=0))		#aarea=0,
	}
	#5 for LEADER, 0.5 for clus, 0.2 for donut
	#---scagnostics picked alpha = avg(mst edge length)
	allpairs = dist(data[,variables])
	m = mst(allpairs)
	minds = (m == 1)
	allpairs = as.matrix(allpairs)
	sml = sum(allpairs[minds])
	aml = mean(allpairs[minds])	
#	ah = ahull(data, alpha=aml)	
	as = ashape(data, alpha=aml)	
#	ar = areaahull(ah)
	return(c(alen=as$length, ml=sml))	#aarea=ar, 
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
	scagLabels <- c("Outlying"  ,"Skewed"  ,  "Clumpy" ,   "Sparse"  ,  "Striated" , "Convex" ,   "Skinny","Stringy" ,  "Monotonic")
	result <- tryCatch({
		res <- scagnostics(d[,variables])
		rr <- NULL
		for (si in 1:length(res))
			rr <- cbind(rr, res[si])
		rr <- as.data.frame(rr)
		colnames(rr) <- scagLabels
		return(as.data.frame(rr))
	}, error=function(err){
		rr <- t(as.data.frame(rep(0,9)))
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
			groupS <- ddply(d, partitioner, .fun=function(x, col=variables){s= getScagnosticsMeasure(x, col)})
		},
		alpha={
			groupS <- ddply(d, partitioner, .fun=function(x,col){ e=getAlphaMSTMeasure(x,col)}, variables)
		},
		stop("Enter something that switches me!")
	)
	return(groupS)
}
#--------
savePartitionData <- function(d, partitioner, variables, scoreName, groupS){
	formula <- paste(partitioner, "~ .",sep=" ")
	allPlots <- ggplot(d, aes_string(x=variables[1], y=variables[2])) + geom_point(shape=1) + theme(aspect.ratio = 1)
	subPlots <- allPlots + facet_grid(formula)

	allHists <- ggplot(groupS, aes(x=variable, y=value)) + geom_bar(stat="identity") + geom_text(aes(y=value, ymax=value, label=value), position= position_dodge(width=0.9), vjust=-.5, color="black", size=2)  +ggtitle(paste("Partitioned by: ",partitioner,sep=""))+ theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))  
	#+ scale_y_discrete(limits = c(-1,5))
	subHists <- allHists + facet_grid(formula)

	ht <- 2
	g <-arrangeGrob(allPlots, allHists, ncol=2, widths=c(4,4), heights=c(ht,ht))
	ggsave(paste(path,"sims/",scoreName,"/ALL-", partitioner, "-", variables[1], "-", variables[2], ".pdf", sep=""), g, width=8,height=ht)
	
	ht <- 2*length(levels(d[[partitioner]]))
	g <-arrangeGrob(subPlots, subHists, ncol=2, widths=c(4,4), heights=c(ht,ht))
	ggsave(paste(path,"sims/",scoreName,"/", partitioner, "-", variables[1], "-", variables[2], ".pdf", sep=""), g, width=8,height=ht)
}
#----------------
saveSeparatePartitionScores <- function(d, partitioner, vs, allS){
	for (v in vs){
		pdata <- allS[,c(partitioner,v)]
		pdata <- cbind(pdata, variable = rep(v, dim(allS)[1]))
		pdata$variable <- as.factor(pdata$variable)
		names(pdata)[which(names(pdata) == v)] <- "value"
		savePartitionData(d, partitioner, variables, v, pdata)
	}
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
computeAllScores <- function(dall, d, partitioner, variables){
	scNames = c("jent")#,"alpha", "scag" )
	NumSamples <- 50
	
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
		allS <- getPartitionScores(d, partitioner, variables, r1, r2, N, scoreName)
		if (scoreName == "alpha"){
			vs = c( "alen", "ml") #"aarea",
		}else if (scoreName == "scag"){
			vs = c("Outlying"  ,"Skewed"  ,  "Clumpy" ,   "Sparse"  ,  "Striated" , "Convex" ,   "Skinny","Stringy" ,  "Monotonic")
		}else{
			vs = c("jent")
		} 
		saveSeparatePartitionScores(d, partitioner, vs, allS)
		allScores <- allS
	
		#---get scores from simulations
		scores <-list()		#---each row is one split's worth of sim scores
		for (i in 1:NumSamples){
			eg <- assignToRandomKGroupsBySize(d, sizes)
			eg <- eg[,-which(names(eg)== partitioner)]
			names(eg)[which(names(eg)=="as.factor(gnum)")] <- "randCluster"
			partScores <- getPartitionScores(eg, "randCluster", variables, r1, r2, N, scoreName)
			for (v in vs){
				scores[[v]] <- cbind(scores[[v]], partScores[[v]])
			}
		}
		cat("Saving histograms...\n")
		#---save histograms for scores from simulations
		for (v in vs){
			for (i in 1:dim(scores[[v]])[1]){
				pdf(paste(path,"sims/",v,"/", partitioner, "-", variables[1], "-", variables[2],"split-",sizes[[i]],".pdf",sep=""))
				mmin = min(allScores[[v]][i], min(scores[[v]][i,]))
				mmax = max(allScores[[v]][i], max(scores[[v]][i,]))
				hist(scores[[v]][i,], xlim=c(mmin, mmax))
				abline(v=allScores[[v]][i],col="red")
				text(x=c(allScores[[v]][i]),y=c(0),labels=c(round(allScores[[v]][i], digits=3)),col="red")
				dev.off()
			}
		}
	}	#---for scoreName
}

#######################################
#######################################
fname <- "/Users/aanand/Documents/RESEARCH Notes/RobustScagnostics/data/scagGenData.csv"
dall <- read.table(fname, sep=",", header=TRUE)
names(dall)
partitioner <- "cluster"

#----------------
variables <- c("clus1" ,"clus2")
d <- dall[,variables]
kfit <- kmeans(d, centers=3)
d <- cbind(d, kfit$cluster)
names(d)[which(names(d)=="kfit$cluster")] <- partitioner

#----------------
variables <- c("ser1" ,"ser2")
d <- dall[,variables]
cluster <- ifelse(d$ser2 <0, 0, 1)
d <- cbind(d, cluster)

#----------------
variables <- c("stri1" ,"stri2")
d <- dall[,variables]
cluster <- ifelse(d$stri2 < 0.5, 0, 1)
d <- cbind(d, cluster)

#----------------
variables <- c("stri1" ,"stri2")
d <- dall[,variables]
cluster <- ifelse(d$stri1 < 0.5, 0, ifelse(d$stri1 < 1.5, 1, ifelse(d$stri1 < 2.5, 2, 3)))
d <- cbind(d, cluster)

#----------------
variables <- c("mono1" ,"mono2")
d <- dall[,variables]
cluster <- ifelse(d$mono1 < 0.4, 0, ifelse(d$mono1 < 0.8, 1, 2))
d <- cbind(d, cluster)

#----------------
variables <- c("quad1" ,"quad2")
d <- dall[,variables]
cluster <- ifelse(d$quad1 < 0, 0, 1)
d <- cbind(d, cluster)

#----------------
variables <- c("quad1" ,"quad2")
d <- dall[,variables]
cluster <- ifelse(d$quad1 < 0, 0, 1)
d <- cbind(d, cluster)

#----------------
variables <- c("donut1" ,"donut2")
d <- dall[,variables]
cluster <- ifelse(d$donut1 >= -0.6 & d$donut1 <= 0.6 & d$donut2 >= -0.6 & d$donut2 <= 0.6, 1, 2)
d <- cbind(d, cluster)

#######################################
#######################################
#---where to store data
path <- "/Users/aanand/Documents/RESEARCH Notes/RobustScagnostics/"
datafiles <- c("ourworld.csv", "IPEDS_data.csv", "scagGenData.csv")

datafiles <- c("IPEDS_data.csv")
for (df in datafiles){
	fname <- paste(path, "data/", df, sep="")
	dall <- read.csv(fname, header=TRUE)

	#---make binned fields factors
	bs=ifelse(grepl("(bin)",names(dall)), TRUE, FALSE)
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

