# Experimental analysis and figure generation code for:

# Automatic Selection of Conditioning Dimensions for Trellis Displays
# Anushka Anand and Justin Talbot
# Tableau Research
# InfoVis 2015

# To regenerate the figures from the paper, with R installed, run:
# Rscript analysis.R


rm(list=ls())
library(scagnostics)
library(entropy)
library(reshape2)
library(plyr)
library(ggplot2)
library(gridExtra)

set.seed(1357)

#--------
standard_theme <- function() {
    theme_set(theme_bw(base_size = 8))
    theme_update(
    panel.background= element_blank(),
    panel.grid.major.y=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.border=element_rect(color="#DDDDDD", fill=NA),
    plot.title=element_text( family="Helvetica", size=11, vjust=1),
    axis.title.x= element_text(vjust=-0.25, family="Helvetica", size=8),
    axis.text.x= element_text( family= "Helvetica", size=7),
    axis.title.y= element_text( family="Helvetica", size=8, angle=90, vjust=0.2),
    axis.text.y= element_text( family= "Helvetica", hjust=1, size=7),
    strip.text.x= element_text( family= "Helvetica", size=8),
    strip.text.y= element_text( family= "Helvetica", size=8),
    legend.key = element_rect(color="#FFFFFF", fill="#FFFFFF"),
    legend.text = element_text( family= "Helvetica", size=9),
    legend.title = element_text( color="#FFFFFF", family= "Helvetica", size=9),
    strip.background= element_rect(color="#DDDDDD", fill="#F3F3F3"),
    plot.margin = unit(c(0.01,0.01,0.1,-0.1), "in")
    )
}
#--------
standard_theme_no_grid <- function() {
    theme_set(theme_bw(base_size = 8))
    theme_update(
    panel.background= element_blank(),
    panel.grid.major.y=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.border=element_rect(color="#DDDDDD", fill=NA),
    plot.title=element_text( family="Helvetica", size=11, vjust=1),
    axis.title.x= element_text(vjust=-0.25, family="Helvetica", size=8),
    axis.text.x= element_text( family= "Helvetica", size=7),
    axis.title.y= element_text( family="Helvetica", size=8, angle=90, vjust=0.2),
    axis.text.y= element_text( family= "Helvetica", hjust=1, size=7),
    strip.text.x= element_text( family= "Helvetica", size=8),
    strip.text.y= element_text( family= "Helvetica", size=8),
    legend.key = element_rect(color="#FFFFFF", fill="#FFFFFF"),
    legend.text = element_text( family= "Helvetica", size=9),
    legend.title = element_text( color="#FFFFFF", family= "Helvetica", size=9),
    strip.background= element_rect(color="#DDDDDD", fill="#F3F3F3"),
    plot.margin = unit(c(0.01,0.01,0.1,-0.1), "in")
    )
}

#--------
getMSRMeasure <- function(d, variables, partitioner){
	tmp = d[!names(d)==partitioner]
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
	return(c(R2=r, Pvalue=p, Slope=sl))
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
savePartitionData <- function(d, partitioner, variables, scoreName, rankValue, dscores, oscores, varpath){
	formula <- paste(partitioner, "~ .",sep=" ")
	allPlots <- ggplot(d, aes_string(x=variables[1], y=variables[2])) + geom_point(shape=1) + theme(aspect.ratio = 1) + theme(axis.title=element_text(size=5)) + theme(axis.text=element_text(size=4))
	if (scoreName %in% c("R2","Pvalue","Slope"))
		allPlots <- allPlots  + geom_smooth(method='lm') 
	standard_theme()
	ggsave(paste0(outputPath,folder,"/", varpath, ".pdf"), allPlots, width=2,height=2)	

	subPlots <- allPlots + facet_grid(formula) + theme(plot.margin= unit(c(0, -1, 0.5, -1), "lines"))	+ theme(strip.background = element_blank(), strip.text.y = element_blank())
	#margin = top, right, bottom, left

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
	cnts=ddply(a, partitioner, summarise,N = max(hist(value, plot=F)$counts))
	maxCnt = max(cnts$N)

	allHists <- ggplot(a, aes(x=value)) + geom_histogram() + xlab(paste0("Distribution of ", scoreName," Measure")) + geom_vline(data=oscores, aes_string(xintercept=scoreName), linetype=3, col="blue") + geom_text(data=oscores, aes_string(x=scoreName,y= maxCnt,label=paste0("round(",scoreName,",3)")), col="blue", size=3) + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())
	subHists <- allHists + facet_grid(formula) + theme(plot.margin= unit(c(0, 0, 0.5, 0), "lines")) 
	
	ht <- 2*length(levels(d[[partitioner]]))
	g <-arrangeGrob(subPlots, subHists, ncol=2, widths=c(2,3), heights=c(ht,ht))
	ggsave(paste0(outputPath,folder,"/",scoreName,"/", rankValue, "-", partitioner, ".pdf"), g, width=6,height=ht)	
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
	scNames = c("msr", "jent", "scag" )
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
			vs = c("R2","Pvalue","Slope") 
		} else if (scoreName == "scag"){
			vs = c("Outlying"  ,"Skewed"  ,  "Clumpy" ,   "Sparse"  ,  "Striated" , "Convex" ,   "Skinny","Stringy" ,  "Monotonic","AlphaArea","MstLength", "HullPerimeter")
		} else{
			vs = c("jent")
		} 
		#---get scores from simulations
		#---each row is one split's worth of sim scores
		scores <-list()
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
		zscores = list()
		#---loop through the score histograms
		for (v in vs){
			if (! file.exists(paste(outputPath,folder,v,sep="/")))
				dir.create(file.path(outputPath, folder, v), showWarnings = FALSE)
			mu = apply(scores[[v]],1, mean, na.rm=TRUE)
			sdev = apply(scores[[v]],1, sd, na.rm=TRUE)
			if (v == "pval")
				zscores[[v]] <- mean(abs(allScores[[v]] - mu)/sdev, na.rm=TRUE)
			else
				zscores[[v]] <- max(abs(allScores[[v]] - mu)/sdev, na.rm=TRUE)
			cat("VAR: ",partitioner,"\nv: " , v, " mean: " ,mu, " sdev ", sdev, " zscore: " , zscores[[v]], "  ", allScores[[v]],"\n")
			#---save partition + ranking
			savePartitionData(d, partitioner, variables, v, zscores[[v]], scores[[v]], allScores[c(partitioner,v)], varpath)
		}
	}	#---for scoreName
}


#######################################
#######################################
#---where to store data
setwd("/Users/aanand/Documents/RESEARCH Notes/RobustScagnostics/")
outputPath <- "output/"
inputPath <- "input/"
#---set up interesting bivariate relationships for datasets
vars <- list()
vars[["ourworld.csv"]] <-c("DEATH_RT", "BIRTH_RT")
vars[["IPEDS_data.csv"]] <- c("Percent.admitted...total", "Graduation.rate...Bachelor.degree.within.6.years..total")
vars[["genDonut.csv"]] <-c("donut1", "donut2")
vars[["genDonut.csv"]] <-c("capital-gain", "fnlwgt")
vars[["Sample - Superstore Sales.csv"]] <- c("Profit","Shipping_Cost") 
datafiles = names(vars)[2]

for (df in datafiles){
	fname <- paste0(inputPath, df)
	dall <- read.csv(fname, header=TRUE)
	folder <- strsplit(x=df,"\\.")[[1]][1]
	dir.create(file.path(outputPath, folder), showWarnings = FALSE)

	#---make binned fields factors
	bs=ifelse(grepl("(bin)",names(dall)), TRUE, ifelse(grepl("(cluster)",names(dall)), TRUE, FALSE))
	if (length(bs[which(bs==TRUE)]) > 0)
		dall[[names(dall)[bs]]] <- as.factor(dall[[names(dall)[bs]]])
	#---get list of partitioners
	dfactors <- sapply(dall, is.factor)
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
	#allVars <- names(dall)[!(names(dall) %in% ps)]
	allVars <- vars[[df]]
	
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
