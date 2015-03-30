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
library(lattice)
library(reshape2)
library(plyr)
library(ggplot2)
library(gridExtra)
library(gtools)

set.seed(1357)

#--------
standard_theme <- function() {
    theme_set(theme_bw(base_size = 8))
    theme_update(
    panel.background= element_blank(),
#    panel.grid.major.y=element_blank(),
 #   panel.grid.major.x=element_blank(),
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
	corr = cor(d[,variables], use = "na.or.complete", method="kendall")
	return(c(Rsquared=r, Pvalue=p, Slope=sl, Correlation=corr[1,2]))
}

#--------
getEntropyMeasure <- function(d, variables, r1, r2, N){
	y2d = discretize2d(d[[variables[1]]], d[[variables[2]]], numBins1=5, numBins2=5, r1=r1, r2=r2) 
	freqs = y2d/N
	H = -sum(ifelse(freqs > 0, freqs * log(freqs), 0))
	return(c(Jent=H))
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
getShingledPartitionScores <- function(d, pIndices, pLabels, variables, r1, r2, N, scoreName, partitioner){
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
				groupS <- rbind(groupS, getMSRMeasure(sset,variables, partitioner))
			},
			stop("Enter something that switches me!")
		)
	}
	rownames(groupS) = pLabels
	return(groupS)
}
#----------------
#--------
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
  #bottom, left, top, right
#  par(mar= c(0,0,0,0))
  
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col, width=1.75, height=1.75))
    }
  }
}

#--------
saveTeaserImage <- function(d, partitioner, variables){
	d <- d[complete.cases(d),]
	d <- cbind(d, dummy=rep(1,nrow(d)))
	allPlots <- ggplot(d, aes_string(x=variables[1], y=variables[2])) + geom_point(shape=1, size=1, alpha=0.5) +  theme(aspect.ratio = 1, plot.margin= unit(c(0, 0.5, 1, 0), "lines"), axis.title=element_text(size=5), axis.text=element_text(size=4), plot.title = element_text(size = 8, vjust = -20)) + coord_cartesian(xlim=c(0, 110), ylim=c(0, 110)) + scale_y_continuous(breaks=c(0, 50, 100), labels=c("0%", "50%", "100%"), name="Graduation Rate") + scale_x_continuous(breaks=c(0, 50, 100), labels=c("0%", "50%", "100%"), name="Admission Rate")
	standard_theme()
	formula <- "~ dummy"
	subPlots1 <- allPlots + facet_grid(formula) + theme(plot.margin= unit(c(0.5, 0.5, 1, 0), "lines"), plot.title = element_text(size=8, vjust=-20), panel.margin = unit(0.75, "lines"), strip.background=element_rect(fill='white', color='white'), strip.text.x = element_text(color='white')) 	
	ggsave(paste0(outputPath,folder,"/teaser1.pdf"), subPlots1, width=2, height=2)	

	levels(d[[partitioner]]) <- c("10-20", "20-30", "30-36")
	formula <- paste("~ ",partitioner,sep=" ")
	subPlots <- allPlots + facet_grid(formula) + theme(plot.margin= unit(c(0.5, 0.5, 1, 0), "lines"), plot.title = element_text(size=8, vjust=-20), panel.margin = unit(0.75, "lines")) 
#	g <-arrangeGrob(allPlots, subPlots, ncol=2, widths=c(1.75,4.5), heights=c(1.5,2))	
	ggsave(paste0(outputPath,folder,"/teaser2.pdf"), subPlots, height=2, width=5)
}
#--------
saveShingledPartitionData <- function(d, partitioner, variables, varpath, scoreName,rankValue, dscores, oscores, pIndices, pLabels){
	
	allPlots <- ggplot(d, aes_string(x=variables[1], y=variables[2])) + geom_point(shape=1, size=0.75, alpha=0.5)+ theme(aspect.ratio = 1, axis.title.y=element_text(vjust=1,size=6), axis.text.y=element_text(size=4), axis.text.x=element_text(size=4), axis.title.x=element_text(vjust=0.25, size=6), plot.margin= unit(c(0, 0.5, 0, 0), "lines")) + xlab("Proportion of old homes") + ylab("Median Value") + scale_y_continuous(limits=c(0, 55), breaks=seq(0,50,10)) 
#	allPlots <- allPlots  + geom_smooth(method='lm') 
	standard_theme()
	ggsave(paste0(outputPath,folder,"/", varpath, ".pdf"), allPlots, width=1.3, height=1.3)

	nplots <- dim(pIndices)[1]
#	ht <- 1.8*length(nplots)
	pplots = list()
	for (i in 1:(nplots-1)){
		sset = d[pIndices[i,1]:pIndices[i,2],]
		name = pLabels[i]
		pplots[[name]] <- ggplot(sset, aes_string(x=variables[1], y=variables[2])) + geom_point(shape=1, size=0.75, alpha=0.5) + theme(aspect.ratio = 1, plot.margin = unit(c(0.05,-1,0,-1), "cm"),  axis.title.y=element_text(vjust=1,size=6), axis.text.y=element_text(size=4), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + scale_y_continuous(limits=c(0, 55), breaks=seq(0,50,10)) + ylab("Median Value") + xlab(" ")
		if (scoreName %in% c("Rsquared","Pvalue","Slope","Correlation"))
			pplots[[name]] <- pplots[[name]]  + geom_smooth(method='lm') 
	}
	sset = d[pIndices[nplots,1]:pIndices[nplots,2],]
	name = pLabels[nplots]
	pplots[[name]] <- ggplot(sset, aes_string(x=variables[1], y=variables[2])) + geom_point(shape=1, size=0.75, alpha=0.5) + theme(aspect.ratio = 1, plot.margin = unit(c(0,-1,0,-1), "cm"), axis.title.y=element_text(vjust=1,size=6), axis.text.y=element_text(size=4), axis.text.x=element_text(size=4), axis.title.x=element_text(vjust=0.5, size=6)) + scale_y_continuous(limits=c(0, 55), breaks=seq(0,50,10))+ ylab("Median Value") + xlab("Proportion of old homes")
	if (scoreName %in% c("Rsquared","Pvalue","Slope","Correlation"))
		pplots[[name]] <- pplots[[name]]  + geom_smooth(method='lm') 

	standard_theme()
	ht <- 5	#1.75*length(pLabels)
	lout=matrix(1:nplots, ncol=1, byrow=FALSE)
	fpath=paste0(outputPath,folder,"/",scoreName,"/", partitioner, ".pdf")
#	fpath=paste0(outputPath,folder,"/",scoreName,"/randCluster.pdf")
	pdf(fpath, width=1.5, height=ht)
	multiplot(plotlist=pplots, cols=1,layout=lout)
	dev.off()

	a = data.frame(cbind(dscores, partitioner=as.vector(pLabels)))
	mmin = NULL
	mmax = NULL
	for (i in 1:nrow(dscores)){
		mmin = rbind(mmin, min(oscores[i], min(dscores[i,], na.rm=T), na.rm=T))
		mmax = rbind(mmax, max(oscores[i], max(dscores[i,], na.rm=T), na.rm=T))
	}
	a= melt(a, id="partitioner")
	names(a)[which(names(a)=="partitioner")] = partitioner
	a$value = as.numeric(a$value)
	oscores = data.frame(cbind(oscores, partitioner=pLabels))
	names(oscores) = c("values", partitioner)
	oscores$values = as.numeric(levels(oscores$values)[oscores$values])
	oscores[is.na(oscores[,"values"]), "values"] = 0
	oscores[is.null(oscores[,"values"]), "values"] = 0
	maxs=ddply(a, partitioner, summarise,N = max(value))
	notNAs = a[a[[partitioner]] %in% maxs[!is.na(maxs$N),partitioner],]
	cnts=ddply(notNAs, partitioner, summarise,N = mean(hist(value, plot=F)$counts))
	maxCnt = max(cnts$N)+50
	cat("maxCnt: " , maxCnt)
	
	rLabels = gsub(", ", ",\n", pLabels)
	levels(a[[partitioner]]) = rLabels
	levels(oscores[[partitioner]]) = rLabels

	formula =paste(partitioner, "~ .",sep=" ")
	allHists <- ggplot(a, aes(x=value)) + geom_histogram() 
	subHists <- allHists + facet_grid(formula, scales="free")
	subHists <- subHists + geom_vline(data=oscores, aes(xintercept=values), linetype=3, col="blue") + geom_text(data=oscores, aes(x=values,y=maxCnt,label=as.character(round(values,3))), col="blue", size=2) +  theme(strip.text.y=element_text(size=5), strip.text.x=element_text(size=5), axis.text.x=element_text(size=5), axis.title.x=element_text(vjust=0.5, hjust=0.5, size=6), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank(), plot.margin= unit(c(0.5, 0, 0.5, 0), "lines"), panel.margin=unit(1.5,"lines")) + xlab(paste0("Distribution of ",scoreName, " measure"))

	ggsave(paste0(outputPath,folder,"/",scoreName,"/hist-", partitioner, ".pdf"), subHists, width=2.5, height=5)	

	pplots[["hists"]] <- subHists
	lout=matrix(c(1:nplots,rep(nplots+1,nplots)), ncol=2, byrow=FALSE)
	fpath=paste0(outputPath,folder,"/",scoreName,"/",rankValue,"-", partitioner, ".pdf")
	pdf(fpath,width=3.5,height=ht)
	multiplot(plotlist=pplots, cols=2,layout=lout)
	dev.off()
}

#--------
savePartitionData <- function(d, partitioner, variables, scoreName, rankValue, dscores, oscores, varpath, index=""){
	d[[partitioner]] = as.factor(d[[partitioner]])
	d[[partitioner]] = factor(d[[partitioner]], levels = mixedsort(levels(d[[partitioner]])))
	#levels(d[[partitioner]])=gsub("\\.", "\n ", levels(d[[partitioner]]))
#	levels(d[[partitioner]]) = c("Colleges\nArts & Sciences", "Colleges\nDiverse fields", "Associates\nColleges", "Research\nUniversities", "Universities\nLarge", "Universities\nMedium", "Universities\nSmall", "Universities\nHigh Research", "Universities\nV. High Research")
#	d[[partitioner]] = as.factor(d[[partitioner]])
	
	allPlots <- ggplot(d, aes_string(x=variables[1], y=variables[2])) + geom_point(shape=1, size=1, alpha=0.7)+ theme(aspect.ratio = 1, axis.title=element_text(size=6), axis.text=element_text(size=4), axis.title.x=element_text(vjust=0.5), plot.margin= unit(c(0, 0.5, 0, 0), "lines")) + xlab("Death rate") + ylab("Birth rate")
	#+ xlab("Admission rate") + ylab("Graduation rate")

	#---parsimony
	#allPlots <- ggplot(d, aes_string(x=variables[1], y=variables[2])) + geom_point(shape=1, size=0.5, alpha=0.7)+ theme(aspect.ratio = 1, axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.title.x=element_text(vjust=0.5), plot.margin= unit(c(0, 0, 0, 0), "lines"))
	if (scoreName %in% c("Rsquared","Pvalue","Slope","Correlation"))
		allPlots <- allPlots  + geom_smooth(method='lm') 
	standard_theme()
	ggsave(paste0(outputPath,folder,"/", varpath, ".pdf"), allPlots, width=1.75, height=1.75)
#	ggsave(paste0(outputPath,folder,"/", varpath, ".pdf"), allPlots, width=0.8, height=0.8)

#	formula <- paste0( "~ ",partitioner)
	formula <- paste0( partitioner,"~ .")
	subPlots <- allPlots + facet_grid(formula) + geom_point(shape=1, size=1, alpha=0.7)+ theme(plot.margin= unit(c(0.5, 0, 0.5, 0), "lines"), strip.background =element_blank(), strip.text.y = element_blank(), aspect.ratio = 1, axis.title.x=element_text(size=5), axis.text.x=element_text(size=3),  axis.title.y=element_text(size=5), axis.text.y=element_text(size=3)) 
#	subPlots <- allPlots + facet_wrap(as.formula(formula),ncol=2) + geom_point(shape=1, size=1, alpha=0.7)+ theme(plot.margin= unit(c(0, 0.5, 0, 0.5), "lines"), strip.text.y = element_text(size=3), strip.text.y = element_blank(), aspect.ratio = 1, axis.title.x=element_text(size=5), axis.text.x=element_text(size=3),  axis.title.y=element_text(size=5), axis.text.y=element_text(size=3)) 

	#---parsimony
	#subPlots <- allPlots + facet_grid(formula) + geom_point(shape=1, size=0.5, alpha=0.5)+ theme(aspect.ratio = 1, axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.title.x=element_text(vjust=0.5), plot.margin= unit(c(0, 0, 0, -0.75), "lines"), strip.background =element_blank(), strip.text.y = element_blank())
	ht <- 3.5 #1.5*length(levels(d[[partitioner]]))
	ggsave(paste0(outputPath,folder,"/",scoreName,"/", partitioner, ".pdf"), subPlots, width=3.5,height=ht)	

	a = data.frame(cbind(dscores, partitioner=as.character(levels(d[[partitioner]]))))
	mmin = NULL
	mmax = NULL
	for (i in 1:nrow(dscores)){
		mmin = rbind(mmin, min(oscores[i,scoreName], min(dscores[i,])))
		mmax = rbind(mmax, max(oscores[i,scoreName], max(dscores[i,])))
	}
	a= melt(a, id="partitioner")
	names(a)[which(names(a)=="partitioner")] = partitioner
	a$value = as.numeric(a$value)

	oscores[is.na(oscores[,scoreName]),scoreName] = 0
	oscores[is.null(oscores[,scoreName]),scoreName] = 0
	maxs=ddply(a, partitioner, summarise,N = max(value))
	notNAs = a[a[[partitioner]] %in% maxs[!is.na(maxs$N),partitioner],]
	cnts=ddply(notNAs, partitioner, summarise,N = mean(hist(value, plot=F)$counts))
	maxCnt = max(cnts$N)+150

	levels(oscores[[partitioner]])=gsub("\\.", "\n", levels(oscores[[partitioner]]))
	levels(oscores[[partitioner]])=gsub(" ", "\n", levels(oscores[[partitioner]]))
	levels(a[[partitioner]])=gsub("\\.", "\n", levels(a[[partitioner]]))
	levels(a[[partitioner]])=gsub(" ", "\n", levels(a[[partitioner]]))

#	levels(oscores[[partitioner]]) = c("Colleges\nArts & Sciences", "Colleges\nDiverse fields", "Associates\nColleges", "Research\nUniversities", "Universities\nLarge", "Universities\nMedium", "Universities\nSmall", "Universities\nHigh Research", "Universities\nV. High Research")

	allHists <- ggplot(a, aes(x=value)) + geom_histogram() + xlab(paste0("Distribution of ", scoreName," Measure")) + geom_vline(data=oscores, aes_string(xintercept=scoreName), linetype=3, col="blue") + geom_text(data=oscores, aes_string(x=scoreName, y= maxCnt,label=paste0("round(",scoreName,",3)")), col="blue", size=2) + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), strip.text.x = element_text(size=5), strip.text.y = element_text(size=5), axis.title.x=element_text(size=5), axis.text.x=element_text(size=3))
	subHists <- allHists + facet_grid(formula) + theme(plot.margin= unit(c(0.5, 0, 0.5, 0), "lines")) 
	
#	ht=3
	ggsave(paste0(outputPath,folder,"/",scoreName,"/hist-", partitioner, ".pdf"), subHists, width=2.5,height=ht)	
	g <-arrangeGrob(subPlots, subHists, ncol=2, widths=c(1.25,1.75), heights=c(ht,ht))
	ggsave(paste0(outputPath,folder,"/",scoreName,"/",rankValue,"-", partitioner, ".pdf"), g, width=3,height=ht)	

	#---parsimony
#	allHists <- ggplot(a, aes(x=value)) + geom_histogram() + xlab(paste0("Distribution of ", scoreName," Measure")) + geom_vline(data=oscores, aes_string(xintercept=scoreName), linetype=3, col="blue")+ theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank(), axis.title.x = element_text(size=1), strip.text.x = element_text(size=1))
#	subHists <- allHists + facet_grid(formula) + theme(plot.margin= unit(c(0, 0, 0, -1), "lines")) 

	#---cluster
#	g <-arrangeGrob(subPlots, subHists, ncol=2, widths=c(1.5,1.5), heights=c(1.5,1.5))
#	ggsave(paste0(outputPath,folder,"/",scoreName,"/",rankValue,"-", partitioner, ".pdf"), g, width=1.5,height=1.5)
	#---cluster1
#	g <-arrangeGrob(subPlots, subHists, ncol=2, widths=c(1.5,1.5), heights=c(2.5,2.5))
#	ggsave(paste0(outputPath,folder,"/",scoreName,"/",rankValue,"-", partitioner, ".pdf"), g, width=1.5,height=2.5)
	#---cluster2
#	g <-arrangeGrob(subPlots, subHists, ncol=2, widths=c(1.5,1.5), heights=c(4.5,4.5))
#	ggsave(paste0(outputPath,folder,"/",scoreName,"/",rankValue,"-", partitioner, ".pdf"), g, width=1.5,height=4.5)

}
#----------------
#----------------
assignToRandomKGroupsBySize <- function(d, ksizes){
	nrows <- nrow(d)
	td<-d[sample(1:nrows,nrows,replace=FALSE),]
	gnum <- NULL
	for (i in 1:length(ksizes)){
		gnum <- c(gnum, rep(names(ksizes)[i], ksizes[[i]]))
	}
	return(cbind(td, as.factor(gnum)))
}

#----------------
computeShingles <- function(d, partitioner){
	d=d[complete.cases(d),]
	z = equal.count(d[[partitioner]],number=4,overlap=0)
#	plot(z)
	t=levels(z)
	s <- sort(d[[partitioner]])
	inds <- NULL
	labels <- NULL
	for (i in 1:length(t)){
		ss = which(s >= t[[i]][1] & s <= t[[i]][2])
		inds <- rbind(inds, c(ss[1], ss[length(ss)]))
	}
	return(list(inds=inds, labels=as.character(t)))
}
#----------------
#----------------
computeAllShingledScores <- function(d, partitioner, variables, res){
	scNames = c("scag")	#"msr", "jent",
	NumSamples <- 1000
	varpath = paste(variables[1],variables[2],sep="-")
	d<- d[order(d[[partitioner]]),]

	for (scoreName in scNames){
		cat("Computing score: ", scoreName, "\n")
		#---entropy related prep
		r1 <- range(dall[[variables[1]]], na.rm=TRUE)
		r2<- range(dall[[variables[2]]], na.rm=TRUE)
		N <- dim(dall)[1]

		#---split by known partitioner and save scores
		pIndices = res$inds
		pLabels = res$labels
		for (i in 1:length(pLabels)){
			ll = pLabels[i]
			q=as.numeric(unlist(strsplit(unlist(ll), "[^0-9.]+")))
			pLabels[i] = paste0("[ ",round(q[2],1), ",",round(q[3],1)," ]")
		}
		allScores <- getShingledPartitionScores(d, pIndices, pLabels, variables, r1, r2, N, scoreName, partitioner)
		allScores = as.data.frame(allScores)
		#---Get known partitioner split sizes
		sizes <- pIndices[,2]-pIndices[,1] + 1
		cat("SIZES: ",sizes, "\n")
		
		if (scoreName == "msr"){
			vs = c("Rsquared","Pvalue","Slope","Correlation") 
		} else if (scoreName == "scag"){
			vs = c("Outlying"  ,"Skewed"  ,  "Clumpy" ,   "Sparse"  ,  "Striated" , "Convex" ,   "Skinny","Stringy" ,  "Monotonic","AlphaArea","MstLength", "HullPerimeter")
		} else{
			vs = c("Jent")
		} 
		
		#---get scores from simulations
		scores <-list()
		#---get scores aligned
		for (i in 1:NumSamples){
			#---randomly permute
			nrows <- dim(d)[1]
			eg<-d[sample(1:nrows,nrows,replace=FALSE),]
			partScores <- getShingledPartitionScores(eg, pIndices, pLabels, variables, r1, r2, N, scoreName, partitioner)
			partScores = as.data.frame(partScores)
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
			if (v == "Pvalue"){
				zscores[[v]] <- mean(abs(allScores[[v]] - mu)/sdev, na.rm=TRUE)
			}else{
				zscores[[v]] <- max(abs(allScores[[v]] - mu)/sdev, na.rm=TRUE)
			}
			zscores[[v]] = replace(zscores[[v]], is.infinite(zscores[[v]]),NA)
			cat("VAR: ",partitioner,"\nv: " , v, " mean: " ,mu, " sdev ", sdev, " zscore: " , zscores[[v]], "  ", allScores[[v]],"\n")
			#---save partition + ranking
			saveShingledPartitionData(d, partitioner, variables, varpath, v, zscores[[v]], scores[[v]], allScores[[v]], pIndices, pLabels)
		}
	}	#---for scoreName	
}

#----------------
computeAllScores <- function(dall, d, partitioner, variables){
	scNames = c("scag")	#"msr", "jent",
	NumSamples <- 500
	varpath = paste(variables[1],variables[2],sep="-")
	
	for (scoreName in scNames){
		cat("Computing score: ", scoreName, "\n")
		#---entropy related prep
		r1 <- range(dall[[variables[1]]], na.rm=TRUE)
		r2<- range(dall[[variables[2]]], na.rm=TRUE)
		N <- nrow(dall)
		#---Get known partitioner split sizes
		d=d[complete.cases(d),]
#		tmp=as.character(d[[partitioner]])
#		d=d[-which(is.na(tmp)),]
#		t=levels(d[[partitioner]])
#		levels(d[[partitioner]]) <- levels(d[[partitioner]])[is.na(t)==F]
		tab <- table(d[[partitioner]])
		if (0 %in% tab){
			tab <- tab[-which(tab==0)]
		}
		sizes <- sapply(tab, function(x) x[1])
		cat("SIZES: ",sizes, "\n")
		#---split by known partitioner and save scores
		d[[partitioner]] <-factor(d[[partitioner]],exclude=NULL)
		allScores <- getPartitionScores(d, partitioner, variables, r1, r2, N, scoreName)
		if (scoreName == "msr"){
			vs = c("Rsquared","Pvalue","Slope","Correlation") 
		} else if (scoreName == "scag"){
			vs = c("Outlying"  ,"Skewed"  ,  "Clumpy" ,   "Sparse"  ,  "Striated" , "Convex" ,   "Skinny","Stringy" ,  "Monotonic","AlphaArea","MstLength", "HullPerimeter")
		} else{
			vs = c("Jent")
		} 
		#---get scores from simulations
		#---each row is one split's worth of sim scores
		scores <-list()
		#---get scores aligned
		allScores = allScores[match(mixedsort(names(sizes)), allScores[[partitioner]]),]
		for (i in 1:NumSamples){
			eg <- assignToRandomKGroupsBySize(d, sizes)
			eg <- eg[,-which(names(eg)== partitioner)]
			names(eg)[which(names(eg)=="as.factor(gnum)")] <- "randCluster"
			partScores <- getPartitionScores(d=eg, partitioner="randCluster", variables, r1, r2, N, scoreName)
			partScores = partScores[match(mixedsort(names(sizes)), partScores[["randCluster"]]),]
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
			if (v == "Pvalue"){
				zscores[[v]] <- mean(abs(allScores[[v]] - mu)/sdev, na.rm=TRUE)
			}else{
				zscores[[v]] <- max(abs(allScores[[v]] - mu)/sdev, na.rm=TRUE)
			}
			zscores[[v]] = replace(zscores[[v]], is.infinite(zscores[[v]]),NA)
			cat("VAR: ",partitioner,"\nv: " , v, " mean: " ,mu, " sdev ", sdev, " zscore: " , zscores[[v]], "  ", allScores[[v]],"\n")
			#---save partition + ranking
			savePartitionData(d, partitioner, variables, v, zscores[[v]], scores[[v]], allScores[c(partitioner,v)], varpath)
		}
	}	#---for scoreName
}
#--------
#--------
computeTeaser <- function(vars){
	df=names(vars)[2]		#---IPEDS + care about Correlation
	variables = vars[[df]]
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
		cat(ps[[i]], " len ",length(splits), "\n")
		avgSplitSize <- mean(sapply(splits, function(x) nrow(x))) 
		prop <- nrow(dall)/10
		#---for factors that have small cardinality
		#---and split the dataset into relatively big chunks
		if (length(levels(dall[,ps[i]])) < nrow(dall)/2 & avgSplitSize >= prop & length(splits) > 1)
			partitioners <- c(partitioners, ps[i])
	}

	for (partitioner in partitioners){
		d <- dall[,c(variables,partitioner)]
		computeAllScores(dall, d, partitioner, variables)
	}	
}
#--------
computeSupport <- function(vars){
	df=names(vars)[2]		#---IPEDS + care about Correlation
	variables = vars[[df]]
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
		cat(ps[[i]], " len ",length(splits), "\n")
		avgSplitSize <- mean(sapply(splits, function(x) nrow(x))) 
		prop <- nrow(dall)/10
		#---for factors that have small cardinality
		#---and split the dataset into relatively big chunks
		if (length(levels(dall[,ps[i]])) < nrow(dall)/2 & avgSplitSize >= prop & length(splits) > 1)
			partitioners <- c(partitioners, ps[i])
	}

	d <- dall
	varpath = paste(variables[1],variables[2],sep="-")
	setName = "scag"
	scoreName = "Monotonic"
	NumSamples <- 500
	nThrows = 10
	res = data.frame(partitioner=character(), percent=double(), zscore=double(), stringsAsFactors=F)
	ns <- NULL
	for (it in 0:nThrows){
		x = ifelse(it > 0, 0.1, 0)
		#---remove 10% of the rows randomly and progressively
		thrown <- sample(1:nrow(d), x*nrow(dall), replace=FALSE)
		if (x > 0){
			if(x >= 10){
				toThrow=nrow(d)/2
				thrown <- sample(1:nrow(d), x*nrow(dall), replace=FALSE)
			}
			d <- d[-thrown,]
		}
		N <- nrow(d)
		cat("Processing # pts: ", N,"\n")
		#---entropy related prep
		r1 <- range(d[[variables[1]]], na.rm=T, showWarnings=F)
		r2<- range(d[[variables[2]]], na.rm=T, showWarnings=F)
		#---loop through the partitioners
		for (partitioner in partitioners){
			dat <- d[,c(variables, partitioner)]
			dat <- dat[complete.cases(dat),]
			#---Get known partitioner split sizes
			dat[[partitioner]] <-factor(dat[[partitioner]],exclude=NULL)
			tab <- table(dat[[partitioner]])
			sizes <- sapply(tab, function(x) x[1])
			if (length(sizes) < 2){
				next
			}
			ns <- c(ns, N)
			#---nodata =names(sizes[which(sizes==0)])
			cat(partitioner, " SIZES: ",sizes, "\n")
			#---split by known partitioner and save scores
			allScores <- getPartitionScores(dat, partitioner, variables, r1, r2, N, setName)
			#---get scores from simulations
			scores <-list()
			#---get scores aligned
			allScores = allScores[match(names(sizes), allScores[[partitioner]]),]
			for (i in 1:NumSamples){
				eg <- assignToRandomKGroupsBySize(dat, sizes)
				eg <- eg[,-which(names(eg)== partitioner)]
				names(eg)[which(names(eg)=="as.factor(gnum)")] <- "randCluster"
				#---because they are the same factors (partitioner levels), the ddply result will be ordered the same 
				partScores <- getPartitionScores(d=eg, partitioner="randCluster", variables, r1, r2, N, setName)
				scores[[scoreName]] <- cbind(scores[[scoreName]], partScores[[scoreName]])
			}
			#---get zscores
			mu = apply(scores[[scoreName]],1, mean, na.rm=TRUE)
			sdev = apply(scores[[scoreName]],1, sd, na.rm=TRUE)
			zscores <- max(abs(allScores[[scoreName]] - mu)/sdev, na.rm=TRUE)
			zscores = replace(zscores, is.infinite(zscores),NA)
			percent = 1-(x* it)
			cat(partitioner, " zscore: ", zscores, " mu " , mu, " sd ", sdev , " percent ", percent,"\n")
			res <- rbind(res, data.frame(list(partitioner=as.character(partitioner), percent=percent, zscore=zscores),stringsAsFactors=F))
			
			savePartitionData(dat, partitioner, variables, scoreName, zscores, scores[[scoreName]], allScores[c(partitioner,scoreName)], varpath,it)
		}#---for each partitioner
	}#---for each throw
	
	res1 <- cbind(res, ns)
	res2=data.frame(lapply(res1, function(x) replace(x, is.infinite(x),NA)))
	annot = res1[res1$ns == range(res1$ns)[2],]
	sannot = annot[order(-annot$zscore),]
	annot = sannot[1:4,]
	annot$partitioner = as.factor(annot$partitioner)
	levels(annot$partitioner) = c("ACT 75th\nPercentile\nScores", "Carnegie\nClassification", "Historically\nBlack College", "Offers Assoc.\nDegree")
	annot[which(annot$partitioner=="Carnegie\nClassification"),3] =annot[which(annot$partitioner=="Carnegie\nClassification"),3]-0.5
	
	g=ggplot(data=res2, aes(x=ns, y=zscore, group=factor(partitioner), col=factor(partitioner))) + geom_line(size=0.25) + geom_point(size=1) + scale_x_continuous(name="Number of Data Points", breaks=seq(0, 1600, 200), labels=seq(0, 1600, 200)) + coord_cartesian(xlim=c(-5,1850), ylim=c(0, 20)) + scale_y_continuous(breaks=seq(0, 20, 5), name="Z-Scores") + theme(plot.margin= unit(c(0.5, 0, 0.5, 0.5), "lines"), legend.position = "none", plot.title=element_text(size=7), axis.title=element_text(size=4), axis.text=element_text(size=3)) + geom_text(data = annot, aes(x=ns, y=zscore,label = partitioner), hjust = -0.1, vjust = -0.1, size=1.75) + ggtitle("Effect of Support (Number of data points)")
#	standard_theme_no_grid()
	standard_theme()
	ggsave(paste0(outputPath,folder,"/",scoreName,"/support", ".pdf"), g, width=3.5,height=3.5)	
}

#----------
computeInformative <- function(vars){
	df=names(vars)[1]		#---ourworld and Monotonic
	variables = vars[[df]]
	fname <- paste0(inputPath, df)
	dall <- read.csv(fname, header=TRUE)
	folder <- strsplit(x=df,"\\.")[[1]][1]
	dir.create(file.path(outputPath, folder), showWarnings = FALSE)

	dfactors <- sapply(dall, is.factor)
	partitioners <- names(dall[dfactors])
	partitioners=partitioners[-1]
	for (partitioner in partitioners){
		d <- dall[,c(variables,partitioner)]
		computeAllScores(dall, d, partitioner, variables)
	}
}
#----------
computeInformative <- function(vars){
	df=names(vars)[1]		#---ourworld + care about Monotonic
	variables = vars[[df]]
	fname <- paste0(inputPath, df)
	dall <- read.csv(fname, header=TRUE)
	folder <- strsplit(x=df,"\\.")[[1]][1]
	dir.create(file.path(outputPath, folder), showWarnings = FALSE)

	dfactors <- sapply(dall, is.factor)
	partitioners <- names(dall[dfactors])
	for (partitioner in partitioners){
		d <- dall[,c(variables,partitioner)]
		computeAllScores(dall, d, partitioner, variables)
	}
}

#----------
computeVisuallyRich <- function(vars){
	df=names(vars)[7]		#---olive oils + care about Striation
	variables = vars[[df]]
	fname <- paste0(inputPath, df)
	dall <- read.csv(fname, header=TRUE)
	folder <- strsplit(x=df,"\\.")[[1]][1]
	dir.create(file.path(outputPath, folder), showWarnings = FALSE)

	dfactors <- sapply(dall, is.factor)
	partitioners <- names(dall[dfactors])
	for (partitioner in partitioners){
		d <- dall[,c(variables,partitioner)]
		computeAllScores(dall, d, partitioner, variables)
	}
}

#----------
computeParsimonius <- function(vars){
	df=names(vars)[3]	#---genDonut -- care about Clumpy
	variables = vars[[df]]
	fname <- paste0(inputPath, df)
	dall <- read.csv(fname, header=TRUE)
	folder <- strsplit(x=df,"\\.")[[1]][1]
	dir.create(file.path(outputPath, folder), showWarnings = FALSE)

	partitioner="cluster"
	d <- dall[,c(variables,partitioner)]
	d[[partitioner]] = as.factor(d[[partitioner]])
	nparts = levels(d[[partitioner]])
	for (ti in 1:3){
		newpart=paste0("cluster",ti)
		d[[newpart]] <- 0
		for (part in nparts){
			pieces = 2*ti
			p1 = which(d[[partitioner]] == part)
			np = length(p1)
			i=1
			locs <- vector()
			for (pi in 1:(pieces-1)){
				inds <- sample(p1,np/(pieces),replace=FALSE)
				locs <- match(inds, p1)
				d[inds,newpart] = paste0(part, letters[i])
				p1 <- p1[-locs]
				i=i+1
			}
			d[p1,newpart] = paste0(part, letters[i])
		}
		d[[newpart]] = as.factor(d[[newpart]])
	}
	dfactors <- sapply(d, is.factor)
	ps <- names(d[dfactors])
	for (partitioner in ps){
		computeAllScores(dall, d, partitioner, variables)
	}
}

#------------
#---uses shingles too
computeRunningExample <- function(vars){
	df=names(vars)[8]	#---boston + care about Monotonic
	variables = vars[[df]]
	
	fname <- paste0(inputPath, df)
	dall <- read.csv(fname, header=TRUE)
	folder <- strsplit(x=df,"\\.")[[1]][1]
	dir.create(file.path(outputPath, folder), showWarnings = FALSE)
	dall[["RAD"]] = as.factor(dall[["RAD"]])
	dall[["CHAS"]] = as.factor(dall[["CHAS"]])

	dfactors <- sapply(dall, is.factor)
	partitioners <- names(dall[dfactors])
	for (partitioner in partitioners){
		d <- dall[,c(variables,partitioner)]
		computeAllScores(dall, d, partitioner, variables)
	}
	qpartitioners <- names(dall)[-which(names(dall) %in% c(partitioners, variables))]
	for (partitioner in qpartitioners){
		d <- dall[,c(variables,partitioner)]
		res = computeShingles(d, partitioner)
		computeAllShingledScores(d, partitioner, variables, res)
	}
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
vars[["Sample - Superstore Sales.csv"]] <- c("Profit","Shipping_Cost") 
vars[["adult-test.csv"]] <- c("age","capital.gain") 
vars[["iris.csv"]] <- c("petal.width","sepal.length") 
vars[["oliveoil.csv"]] <- c("linolenic","linoleic") 
vars[["boston.csv"]] <- c("AGE", "MEDV") 

computeTeaser(vars)

computeRunningExample(vars)

computeParsimonius(vars)

computeVisuallyRich(vars)

computeSupport(vars)

computeInformative(vars)

}
