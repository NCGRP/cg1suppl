#install.packages("colorRamps")
#install.packages("plot3D")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Populus/Experiment2/+revision\ with\ hclust/+results/hclust/+PopGenomicGeography")
library(reshape2)
library(plot3D)

#SUBROUTINES
shapematrix <- function(g,meas,currchr)
{
	g1 <- subset(g, chr==currchr) #extract only the values that correspond to the current chromosome from the matrix g
	#g2 <- g1[,c("blocklength","chrindex", aatype)] #extract the relevant columns to plot with sequential block position (1,2,3,4,...)
	locpos <- g1[,c("locus2")]
	lp2 <- gsub(paste(currchr,"_",sep=""),"",locpos) #remove underscore separator in locus name
	g1$locus2 <- as.numeric(lp2) #assign modified test to data frame as a number
	g2 <- g1[,c("blocklength","locus2", meas)] #extract the relevant columns to plot with genomic location
	
	return(g2)
}


fingerprintplot <- function(g2,infilename,meas,currchr)
{
	t=gsub("\\+v2STATS", "Fplot", infilename)
	pdf(paste(t,meas,"Chr",currchr,".pdf",sep=""))

	#determine breakpoints for color ramp so that zero is white
	ns<-1000 #total number of shades in color ramp
	zpt<-0 #the zero point
	rbwPal <- colorRampPalette(c('blue','white','red'))
	#test whether lower limit is further than upper limit from zero point, to set equal dimensions around zpt
	if ( abs(max(g2$sum)-zpt) >= abs(min(g2$sum)-zpt) ) {
		maxliml<-max(g2$sum)
		minliml<-zpt-max(g2$sum)
	} else {
		maxliml<-zpt+abs(min(g2$sum))
		minliml<-min(g2$sum)
	}
	breakpoints <-seq(from=minliml, to=maxliml, length.out=ns+1)

	#add a color column
	#g2$Colorrbw <- rbwPal(10)[as.numeric(cut(g2$sum,breaks = 10))]  #add column with colors to data
	g2$Colorrbw <- rbwPal(ns)[as.numeric(cut(g2$sum,breaks = breakpoints))]  #add column with colors to data
		
	#plot
	plot(g2$blocklength, g2$locus2, col=g2$Colorrbw, pch=15, main=paste(t,"Chr",currchr,sep=""), xlab="Block length", ylab="Genomic position")
	
	dev.off()
}



#MAIN
#make heatmaps of retention for optimized or random cores
maxchr=19 #the number of chromosomes
infilenames<-c("+v2STATS.Populus.RgeoTgen.txt","+v2STATS.Populus.RpcoTgen.txt")

for (infilename in infilenames)
{
	g <- read.table(infilename, header=TRUE, sep="\t")
	#read so that column with locus identifier is treated as a complete string, and no trailing zeroes are dropped
	
	measures <- c("sum","mean")
	meas = "sum"
	for (currchr in 1:maxchr) #make a separate heatmap for each chromosome
	{
		g2 <- shapematrix(g, meas, currchr)
		fingerprintplot(g2,infilename,meas,currchr) 
	}
}



#library(colorRamps)
#col5 <- colorRampPalette(c('blue', 'gray96', 'red'))  #create color ramp starting from blue to red
#color_levels=20 #the number of colors to use
#max_abolute_value=0.4 #what is the maximum absolute value of raster?
#plot(img, col=col5(n=color_levels), breaks=seq(-max_abolute_value,max_abolute_value,length.out=color_levels+1) , axes=FALSE)
#
#plot(g2$blocklength, g2$locus2, col=rbwPal(n=ns), breaks=seq(from=minliml, to=maxliml, length.out=ns+1), pch=15, main=paste(t,"Chr",currchr,sep=""), xlab="Block length", ylab="Genomic position")
#
#
#
#
#	colorsj <- colorRampPalette(c("blue", "white", "red"))
#	#test whether lower limit is further than upper limit from zero point, to set equal dimensions around zpt
#	if ( abs(max(gm)-zpt) >= abs(min(gm)-zpt) ) {
#		maxliml<-max(gm)
#		minliml<-zpt-max(gm)
#	} else {
#		maxliml<-zpt+abs(min(gm))
#		minliml<-min(gm)
#	}
#	breakpoints <-seq(from=minliml, to=maxliml, length.out=ns+1)
#	image.plot(x=cn, y=rn, z=gm, col=colorsj(ns), breaks=breakpoints, xlab="Haplotype block length", ylab="Core size", main=infilename)
