options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+ Study/Analysis/Populus/M+ runs/SlidingWindow/heatmap")
library(reshape2)
library(plot3D)

#SUBROUTINES
shapematrix <- function(g,aatype,currchr)
{
	g1 <- subset(g, chr==currchr) #extract only the values that correspond to the current chromosome from the matrix g
	g2 <- g1[,c("blocklength","chrindex", aatype)] #extract the relevant columns to plot with sequential block position (1,2,3,4,...)
	#locpos <- g1[,c("locus")]
	#lp2 <- gsub(paste(currchr,"_",sep=""),"",locpos) #remove underscore separator in locus name
	#g1$locus <- lp2 #assign modified test to data frame
	#g2 <- g1[,c("blocklength","chrindex", aatype)] #extract the relevant columns to plot with genomic location
	g_matrix <- acast(g2, chrindex~blocklength, value.var=aatype) #reshape the data frame into a matrix

	cn=as.numeric(colnames(g_matrix))
	rn=as.numeric(rownames(g_matrix))

	rbot <- list("g_matrix"=g_matrix,"cn"=cn,"rn"=rn)
	return(rbot)
}

heatmap <- function(sub_matrix, infilename, suffix, cn, rn, currchr, aatype)
{
	library(plot3D)
	gm <-t (sub_matrix) #transpose matrix for plotting
	pdf(paste(infilename,suffix,sep=""))
	jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
						 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
	image2D(z=gm, x=cn, y=rn, col = jet.colors(1000), clim = c(0, 1), xlab="Haplotype block length", ylab=paste("Genomic position, chr ",currchr, sep=""), clab=paste(infilename,aatype,"                                           ",sep=" "))
	dev.off()
}


#MAIN
#make heatmaps of retention for optimized or random cores
maxchr=19 #the number of chromosomes
infilenames<-c("PopLocusWeights.txt")

for (infilename in infilenames)
{
	g <- read.table(infilename, header=TRUE, sep="\t")
	#read so that column with locus identifier is treated as a complete string, and no trailing zeroes are dropped
	#g <- read.table(infilename, header=TRUE, sep="\t", colClasses = c(rep("numeric",3),"factor",rep("numeric",8)))
	
	aatypes <- c("NGpr","SSpr","NSpr")
	for (aatype in aatypes)
	{
		for (currchr in 1:maxchr) #make a separate heatmap for each chromosome
		{
			rbot <- shapematrix(g, aatype, currchr)
			heatmap(rbot$g_matrix,infilename,suffix=paste(aatype,"chr",currchr,"pdf",sep="."),rbot$cn,rbot$rn,currchr,aatype)
		}
	}
}


