Purpose of round 2 of Postprocessing is to determine similarity in outcomes between
Rcut, Rcut2, Hclust, and Kmeans methods for selecting subsets of populations that are
"far" apart.

###Example of genotypic coding with 12 SNPs from Arabidopsis###

1			3		X996	1	1	1	1	1	0	0	1	1	0	0	1	0	0	0	0	1	0	1	1	0	0	0	0	0	0	1	1	0	0	0	0	0	1	0	0	0	0	0	1	0	0	1	0	
1			4		X997	1	1	1	0	0	1	0	1	1	0	0	1	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	
2			1		X1153	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
2			2		X1163	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
3			6		X12		1	1	1	0	0	1	1	1	1	0	0	1	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
4			1		X1312	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0

b1
X996	1	1	1	1	1	0	0	1	1	0	0	1
X997	1	1	1	0	0	1	0	1	1	0	0	1
X1153	1	0	0	0	0	1	0	0	0	0	0	0
X1163	0	0	0	0	1	1	0	0	0	0	0	0
X12		1	1	1	0	0	1	1	1	1	0	0	1
X1312	0	0	0	0	1	1	0	0	0	0	0	0

b2
X996		11		11		10		01		10		01
X997		11		10		01		01		10		01
X1153		10		00		01		00		00		00
X1163		00		00		11		00		00		00
X12			11		10		01		11		10		01
X1312		00		00		11		00		00		00

b3
X996			111			110			011			001
X997			111			001			011			001
X1153			100			001			000			000
X1163			000			011			000			000
X12				111			001			111			001
X1312			000			011			000			000

b4
X996				1111			1001			1001		
X997				1110			0101			1001
X1153				1000			0100			0000
X1163				0000			1100			0000
X12					1110			0111			1001
X1312				0000			1100			0000



###Generating a histogram of enrichment index, the potential for enhancing diversity at###
###a single SNP, measured across all associated blocklengths.                          ###

	Calculate values for the histogram, start in folder /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/PostProcessing2:
p=$(pwd)"/";
for d in hclust kmeans Rcut2 Rcut;
  do for j in Pop At Sor;
    do for i in geo pco;
      do printf "%0.2f\n" $(cut -d' ' -f5 "$p$d"/+"$j"SUMEnrichAcrossBlocks.R"$i"Tgen.txt | tail -n +2) | sort -n | sed 's/-0.00/0.00/g' | uniq -c | awk '{print $1"\t"$2}' > "$p$d"/"$j""$i"EnrichmentHistogram.txt;
      done;
    done;
  done;
  
	Plot the histogram:
###Rscript###
wdroot="/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/PostProcessing2"
setwd(wdroot)

t=c("At","Pop","Sor")
e=c("geo","pco")
d=c("hclust","kmeans","Rcut2","Rcut")

pdf(file="EnrichmentHistograms.pdf")
par(mfrow=c(4,3))

for (dx in d)
{
  setwd(file.path(getwd(),dx))
  for (ex in e)
  {
    for (tx in t)
    {
      a=read.table(paste(tx,ex,"EnrichmentHistogram.txt", sep=""))
      plot(a$V2/max(abs(a$V2)),a$V1/max(abs(a$V1)), main=paste(tx,ex,dx,sep=""), type="l", xlab="", ylab="", lwd=0.5, xlim=c(-1,1), ylim=c(0,1))
    }
  }
  setwd(wdroot)
}
dev.off()
######


	Calculate the correlation between enrichment values using geo vs pco data for hclust,kmeans
	Rcut2, and Rcut:
###Rscript###
wdroot="/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/PostProcessing2"
setwd(wdroot)

#t=c("At","Pop","Sor")
d=c("hclust","kmeans","Rcut2","Rcut")

pdf(file="geoXpcoCorrelations.pdf")
par(mfrow=c(4,3))

for (dx in d)
{
  setwd(file.path(getwd(),dx))
  for (tx in t)
  {
    ageo=read.table(paste("+",tx,"SUMEnrichAcrossBlocks.RgeoTgen.txt", sep=""), header=TRUE)
    apco=read.table(paste("+",tx,"SUMEnrichAcrossBlocks.RpcoTgen.txt", sep=""), header=TRUE)

	x=ageo$sum
	y=apco$sum
	rs=cor(x,y,method="spearman")
    plot(x,y,main=paste(tx,dx,sep=""),xlab="Rgeo",ylab="Rpco",bg=rgb(0.35,0.35,0.35,0.06),col=0,pch=21,cex=0.6)
    abline(lm(y~x), col="black", lwd=3)
    text(max(x),min(y),pos=2,paste("rs=",signif(rs,digits=4),sep=""))
  }
  setwd(wdroot)
}
dev.off()
######



	Calculate the correlation between enrichment values using geo vs pco data for Rcut and Rcut2:
###Rscript###
wdroot="/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/PostProcessing1"
setwd(wdroot)

t=c("At","Pop","Sor")
d=c("Rcut","Rcut2")

pdf(file="geoXpcoCorrelations.pdf")
par(mfrow=c(3,2))

for (dx in d)
{
  setwd(file.path(getwd(),dx))
  for (tx in t)
  {
    ageo=read.table(paste("+",tx,"SUMEnrichAcrossBlocks.RgeoTgen.txt", sep=""), header=TRUE)
    apco=read.table(paste("+",tx,"SUMEnrichAcrossBlocks.RpcoTgen.txt", sep=""), header=TRUE)

	x=ageo$sum
	y=apco$sum
	rs=cor(x,y,method="spearman")
    plot(x,y,main=paste(tx,dx,sep=""),xlab="Rgeo",ylab="Rpco",bg=rgb(0.35,0.35,0.35,0.06),col=0,pch=21,cex=0.6)
    abline(lm(y~x), col="black", lwd=3)
    text(max(x),min(y),pos=2,paste("rs=",signif(rs,digits=4),sep=""))
  }
  setwd(wdroot)
}
dev.off()
######


	Above analysis shows a high level of correlation for AtHclust, ~0.92.  Verify this by repeating
	a subset of blocklengths (b150-b200) for AtHclust Rgeo/Rpco.  The correlation should be about the
	same as before, ~0.92.  If not, it suggests that the original analysis may have been messed up.
	The only way I can figure it would get messed up is if Rgeo and Rpco were in fact duplicate
	runs of one or the other.
	
	You need to process the test SUM files up to the point of +AtSUMEnrichAcrossBlocks.*.txt.  Do all
	of this on the cluster then transfer those data to Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/PostProcessing2/hclust/HclustRgeob150-200repeat.
	From that folder do:
###Rscript###
wdroot="/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/PostProcessing2/hclust/HclustRgeob150-200repeat/test6"
setwd(wdroot)

t=c("At")
d=c("hclust")

pdf(file="geoXpcoCorrelationsAtHclustRepeat.pdf")
par(mfrow=c(3,2))

for (dx in d)
{
  for (tx in t)
  {
    ageo=read.table(paste("+",tx,"SUMEnrichAcrossBlocks.RgeoTgen.txt", sep=""), header=TRUE)
    apco=read.table(paste("+",tx,"SUMEnrichAcrossBlocks.RpcoTgen.txt", sep=""), header=TRUE)

	x=ageo$sum
	y=apco$sum
	rs=cor(x,y,method="spearman")
    plot(x,y,main=paste(tx,dx,sep=""),xlab="Rgeo",ylab="Rpco",bg=rgb(0.35,0.35,0.35,0.06),col=0,pch=21,cex=0.6)
    abline(lm(y~x), col="black", lwd=3)
    text(max(x),min(y),pos=2,paste("rs=",signif(rs,digits=4),sep=""))
  }
  setwd(wdroot)
}
dev.off()
######

	If you perform this analysis two times (see folders test1 and test2) you can see that there
	is much correspondence in results (compare files test[1,2]/geoXpcoCorrelationsAtHclustRepeat.pdf).
	The correlation is strong (~0.6), but is not exactly 0.91, possibly because it only samples a portion
	of the full range of blocklengths.  To compare exactly, you need to recalculate with only
	b150-b200 from the original run (folder test3)
	
	When you do this, the correlation is 0.87 and the pattern and range of values look quite 
	different from the repeats performed on glitch (HPC).  It seems that the ceres (another HPC) arabidopsis data is 
	sketchy.  Repeat the analysis from b150-200 on ceres (folder test4).
	
	Doing this returns almost exactly the same distribution as before on ceres, suggesting that
	there is a random number seeding problem.  Upon review of the code I found an error that 
	causes the seed to default to 0 when the program is compiled with no OMP, i.e. in single
	processor mode.  This has been corrected and a new test will be done on glitch and compared
	with the prior distribution from glitch (test 5 vs test1/2).
	
	Test 5 is not really that different from tests1/2, all run on glitch.  Thus the random number
	seeding problem (which -was- actually a real problem, potentially) has no effect on this
	analysis.  If you seed with 0 or you seed with time, you get basically the same result.  That
	provides some confidence that the search algorithm is actually converging on the same solution
	set, regardless of seed.
	
	But there remains a major difference between the results found on glitch and on ceres (see test 6,
	from ceres, with new random number generator seeding mechanism--the same as with seed=0, incidentally).
	This suggests platform-dependent differences in the mechanism of generating random numbers, that
	affect the outcome of the analysis.
	
	
	
	
	
	
	
	

###Notes on quantization methods tested, Rcut, Rcut2, kmeans, hclust
###Rcut and Rcut2 give strange results:

###Figure out whether Rcut or Rcut2 or Random selection produces longer geographic distances.
###Use core size = 2 only, cut all optimized core sets out of RAW files, calculate mean
###walking distance between pairs.


#######BASHPARALLEL#######
myp() {
       x=$1; #get input string, then process into single elements
       p=$(awk '{print $1}' <<<"$x"); #pop1
       q=$(awk '{print $2}' <<<"$x"); #pop2
       
       echo $p"--"$q;
       po=$(grep "^$p " <(echo "$key") | awk '{print $2","$3}');
       qo=$(grep "^$q " <(echo "$key") | awk '{print $2","$3}');
       a=$(Rscript -e 'library("geosphere"); distVincentyEllipsoid(c('$po'),c('$qo'))');
       
       b=$(mktemp ./tmp.XXXXXXXX);
       xo=$(echo -n "$x" | tr " " "," | sed 's/^/(/g' | sed 's/$/)/g');
       echo -n "$xo $po $qo " >> "$b";
       awk '{print $2}' <<<$a >> "$b";
}
export -f myp;


rr="Rcut Rcut2";
ee="geo";
tt="Pop At Sor";

echo "binmethod taxon data meandistance" > DistanceBetweenPops.txt;
for r in $rr;
  do for t in $tt;
    do for e in $ee;
      do echo $r $t $e;
        #define paths
        if [[ $t == "Pop" ]]; then tp="Populus";
        elif [[ $t == "At" ]]; then tp="Arabidopsis";
        elif [[ $t == "Sor" ]]; then tp="Sorghum";
        fi;
        path="/Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/"$tp"/Experiment2/+revision with "$r"/+results/"
        
        #process lat/long file into a key
        key=$(cut -d$'\t' -f6,2,3 "$t"LatLongsp4MOD.txt | tail -n +2 | sort -n -u -t$'\t' -k3,3 | awk '{print $3" "$2" "$1}');
        export key;
        
        #get frequencies of all 2 population optimal sets for blocklength 10, a compromise on wait time and good sampling of complete distribution
        ppath="$path"RAW.b10.R"$e"Tgen.txt;
        freqk=$(grep "([0-9]*,[0-9]*)" "$ppath" | rev | awk '{print $1}' | rev | sort | uniq -c);
        #extract 2 population optimal sets from $freqk
        pops=$(awk '{print $2}' <(echo "$freqk") | sed 's/(//g' | sed 's/)//g' | sed 's/,/ /g');
        #generate a list of paired lat/longs, between which to calculate distance
        
        rm ./tmp.*; #clean up any existing tmp files
        echo "$pops" | parallel --env key myp;
        dists=$(cat ./tmp.* | sort -t' ' -k1,1); 
        rm ./tmp.*;
        
        #multiply distances by frequency of observation
        f=$(awk '{print $1}' <<<"$freqk"); #freqs of observations, in same sort order as dists --verify: diff <(awk '{print $2}' <<<"$freqk") <(awk '{print $1}' <<<"$dists")   
        d=$(awk '{print $4}' <<<"$dists"); #distances between paired populations
        n=$(echo "$f" | awk '{s+=$1}END{print s}'); #total number of pop pairs
        s=$(paste <(echo "$f") <(echo "$d") | awk '{print $1 * $2}' | awk '{s+=$1}END{print s}'); #multiply freq * dist, sum column
        result=$(echo "$s $n 1609" | awk '{print ($1/$2)/$3}'); #calculate average distance between pop pairs, convert meters to miles
        echo "$r $t $e $result" >> DistanceBetweenPops.txt;

        #calculate distance between population pairs chosen at RANDOM
        npops=$(wc -l <<<"$key" | awk '{print $1}'); 
        nreps=$(wc -l <<<"$pops" | awk '{print $1}');
        ranz="";
        for ((zz=0;zz<"$nreps";zz++)); 
        do p1=0; p2=0;
          while [ $p1 -eq $p2 ]; #select two different populations at random
            do p1=$(gshuf -i1-$npops -n1);
               p2=$(gshuf -i1-$npops -n1);
            done;
          ranz+="$p1 $p2"$'\n';
        done;
        ranz=$(sed 's/^$//g' <<<"$ranz"); remove trailing empty line from $ranz
        ranfreqk=$(echo "$ranz" | sort | uniq -c); #remove trailing empty line from $ranz, then count frequency
        ranpops=$(awk '{print $2" "$3}' <<< "$ranfreqk");  
        rm ./tmp.*; #clean up any existing tmp files
        echo "$ranpops" | parallel --env key myp;
        randists=$(cat ./tmp.* | sort -t' ' -k1,1); 
        rm ./tmp.*;
        f=$(awk '{print $1}' <<<"$ranfreqk"); #freqs of observations, in same sort order as dists --verify: diff <(awk '{print $2}' <<<"$freqk") <(awk '{print $1}' <<<"$dists")   
        d=$(awk '{print $4}' <<<"$randists"); #distances between paired populations
        n=$(echo "$f" | awk '{s+=$1}END{print s}'); #total number of pop pairs
        s=$(paste <(echo "$f") <(echo "$d") | awk '{print $1 * $2}' | awk '{s+=$1}END{print s}'); #multiply freq * dist, sum column
        result=$(echo "$s $n 1609" | awk '{print ($1/$2)/$3}'); #calculate average distance between pop pairs, convert meters to miles
        echo "RAN $t $e $result" >> DistanceBetweenPops.txt;
      done;
    done;
  done;



#######ENDBASHPARALLEL#######


		The result is:
	binmethod taxon data meandistance
	Rcut Pop geo 278.903
	Rcut2 Pop geo 240.391
	RAN Pop geo 224.711
	RAN Pop geo 225.634

	Rcut At geo 704.802
	Rcut2 At geo 688.834
	RAN At geo 601.201
	RAN At geo 602.37

	Rcut Sor geo 1929.51
	Rcut2 Sor geo 2178.83
	RAN Sor geo 1884.96
	RAN Sor geo 1898.98

	This shows that Rcut produces pairs of populations that are farther apart for Populus and
	Arabidopsis, but closer together for Sorghum.  For Sorghum Rcut2 works better.  Both methods
	are better than random pairing, producing pairs that are farther apart.
	
	But they are still pretty crappy. The reason is that as you add variables with sqrt(n) categories,
	you soon create a situation where all populations are unique, and M+ enrichment is meaningless
	compared to random.  Probably, all reference variables (2 for geo, 120 for env, 4 for pco) 
	should be converted into a single variable using a multivariate clustering procedure like Kmeans.
	This would ensure that physically proximal populations receive the same clusterID.  With independent
	binning of each variable (e.g. lat AND long each get 7 bins), you can end up with nearby populations
	receiving different lat-long bins.  Test this by creating single summary bins for env, gen, pco
	multivariate data.


	Try some new methods for quantizing the continuous geographical and environmental data.
	Kmeans and hierarchical clustering.  Note that this is just a test, and is fudged a little
	bit in the case of Arabidopsis, where I discovered that pops 15 and 33 have same coordinates.
	This has to be redone, but here, just allow the two populations to not be uniquely specified
	by the bins.  Done by setting break condition to nrows(g)-1 instead of nrows(g).
	
	
###RSCRIPT###
options(error = recover)
rm(list=ls()) 
library(ggplot2)

for (wd in c("sor","pop","at"))
{
	setwd(paste("/Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/PostProcessing/kmeans/",wd,sep=""))

	#c("geo","env","pco")
	for (infile in c("geo"))
	{
		infilename=paste(infile,"DIVAtable.txt",sep="")
		g <- read.table(infilename, header=TRUE, sep="\t")
		x=g[,2:ncol(g)] #get column 2 thru last column

		#calculate Kmeans clusters, from k=2 to k=numind
		kout=g
		for (i in 2:nrow(g))
		{
		  ko=kmeans(x,i)
		  kout[,paste("k",i,sep="")]<-ko$cluster #add cluster coding for k=i to output
		  #test whether each population can be uniquely identified, fudged a bit to accomodate
		  #problem with arabidopsis identical populations
		  if (wd == "at" ) {
		  	nahnah=nrow(g)-1
		  	}
		  	else {
		  	nahnah=nrow(g)
		  	}
		  if ( nrow(unique(kout[(ncol(g)+1):ncol(kout)])) == nahnah ) { break } 
		}
		#write kout to a file
		write.table(kout,file=paste(wd,infile,"kout.txt",sep=""), row.names=FALSE)
		

		#calculate a dendrogram from hierarchical clusters, plot each level
		d=dist(x) #calculate a euclidean distance matrix
		c=hclust(d)
		#plot(c, hang=-1) #plot the dendrogram
		hout=g
		for (i in 2:nrow(g))
		{
		  hout[,paste("h",i,sep="")]<-cutree(c,i)

		  #test whether each population can be uniquely identified
		  if ( nrow(unique(hout[(ncol(g)+1):ncol(hout)])) == nrow(g) ) { break }
		}
		#write hcout to a file
		write.table(hout,file=paste(wd,infile,"hout.txt",sep=""), row.names=FALSE)

	}
}
###ENDRSCRIPT###

	Visualize these clustering schemes for K=2-10:
###BEGINRSCRIPT###
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/PostProcessing/kmeans/sor")
library(ggplot2)

infilename="geoDIVAtable.txt"
g <- read.table(infilename, header=TRUE, sep="\t")
x=g[,2:3] #get Lat/Long columns

#calculate Kmeans clusters for K=2-10
par(mfrow=c(3,3))
for (i in 2:10)
{
  ko=kmeans(x,i)
  y=x
  y[,"cluster"]<-ko$cluster #add cluster designation to lat/long table
  plot(y$Longitude, y$Latitude, col=y$cluster, main=paste("kmeans, K=",i))
}

#calculate a dendrogram from hierarchical clusters, plot each level
par(mfrow=c(1,1))
d=dist(x) #calculate a euclidean distance matrix
c=hclust(d)
plot(c, hang=-1) #plot the dendrogram
hcout=g
par(mfrow=c(3,3))
for (i in 2:10)
{
  z=x
  z[,"cluster"]<-cutree(c,i) #add the groups at each hierarchical level
  plot(z$Longitude, z$Latitude, col=z$cluster, main=paste("hclust, K=",i))
  hcout[,paste("c",i,sep="")]<-cutree(c,i)

  #test whether each population can be uniquely identified
  #cat(nrow(unique(hcout[4:ncol(hcout)])),",", sep="")
  
  if ( nrow(unique(hcout[4:ncol(hcout)])) == nrow(g) ) {
    cat("exiting...")
    break
  }
}

#calculate Rcut breaks
nt=x
nbins=round(sqrt(nrow(x))) #rule of thumb for nbins = square root of n obs
for (i in 1:ncol(x))
{
	f<-cut(x[,i], breaks=nbins, labels=FALSE)
	options(width=10000)
	nt[,i]<-f #replace the original values with the binned values
}
print(nt)

#calculate Rcut2 breaks
nt=x
nbins=round(sqrt(nrow(x))) #rule of thumb for nbins = square root of n obs
for (i in 1:ncol(x))
{
	f<-as.integer(Hmisc::cut2(x[,i],g=nbins))
	options(width=10000)
	nt[,i]<-f #replace the original values with the binned values
}
print(nt)


###ENDRSCRIPT###


	Run M+ to create approximately nloci optimizations on the Kmeans and Hclust bins. First,
	create M+ .dat and .var files:

	Modify the x.b10.RgenTgeo.var file to x.b10.RgeoTgen.[H,K].var. Do this manually.
	Modify the x.b10.envgeopco.dat files, removing the old env/geo/pco bins and 
	adding the newly coded geo bins. Execute below from /Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/PostProcessing/kmeans:
	
for t in at pop sor;
do if [[ $t == "pop" ]]; then tp="PopM";
     elif [[ $t == "at" ]]; then tp="AtSNP4";
     elif [[ $t == "sor" ]]; then tp="SorM";
   fi;
  for c in h k;
    do rev ./"$t"/"$tp".b10.envgeopco.dat | cut -d' ' -f27- | rev > ./"$t"/tmp.txt; #cut off the old bins
      >./"$t"/tmp2.txt;
      while read -r line;
      do m=$(echo "$line" | awk '{print $1}'); #get the first column, the MINDEX, of each incoming line
        binstoadd=$(grep "^$m " ./"$t"/"$t"geo"$c"out.txt | cut -d' ' -f4-); #find the corresponding bins for this MINDEX
        paste -d' ' <(echo "$line") <(echo "$binstoadd") >> ./"$t"/tmp2.txt; #add bins onto the end of the line
      done < ./"$t"/tmp.txt;
      mv ./"$t"/tmp2.txt ./"$t"/"$tp".b10.envgeopco."$c".dat;
      rm ./"$t"/tmp.txt;
    done;
done;

#run M+ with some fair number of replicates (32860) to sample the distribution of best 2pop combos for geo bins:
mpirun -np 24 ./m+ PopM.b10.RgeoTgen.h.var PopM.b10.envgeopco.h.dat -m 2 2 1 32860 ./PopM.b10.h.out.txt;
mpirun -np 24 ./m+ PopM.b10.RgeoTgen.k.var PopM.b10.envgeopco.k.dat -m 2 2 1 32860 ./PopM.b10.k.out.txt;
cd ../at;
mpirun -np 16 ./m+ AtSNP4.b10.RgeoTgen.h.var AtSNP4.b10.envgeopco.h.dat -m 2 2 1 32860 ./AtSNP4.b10.h.out.txt;
mpirun -np 16 ./m+ AtSNP4.b10.RgeoTgen.k.var AtSNP4.b10.envgeopco.k.dat -m 2 2 1 32860 ./AtSNP4.b10.k.out.txt;
cd ../sor;
nice mpirun -np 16 ./m+ SorM.b10.RgeoTgen.h.var SorM.b10.envgeopco.h.dat -m 2 2 1 32860 ./SorM.b10.h.out.txt;
mpirun -np 24 ./m+ SorM.b10.RgeoTgen.k.var SorM.b10.envgeopco.k.dat -m 2 2 1 32860 ./SorM.b10.k.out.txt;
#mpirun -np 24 ./m+ SorM.b10.RgeoTgen.k.var SorM.b10.envgeopco.k.dat -m 3 3 1 32860 ./SorM.b10.k.3.out.txt;


	Calculate the mean distance between points calculated using hclust and kmeans, compare to Rcut,
	Rcut2, and RAN:
	
#######BASHPARALLEL#######
myp() {
       x=$1; #get input string, then process into single elements
       p=$(awk '{print $1}' <<<"$x"); #pop1
       q=$(awk '{print $2}' <<<"$x"); #pop2
       
       echo $p"--"$q;
       po=$(grep "^$p " <(echo "$key") | awk '{print $2","$3}');
       qo=$(grep "^$q " <(echo "$key") | awk '{print $2","$3}');
       a=$(Rscript -e 'library("geosphere"); distVincentyEllipsoid(c('$po'),c('$qo'))');
       
       b=$(mktemp ./tmp.XXXXXXXX);
       xo=$(echo -n "$x" | tr " " "," | sed 's/^/(/g' | sed 's/$/)/g');
       echo -n "$xo $po $qo " >> "$b";
       awk '{print $2}' <<<$a >> "$b";
}
export -f myp;


rr="h k";
ee="geo";
tt="Pop At Sor";

echo "binmethod taxon data meandistance" > DistanceBetweenPops.txt;
for r in $rr;
  do for t in $tt;
    do for e in $ee;
      do echo $r $t $e;
        #define paths
        if [[ $t == "Pop" ]]; then tp="PopM";
        elif [[ $t == "At" ]]; then tp="AtSNP4";
        elif [[ $t == "Sor" ]]; then tp="SorM";
        fi;
        path="/Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/PostProcessing/kmeans/"
        
        #process lat/long file into a key
        key=$(cut -d$'\t' -f6,2,3 "$t"LatLongsp4MOD.txt | tail -n +2 | sort -n -u -t$'\t' -k3,3 | awk '{print $3" "$2" "$1}');
        export key;
        
        #get frequencies of all 2 population optimal sets for blocklength 10, a compromise on wait time and good sampling of complete distribution
        ppath="$path""$t"/"$tp".b10."$r".out.txt;
        freqk=$(grep "([0-9]*,[0-9]*)" "$ppath" | rev | awk '{print $1}' | rev | sort | uniq -c);
        #extract 2 population optimal sets from $freqk
        pops=$(awk '{print $2}' <(echo "$freqk") | sed 's/(//g' | sed 's/)//g' | sed 's/,/ /g');
        #generate a list of paired lat/longs, between which to calculate distance
        
        rm ./tmp.*; #clean up any existing tmp files
        echo "$pops" | parallel --env key myp;
        dists=$(cat ./tmp.* | sort -t' ' -k1,1); 
        rm ./tmp.*;
        
        #multiply distances by frequency of observation
        f=$(awk '{print $1}' <<<"$freqk"); #freqs of observations, in same sort order as dists --verify: diff <(awk '{print $2}' <<<"$freqk") <(awk '{print $1}' <<<"$dists")   
        d=$(awk '{print $4}' <<<"$dists"); #distances between paired populations
        n=$(echo "$f" | awk '{s+=$1}END{print s}'); #total number of pop pairs
        s=$(paste <(echo "$f") <(echo "$d") | awk '{print $1 * $2}' | awk '{s+=$1}END{print s}'); #multiply freq * dist, sum column
        result=$(echo "$s $n 1609" | awk '{print ($1/$2)/$3}'); #calculate average distance between pop pairs, convert meters to miles
        echo "$r $t $e $result" >> DistanceBetweenPops.txt;

        #calculate distance between population pairs chosen at RANDOM
        npops=$(wc -l <<<"$key" | awk '{print $1}'); 
        nreps=$(wc -l <<<"$pops" | awk '{print $1}');
        ranz="";
        for ((zz=0;zz<"$nreps";zz++)); 
        do p1=0; p2=0;
          while [ $p1 -eq $p2 ]; #select two different populations at random
            do p1=$(gshuf -i1-$npops -n1);
               p2=$(gshuf -i1-$npops -n1);
            done;
          ranz+="$p1 $p2"$'\n';
        done;
        ranz=$(sed 's/^$//g' <<<"$ranz"); remove trailing empty line from $ranz
        ranfreqk=$(echo "$ranz" | sort | uniq -c); #remove trailing empty line from $ranz, then count frequency
        ranpops=$(awk '{print $2" "$3}' <<< "$ranfreqk");  
        rm ./tmp.*; #clean up any existing tmp files
        echo "$ranpops" | parallel --env key myp;
        randists=$(cat ./tmp.* | sort -t' ' -k1,1); 
        rm ./tmp.*;
        f=$(awk '{print $1}' <<<"$ranfreqk"); #freqs of observations, in same sort order as dists --verify: diff <(awk '{print $2}' <<<"$freqk") <(awk '{print $1}' <<<"$dists")   
        d=$(awk '{print $4}' <<<"$randists"); #distances between paired populations
        n=$(echo "$f" | awk '{s+=$1}END{print s}'); #total number of pop pairs
        s=$(paste <(echo "$f") <(echo "$d") | awk '{print $1 * $2}' | awk '{s+=$1}END{print s}'); #multiply freq * dist, sum column
        result=$(echo "$s $n 1609" | awk '{print ($1/$2)/$3}'); #calculate average distance between pop pairs, convert meters to miles
        echo "RAN $t $e $result" >> DistanceBetweenPops.txt;
      done;
    done;
  done;



#######ENDBASHPARALLEL#######
	
H and K methods are clearly better than Rcut and Rcut2.  The result for Sorghum is exaggerated a bit due to India outlier,
but Pop and At are solid.


binmethod taxon data meandistance
h Pop geo 393.813
k Pop geo 372.061
Rcut Pop geo 278.903
Rcut2 Pop geo 240.391
RAN Pop geo 224.711
RAN Pop geo 225.634
RAN Pop geo 229.537
RAN Pop geo 214.488

h At geo 836.669
k At geo 846.194
Rcut At geo 704.802
Rcut2 At geo 688.834
RAN At geo 601.201
RAN At geo 602.37
RAN At geo 633.896
RAN At geo 573.775

h Sor geo 3583.93
k Sor geo 2781.42
Rcut Sor geo 1929.51
Rcut2 Sor geo 2178.83
RAN Sor geo 1884.96
RAN Sor geo 1898.98
RAN Sor geo 1887.07
RAN Sor geo 1760.39







###Notes on Rscript for Kmeans and Hclust quantization to be integrated into applescript "5b. Calc generalized variance for env data.scpt"

		Kmeans
#!/usr/bin/Rscript
g <- read.table(file=latlongsourcefile, header=TRUE, sep="\t", row.names="MINDEX")
x=g
kout=g
for (i in 2:nrow(g))
{
ko=kmeans(x,i) #calculate kmeans clusters where k=i
kout[,paste('k',i,sep='')]<-ko$cluster #add cluster coding for k=i to output
#test whether each population can be uniquely identified
if ( nrow(unique(kout[(ncol(g)+1):ncol(kout)])) == nrow(g) ) { break } 
}
#print bins, omitting original input data
print(kout[,(ncol(g)+1):ncol(kout)])
		






		Hierarchical clusters
#!/usr/bin/Rscript
g <- read.table(file=latlongsourcefile, header=TRUE, sep="\t", row.names="MINDEX")
x=g
#calculate a dendrogram from hierarchical clusters, plot each level
d=dist(x) #calculate a euclidean distance matrix
c=hclust(d)
hout=g
#add the cluster coding for each depth in the dendrogram
for (i in 2:nrow(g))
{
hout[,paste('h',i,sep='')]<-cutree(c,i)
#test whether each population can be uniquely identified
if ( nrow(unique(hout[(ncol(g)+1):ncol(hout)])) == nrow(g) ) { break }
}
#write hcout to a file
print(hout[,(ncol(g)+1):ncol(hout)])
#write.table(hout,file=paste(wd,infile,"hout.txt",sep=""), row.names=FALSE)


###ENDRSCRIPT###

