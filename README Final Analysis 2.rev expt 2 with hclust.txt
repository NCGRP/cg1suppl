	Starting 2/28/2017.
	Experiment 2 is to visualize the genomic distribution of diversity enhancement in core
	via geographical or environmental variables.  Goal is to produce "heatmaps" showing
	the enhancement across each chromosome for all blocklengths.
	
	This revision uses R's hclust to form geographical and environmental bins. Rcut and Rcut2
	gave strange results. See README Postprocessing expt2.txt for some notes on this.  Note 
	also that the Arabidopsis data has been modified to 40 populations from the 41 before
	(two had identical coordinates, which I missed previously).  Also that the outlying
	Sorghum population in India has been removed, reducing the number of populations from 23
	to 22.
	

***POPULUS***
	Make a new set of haplotypista output files that contain start and end points for each block
	
ssh compute-0-1
./haplotypista -i PopHaplotypistaIn.txt -o PopSNP -l PopSNPlog1.txt -b 1 1 -m ? -p 2 &
./haplotypista -i PopHaplotypistaIn.txt -o PopSNP -l PopSNPlog2.txt -b 2 2 -m ? -p 2 &
./haplotypista -i PopHaplotypistaIn.txt -o PopSNP -l PopSNPlog3.txt -b 3 6 -m ? -p 2 &
./haplotypista -i PopHaplotypistaIn.txt -o PopSNP -l PopSNPlog4.txt -b 7 14 -m ? -p 2 &
./haplotypista -i PopHaplotypistaIn.txt -o PopSNP -l PopSNPlog5.txt -b 15 30 -m ? -p 2 &
./haplotypista -i PopHaplotypistaIn.txt -o PopSNP -l PopSNPlog6.txt -b 31 60 -m ? -p 2 &
./haplotypista -i PopHaplotypistaIn.txt -o PopSNP -l PopSNPlog7.txt -b 61 200 -m ? -p 2 &
	~40 minutes

	Undouble the haplotypista output using the following bash script:
**********BASH**********
 MyUndub() {
 b=$1;
 #file root (f) and outfile root (o)
 f=PopSNP;
 o=PopM;
 
 # make doubled version of header to reflect diploid genotypes
    echo $f.b$b;
    z="";
    for i in {1..6}; #6 header lines
      do a=$(head -"$i" $f.b$b | tail -n -1 | tr ' ' '\n'); #extract ith line
        a=$(paste -d' ' <(echo "$a") <(echo "$a") | tr '\n' ' '); #double it side by side, then make it one line
        z=$z"$a
";
     done;
    
   # un-double lines containing genotypes, for each individual
    t4="";
    for i in {1..251}; #half of (total line number minus 6)
      do start=$(expr $i \* 2 + 5);
        t1=$(tail -n +$start $f.b$b | head -1 | tr ' ' '\n');
        t2=$(tail -n +$(expr $start + 1) $f.b$b | head -1 | tr ' ' '\n');
        t3=$(paste -d' ' <(echo "$t1") <(echo "$t2") | tr '\n' ' ');
        t4=$t4"$t3
";
      done;

    # assemble the un-doubled matrix with the doubled header lines
    { echo -n "$z"; ghead -c -1 <<<"$t4" | cut -d' ' -f2-; } > $o.b$b;
 }
export -f MyUndub;
 
seq 1 200 | parallel --env MyUndub MyUndub;
**********BASHEND**********
	
	extract the NG:SS:NS values for each locus from the headers of the PopM.b* files.  Make a table
	called "PopLocusWeights". These values will serve to compare with enrichment values from m+1 of the same loci.

**********BASH**********
f="PopLocusWeights";
header=$'genomeindex\tblocklength\tblocklengthindex\tlocus\tchr\tchrindex\tNGcount\tSScount\tNScount\tNGpr\tSSpr\tNSpr\tblockstart\tblockend';
> $f.txt;
for ((b=1;b<=200;b++));
do echo PopM.b$b;
  chr=$(head -1 PopM.b$b | tr " " "\n" | awk 'NR%2'); #extract line 1 (chr), get odd numbered lines
  loc=$(head -3 PopM.b$b | tail -1 | tr " " "\n" | awk 'NR%2'); #extract line 3 (midpt), get odd numbered lines
  st=$(head -5 PopM.b$b | tail -1 | tr " " "\n" | awk 'NR%2'); #extract blockstart
  en=$(head -6 PopM.b$b | tail -1 | tr " " "\n" | awk 'NR%2'); #extract blockend
  w=$(wc -w <<< $loc); #get the number of loci
  blI=$(jot $w); #compute an index value specific to this blocklength
  chrloc=$(paste -d'_' <(echo "$chr") <(echo "$loc"));#fuse chr and locus midpoint into locus name
  
  #calculate an index by chromosome, within blocklength
chrI=$(i=1;
       prevc=0;
       for c in $chr;
         do if [ $c -eq $prevc ];
            then i=$(($i + 1));
                echo $i;
            else i=1;
                echo $i;
            fi;
         prevc=$c;
         done;
       )
               
  aacats=$(head -4 PopM.b$b | tail -1 | tr " " "\n" | awk 'NR%2' | sed "s/:/"$'\t'"/g"); #extract line 4 (NG:SS:NS), get odd numbered lines
  aapr=$(echo "$aacats" | awk '{sum = $1 + $2 + $3; ng = ($1/sum); ss = ($2/sum); ns = ($3/sum); print ng "\t" ss "\t" ns;}'); #sum NGSSNS counts, and find the proportion of each category

  paste -d$'\t' <(echo "$blI") <(echo "$chrloc") <(echo "$chr") <(echo "$chrI") <(echo "$aacats") <(echo "$aapr") <(echo "$st") <(echo "$en") | sed "s:^:$b"$'\t'":g" >> $f.txt;
done;
awk 'BEGIN { OFS = "\t" } ; {print NR,$0}' $f.txt > tmp.txt; #add a sequential numeric index value (line number) to each row
#nl $f.txt | sed 's/ //g' > tmp.txt; #add a sequential numeric index value (line number) to each row
#cat tmp.txt | sed 's/\./_/g' > tmp.txt; #change . to _ in locus description so r doesn't think its a number
{ echo "$header" ; cat tmp.txt; } > $f.txt;
rm tmp.txt;
**********BASHEND**********

	create a file '+PopSNPtoBlockMap.txt' that shows which SNPs fall in which blocks for
	b2-b200.
**********BASHPARALLEL**********
myp() {
  b=$1;
  #define parameters
  nchr=19; #number of chromosomes
  snps=$(cat "$path""snps.tmp.txt");
  
    echo "blocklength=$b";
    bnames=b"$b"; #bnames contains the haplotype block names associated with each snp for a given blocklength
    bla=$(grep "^[0-9]*"$'\t'"$b"$'\t' "$path""$p"LocusWeights.txt); #get all lines for blocklength b

    #repeat through chromosomes
    for i in $(seq 1 $nchr);
      do echo "  chr""$i";
        snpsizes=$(echo "$snps" | grep ^"$i"_ | sed 's/'$i'_//g'); #get snp sizes for this chromosome and blocklength 1
        
        declare -a cname=( $(echo "$bla" | grep $'\t'"$i_" | awk '{print $4}') ); #locus names for this chr & blocklength
        declare -a cstart=( $(echo "$bla" | grep $'\t'"$i_" | awk '{print $13}') ); #start points for each locus in $cname
        declare -a cend=( $(echo "$bla" | grep $'\t'"$i_" | awk '{print $14}') ); #end points for each locus in $cname
        clen=$(echo ${#cname[@]}); #total length of the arrays for chromosome i
       
        startpt=0;
        for s in $snpsizes;
          do found=0;
            #test whether snp position is found within haplotype block by seeing whether it is both >= start point and <= end point
            for ((c=$startpt;c<$clen;c++)); 
              do if (( $(bc <<< "$s >= ${cstart[$c]}") )) && (( $(bc <<< "$s <= ${cend[$c]}") )); then
                  bnames="$bnames"$'\n'${cname[$c]}; #add the haplotype block name to the list that maps them to SNP name
                  found=1;
                  startpt=$c; #make search more efficient by ignoring earlier blocks in future searches
                fi;
                if [[ $found == 1 ]]; then break; fi;
              done;
            #deal with situation where SNP position is not found in a haplotype block (this will
            #sometimes happen at the end of chromosomes)
            if [[ $found == 0 ]]; then
              bnames="$bnames"$'\n'"--"; 
            fi;
          done; #s
      done; #i
      
      echo "$bnames" > "$path$p""b$b.tmp";
}
export -f myp;

#define parameters
maxbl=200; #maximum blocklength
p="Pop"; export p;
path=$(pwd)"/"; export path; #path to current folder
a=$(grep "^[0-9]*"$'\t'"1"$'\t' "$path""$p"LocusWeights.txt); #get all lines with blocklength 1
snps=$(echo "$a" | awk '{print $4}'); #get a list of all SNP positions
echo "$snps" > "$path""snps.tmp.txt"; #put in a file accessible by all procs so it doesn't have to be passed by parallel or loaded into the environment
outf=$(head -1 "$path""$p"LocusWeights.txt | cut -d$'\t' -f1-6);
outf="$outf"$'\n'$(echo "$a" | cut -d$'\t' -f1-6);

#begin parallel
seq 2 $maxbl | parallel --sshloginfile ~/machines --jobs 24 --env myp --env path --env p myp;

#reassemble
echo "reassembling...";
comm=""; #create a list of tmp files to paste together
for b in $(seq 2 $maxbl);
  do comm="$comm""$p""b$b.tmp ";
  done;
paste -d$'\t' <(echo "$outf") $comm > "$path""+$p""SNPtoBlockMap.txt"; #consecutively paste results from each blocklength

#clean up
rm snps.tmp.txt;
for b in $(seq 2 $maxbl);
  do rm "$p""b$b.tmp"; #clean up
  done;
**********BASHEND**********
	Takes ~45 min.


	make some test plots of the frequency of NG, NS, and SS categories in each haplotype block using
	this data (PopLocusWeights.txt) using R.  Here you use an r script called PopLocusWeightsHeatmap.r
	plots are as expected, warmer values for NG, and SS and NS are about the same and occur at relatively
	low frequency in the haplotype blocks.  You can see the effect of anomalous frequencies propagated thru
	the blocklengths in, e.g., "PopLocusWeights.txtNGpr.chr.7.pdf", as different colored stripes.
	This stuff is here: /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Populus/Experiment2/+results/+test\ heatmaps

###REVISED REGION###
	remake the M+ .dat and .var input files with the new geo, env, pco binned variables, calculated with R's 'hclust'
	and 'kmeans':

	Convert haplotypista data sets to M+:
	Run script "2. Convert haplotypista to M+.scpt"  #set so block midpoint positions end up in var file. 
	Use 1 script instance on MacMini (Mavericks).  The script takes
	as input the PopM files, which have diploid data in a single line per individual.  These input files
	are created by the "un-doubling" process.
	Takes ~1 hr.
	
	Add Geo, Env, PCOEnv data to end of M+ files for various optimizations:
	To get new bins using R's 'hclust' function, run script "5b. Calc generalized variance for env data.scpt" with the settings 
	CenterAndScale="bin", cuttype="hclust", and binonly="yes.  Rename the output files like "Athclustenvinput7.txt", "Popkmeanspcoinput7.txt",
	"Athclustgeoinput7.txt". Then run script "7. Add target data to M+.scpt",
	which will add the new binned values to the .dat and .var files from above.  This was performed 
	with one processor on MacMini (Mavericks).
	Takes ~30 minutes.
	
	run script "8d. Manage M+ runs on clusterPopMOD.sh" on glitch (HPC).  this parallelizes on locus, and post processes the m+1 output,
	leaving a SUM file with stats, and a RAW file which is an archive of all the m+1 output, concatenated:
	
screen;
cd /home/reevesp/mplusrunsPop/blockwise/usingHclust;
./8d*;
cd /home/reevesp/mplusrunsPop/blockwise/usingKmeans;
./8d*;


###

	calculate mean M+ enrichment across core sizes, for each locus, for each blocklength. This uses the
	SUMfiles as input and produces the +STATS files.  The parallel script takes ~30 minutes.


**********BASHPARALLEL-FASTER-RUNINSCREEN**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all core sizes, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum "\t" sum/NF "\t" sqrt(var) "\t" NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

main() {
        l=$1;
        ofile="$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt"; #temporary outfile for each parallel process
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";

        boo=$(grep -F $l $sumfileIN | awk '{ print $13 - $10 }' | tr "\n" " "); #extract all core sizes for this locus, calculate target enrichment
        st=$(mysd "$boo"); #calculate sum, mean, sd, and n using function mysd()
        echo "$b"$'\t'"$l"$'\t'"$st" > $ofile; #write stats to output file
}
export -f main;

#parameters
t="Populus";
elist="geo pco";
bmin=1; #min blocklength
bmax=200;
nnode=6; #number of cluster nodes
blist=$(seq $bmin $bmax); #blocklengthrange
p=$(pwd)"/"; export p;

#parallelize on locus, within loop on blocklength
for e in $elist;
  do echo $e;
    export e;
    for b in $blist;
      do echo "b=$b";
        export b;
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";
        uloc=$(cut -d$'\t' -f2 "$sumfileIN" | tail -n +2 | uniq); #find unique loci from all loci in column 2
        numloc=$(echo "$uloc" | wc -l);
        
        #parallel step
        #sends $gnuN records to each node via parallel step #1, this is piped to parallel step #2, which starts 96 jobs per node from the instruction set it has received
        echo "processing $sumfileIN...";
        gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node
        #echo "$uloc" | parallel --sshloginfile ~/machines --jobs 24 --env main --env mysd --env e --env b --env p main;
        echo "$uloc" | parallel --sshloginfile ~/machines --jobs 1 --env main --env mysd --env e --env b --env p --pipe -N"$gnuN" /home/reevesp/bin/parallel -j96 --env main --env mysd --env e --env b --env p main;

        #consolidate results for each locus into one file
        ofile="+STATS.R"$e"Tgen.b"$b".txt"; #temporary outfile
        > "$p"$ofile;
        for l in $uloc;
          do echo "concatenating loc file $l to +STATS.R"$e"Tgen.b"$b".txt ...";
            cat "$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt" >> "$p"$ofile;
          done;
          #test whether all were concatenated
          #s1=$(wc -l "$p""+STATS.R"$e"Tgen.loc"*".b"$b".txt" | grep total | awk '{print $1}');
          s1=$(ls "$p" | grep -c "loc"); #count the number of files with "loc" in their name using grep. can't wc -l because too many arguments
          s2=$(wc -l "$p"$ofile |  awk '{print $1}');
         if [ $s1 == $s2 ]; 
          then find "$p" -maxdepth 1 -name "*loc*" -delete; #remove files when there are too many for rm
          #rm "$p""+STATS.R"$e"Tgen.loc"*".b"$b".txt"; #remove files if number of lines in concat file is the same as the sum of all input files
          else echo "s1="$s1", s2="$s2". Consolidated file "$ofile" no good. Quitting..." > err.txt;
            exit 1;
        fi;
      done; #$b
      
    #concatenate files for this Rgeo/Rpco
    statsfile="$p""+STATS."$t".R"$e"Tgen.txt";
    echo blocklength$'\t'locus$'\t'sum$'\t'mean$'\t'sd$'\t'n > $statsfile;
    >"$statsfile"TMP;
    for b in $blist;
      do echo "concatenating +STATS.R"$e"Tgen.b"$b".txt ...";
        cat "$p""+STATS.R"$e"Tgen.b"$b".txt" >> "$statsfile"TMP;
      done;
    cat "$statsfile"TMP >> "$statsfile";

    #test whether all individual +STATS files have been added to the final concatenated file using number of lines
    s1=$(wc -l +STATS.R"$e"Tgen.b*.txt | grep total | awk '{print $1}');
    s2=$(wc -l "$statsfile" | awk '{print $1}');
    s2=$(( $s2 - 1 ));
    if [ $s1 == $s2 ]; 
      then rm +STATS.R"$e"Tgen.b*.txt; #remove files if number of lines in concat file is the same as the sum of all input files
        rm "$statsfile"TMP;
      else echo "the number of lines ain't the same. something is 'crewed.";
    fi;

  done; #$e
**********BASHEND**********
	Takes about 30 minutes.

	add the proportion of NG,SS,NS sites from PopLocusWeights.txt to the +STATS summary file,
	generating the +v2STATS file.  Requires the PopLocusWeights.txt file as input.  

**********BASH**********
weightsfile="PopLocusWeights.txt";
elist="geo pco";
for e in $elist;
  do statsfile="+STATS.Populus.R"$e"Tgen.txt";
    s=$(cut -d$'\t' -f2 $statsfile | tail -n +2 | md5sum); #get md5sum of locus ID column from the +STATS file
    w1=$(cut -d$'\t' -f4 $weightsfile | sed 's/_/./g'); #get the locus ID less blocklength index column, replace _ with . for md5 check
    w2=$(cut -d$'\t' -f3 $weightsfile); #get the blocklength index column
    w=$(paste -d'.' <(echo "$w1") <(echo "$w2") | tail -n +2 | md5sum); #get md5sum of locus ID column from the Weights file
    if [ "$s" == "$w" ]; #verify that files are in same order using a checksum of the unique locus ID column
      then boo=$(sed 's/blocklength/blocklength2/g' $weightsfile | sed 's/locus/locus2/g'); #modify the header in the Weights file so that there are no redundant column names after combining
        paste -d$'\t' $statsfile <(echo "$boo") > "+v2STATS.Populus.R"$e"Tgen.txt"; #paste the statsfile and the weights file together
      else echo "md5sums do not match, $e, aborting.";
        kill -INT $$; #terminate the script, return to the shell
    fi;
  done;
**********BASHEND**********
	Takes a couple seconds

	Plot summed enrichment by genomic position across blocklengths as a heatmap. Use Rscript
	PopGenomicGeography.r. Uses the +v2STATS files as input.


	Using the file +*SNPtoBlockMap.txt, calculate the sum of M+ enrichment value across blocklengths,
	for each SNP, i.e. for each position in the genome.  These data will be extracted from the 
	+v2STATS* files in several steps.  The first makes a file containing a table with the max 
	enrichment value at each site for each block, *EnrichAcrossBlocks*. Step 1:

**********BASHPARALLEL**********
myp() {
  e=$1;
  b=$2;
  v2name="$p""+v2STATS."$v2p".R"$e"Tgen.txt";
        
        echo "e=$e  b=$b";
        bout="$p""rr.R"$e"Tgen.b""$b"".tmp"; #temporary outfile

        #cut out the locus id and sum columns for current blocklenth
        gg=$(grep ^$b$'\t' "$v2name" | cut -d$'\t' -f1-3);
        
        #cut column with locus names for current blocklength from SNPtoBlockMap
        if [ $b = 1 ]; then
          cc=$(cut -d$'\t' -f4 "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        else col=$(( $b + 5 )); #calculate the column to cut out for each blocklength
          cc=$(cut -d$'\t' -f$col "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        fi;

        #search for locus names in $gg, retrieve enrichment value
        rr="b$b";

#put another parallel here if necessary
        for l in $cc;
          do if [[ $l = "--." ]]; then
                r="--"; #test for empty enrichment value
              else
                r=$(grep -m 1 ^"$b"$'\t'"$l" <<<"$gg" | awk '{print $3}');
              fi;
            rr+=$'\n'$r; #add current locus enrichment value to list
          done;

         #write out the result
         echo "$rr" > "$bout";
         #truncate -s -1 "$bout"; #remove the hanging newline, use 'gtruncate' only on osx
}
export -f myp;

#parameters
v2p="Populus"; export v2p;
t="Pop";
minbl=1;
maxbl=50; #50,200
blist=$(seq $minbl $maxbl);
elist="geo pco";

p=$(pwd)"/"; export p;
tname="$p""+""$t""SNPtoBlockMap.txt"; export tname;

#begin parallel
parallel --sshloginfile ~/machines --jobs 24 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist;
  
#reconstitute
echo "reconstituting..."
for e in $elist;
  do echo "e=$e";
    comm=""; #create a list of tmp files to paste together
    for b in $blist;
      do comm+="$p""rr.R"$e"Tgen.b$b.tmp ";
    done;
    pos=$(cut -d$'\t' -f4 $tname | sed 's/'^.*_'//g' | sed 's/locus/pos/g');
    h=$(cut -d$'\t' -f1-6 $tname | tr $'\t' ' ');
    h2=$(paste -d' ' <(echo "$h") <(echo "$pos"));
    paste -d' ' <(echo "$h2") $comm > "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt";
    rm $comm; #clean up
  done;#e
  
**********BASHEND**********
	Takes ~20 minutes for maxbl=200.

	Process each line (SNP position) of *EnrichAcrossBlocks*.txt to get sum across all blocklengths
	of cumulative enrichment across all core sizes of max enrichment value from 10 reps of M+. Step 2:

**********BASH**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all blocklengths, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       #returns space delimited string like "sum mean sd n"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum " " sum/NF " " sqrt(var) " " NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

myp() {
    l=$1

    p=$(echo "$l" | cut -d' ' -f8-); #cut out the portion of the line with max enrichment values

    #skip loci that have missing blocks
    if [[ $p != *"--"* ]]; then
      m=$(mysd "$p"); #calculate sum, mean, sd, n
      chr=$(echo "$l" | cut -d' ' -f5);
      pos=$(echo "$l" | cut -d' ' -f7);
      i=$(echo "$l" | cut -d' ' -f1); #get genome index
      #loc=$chr"."$pos"."$(echo "$l" | cut -d' ' -f6);
      loc=$chr"."$pos"."$i;
      echo "snp=$i";

      echo $i" "$chr" "$pos" "$loc" "$m > "$path"$i".tmp"; #add new values to growing string, write to tmp file
    fi;
}
export -f myp;

#parameters
elist="geo pco";
t="Pop";
path=$(pwd)"/"; export path;

#parallel
for e in $elist;
do echo $e;
  nnode=8;
  numloc=$(wc -l +PopEnrichAcrossBlocks.RgeoTgen.txt | awk '{print $1}');
  numloc=$(( $numloc - 1 ));
  gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node

  #tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --sshloginfile ~/machines --jobs 24 --env myp --env mysd --env path myp;
  tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --sshloginfile ~/machines --jobs 1 --env myp --env mysd --env path --pipe --round-robin /home/reevesp/bin/parallel -j96  --env myp --env mysd --env path myp;
  
  #assemble results
  echo "concatenating...";
  ofile="+"$t"SUMEnrichAcrossBlocks.R"$e"Tgen.txt";
  o="genomeindex chr pos loc sum mean sd n";
  echo "$o" > $ofile;
  
  find "$path" -name "*.tmp" -print0 | xargs -0 cat | sort -n -t' ' -k1 >> $ofile;
  
  #clean up
  echo "clean up...";
  rm "$path"*.tmp;
done;
**********ENDBASH**********
	Takes ~3 minutes.

	Use R to make some pretty plots of significant enrichment across the genome. Switch setwd between
	hclust and kmeans.

**********BEGINR**********
#install.packages("fitdistrplus")
#install.packages("plotrix")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Populus/Experiment2/+revision\ with\ hclust/+results/hclust/+PopHclustSummary\ files/ToMaxBL50")
library(fitdistrplus)
library(plotrix)

t<-"Pop"
usezero="no" #this sets whether to include significant sites only when they are also > or < 0
             #to use the zero criterion, say "yes", to just use the ecdf around the mean, say "no"
pvals<-c(0.001,0.05)
devvals<-c("pdf","svg","png")
for (pp in pvals)
{
	for (zg in devvals)
	{
		zgg<-get(zg) #e.g convert string "pdf" to function pdf now named 'zgg'
		if ( zg == "pdf" | zg == "svg" ) { zgg(file=paste(t,"M+GenomeScan",pp, ".", zg, sep="")) }
		else { 
		par(mar=c(1,1,1,1))
		zgg(file=paste(t,"M+GenomeScan",pp, ".", zg, sep=""),width=2000,height=2000,res=300) } #call zgg function with relevant parameters
		
		#zg(file=paste(t,"M+GenomeScan",pp, zg, sep=""))
		#pdf(file=paste(t,"M+GenomeScan",pp,".pdf", sep=""))
		#png(file=paste(t,"M+GenomeScan",pp,".png", sep=""))
		par(mfrow=c(2,1))

			i="geo"
			g <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",i,"Tgen.txt", sep=""), header=TRUE, sep=" ")
			mycdf <- ecdf(g$sum) #calculate the cdf of the empirical distribution
			p <- mycdf(g$sum) #get the probability of each observation (or lower)
			g$pecdf <- p #add new column to g, containing the probability from the empirical cdf
			g$index <- 1:nrow(g) #add a column with the sequential index
			h <- g[which(g$pecdf > 1-pp),] #get snps that are highly enriched
			if ( usezero == "yes" ) {
				hgeo <- h[which(h$sum>0),] #remove snps where sum<0, i.e. less than random
			}
			else {
				hgeo <- h # do not remove snps based on positive or negative
			}
			l <- g[which(g$pecdf < pp),] #get snps that are lowly enriched
			if ( usezero == "yes" ) {
				lgeo <- l[which(l$sum<0),] #remove snps where sum>0, i.e. more than random
			}
			else {
				lgeo <- l # do not remove snps based on positive or negative
			}
			highposgeo<-hgeo$index #get the sequential position of the highly enriched
			lowposgeo<-lgeo$index #get the sequential position of the deficient

			plot(g$index, g$sum, type="l", main=paste(t,i,pp,sep=""))
			color.scale.lines(g$index, g$sum, col=factor(g$chr))
			points(g$index[highposgeo], g$sum[highposgeo], col = "red", cex=0.7) #mark significant values
			points(g$index[lowposgeo], g$sum[lowposgeo], col = "blue", cex=0.7) #mark significant values

			i="pco"
			g <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",i,"Tgen.txt", sep=""), header=TRUE, sep=" ")
			mycdf <- ecdf(g$sum) #calculate the cdf of the empirical distribution
			p <- mycdf(g$sum) #get the probability of each observation (or lower)
			g$pecdf <- p #add new column to g, containing the probability from the empirical cdf
			g$index <- 1:nrow(g) #add a column with the sequential index
			h <- g[which(g$pecdf > 1-pp),] #get snps that are highly enriched
			if ( usezero == "yes" ) {
			hpco <- h[which(h$sum>0),] #remove snps where sum<0, i.e. less than random
			}
			else {
				hpco <- h
			}
			l <- g[which(g$pecdf < pp),] #get snps that are lowly enriched
			if ( usezero == "yes" ) {
			lpco <- l[which(l$sum<0),] #remove snps where sum>0, i.e. more than random
			}
			else {
				lpco <- l
			}
			highpospco<-hpco$index #get the sequential position of the highly enriched
			lowpospco<-lpco$index #get the sequential position of the deficient

			plot(g$index, g$sum, type="l", main=paste(t,i,pp,sep=""))
			color.scale.lines(g$index, g$sum, col=factor(g$chr))
			points(g$index[highpospco], g$sum[highpospco], col = "red", cex=0.7) #mark significant values
			points(g$index[lowpospco], g$sum[lowpospco], col = "blue", cex=0.7) #mark significant values

			#calculate intersection of well-collected significant snps
			isect <- hgeo[is.element(hgeo$pos, intersect(hgeo$pos,hpco$pos)),]
			write.table(isect, file=paste(t,"HighCommonToGeoPco",pp,".txt", sep=""), row.names=FALSE, quote=FALSE)

			#calculate intersection of poorly-collected significant snps
			isect <- lgeo[is.element(lgeo$pos, intersect(lgeo$pos,lpco$pos)),]
			write.table(isect, file=paste(t,"LowCommonToGeoPco",pp,".txt", sep=""), row.names=FALSE, quote=FALSE)
	
			#create BED formatted output for hgeo, hpco. verify first whether there exist any snps that are
			#both significantly enriched and better than random (or significantly suppressed and worse than random.)
			#if not, do not print an output table.
			if ( length(hgeo$pos) != 0 ) {
			hgeobed <- data.frame(chr=paste("Chr", hgeo$chr, sep=""), chromStart=hgeo$pos, chromEnd=hgeo$pos+1)	
			write.table(hgeobed, file=paste(t,"hGeo",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
			}
			if ( length(hpco$pos) != 0 ) {
			hpcobed <- data.frame(chr=paste("Chr", hpco$chr, sep=""), chromStart=hpco$pos, chromEnd=hpco$pos+1)	
			write.table(hpcobed, file=paste(t,"hPco",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
			}
		
			#create BED formatted output for lgeo, lpco
			if ( length(lgeo$pos) != 0 ) {
			lgeobed <- data.frame(chr=paste("Chr", lgeo$chr, sep=""), chromStart=lgeo$pos, chromEnd=lgeo$pos+1)	
			write.table(lgeobed, file=paste(t,"lGeo",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
			}
			if ( length(lpco$pos) != 0 ) {
			lpcobed <- data.frame(chr=paste("Chr", lpco$chr, sep=""), chromStart=lpco$pos, chromEnd=lpco$pos+1)	
			write.table(lpcobed, file=paste(t,"lPco",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
			}
		
		#plot to other formats
		#dev.copy(jpeg,file=paste(t,"M+GenomeScan",pp,".jpg", sep=""))
		#dev.copy(png,file=paste(t,"M+GenomeScan",pp,".png", sep=""))
		#dev.copy2pdf(file=paste(t,"M+GenomeScan",pp,".pdf", sep=""))
		dev.off()
	}
}
**********ENDR**********

	Determine the identity of each important SNP

	Download Populus trichocarpa annotation v2.2
wget http://www.plantgdb.org/download/Download/xGDB/PtGDB/Ptrichocarpa_156_gene.gff3.bz2
bzip -d Ptrichocarpa_156_gene.gff3.bz2

	Convert to BED format using bedops, extract only genes from GFF file
cat Ptrichocarpa_156_gene.gff3 | gff2bed | grep -w gene > Pt22genes.bed

	Create background data file, a bed file containing only those genes for which there are SNPs.
	Create files listing genes enriched by M+ analyses for downstream GO over-representation analysis.
	Uses bedops. For whatever reason, this script has to be run in two parts or it fails.
**********BASH**********
***Part 1***
i="Pop"; #taxon identifier
ag="Pt22genes.bed"; #bed file containing all annotated genes for given genome release

#Create a bed file for each SNP locus where we have a measurement (only needs to be done for geo or pco, since SNPs are the same):
chr=$(cut -d' ' -f2 "+"$i"SUMEnrichAcrossBlocks.RgeoTgen.txt" | tail -n +2 | sed 's/^/scaffold_/g');
pos=$(cut -d' ' -f3 "+"$i"SUMEnrichAcrossBlocks.RgeoTgen.txt" | tail -n +2);
pos2=$(awk '{$1 = $1 + 1; print}' <<<"$pos"); 
paste -d$'\t' <(echo "$chr") <(echo "$pos") <(echo "$pos2") > $i"AllSNPs.bed";

#Sort to make sure it will work for bedops --element-of
/usr/local/bin/sort-bed $i"AllSNPs.bed" > tmp.txt;
mv tmp.txt $i"AllSNPs.bed";

#Calculate the intersection of PopAllSNPs.bed and PT22genes.bed, the list of genes for which there are SNPs in the study (the background):
/usr/local/bin/bedops --element-of 1 "$ag" $i"AllSNPs.bed" > $i"AllGenesWithSNPs.bed";
  
#Make background bed files into a simple gene list.
cut -d$'\t' -f4 $i"AllGenesWithSNPs.bed" | sort -u > $i"BackgroundGenesList.txt";
#clean up
rm $i"AllSNPs.bed";
rm $i"AllGenesWithSNPs.bed";

***Part 2***
#Compute lists of genes enriched by M+
for k in "Geo" "Pco";
  do echo $k;
    #Verify that bed files are sorted correctly using bedops:sort-bed, e.g.
    for j in "0.001" "0.05";
      do echo "  $j";
        for h in "h" "l";
          do echo "    $h";
            /usr/local/bin/sort-bed $i$h$k$j.bed | sed 's/^Chr/scaffold_/g'> tmp.txt;
            mv tmp.txt $i$h$k$j.bed;
    
            #Determine nearest gene to highly enriched snps using closest-features
            /usr/local/bin/closest-features --closest --dist $i$h$k$j.bed "$ag" > $i$h"$k"EnrichedGenes$j.bed;
    
            #Extract hits that are within genes, not just close
            grep "|0$" $i$h"$k"EnrichedGenes$j.bed > $i$h"$k"EnrichedWithinGenes$j.bed;
            rm $i$h"$k"EnrichedGenes$j.bed;
    
            #Get gene identifiers for downstream GO analysis
            cut -d$'\t' -f6 $i$h"$k"EnrichedWithinGenes$j.bed | sort -u > $i$h"$k"EnrichedGenesList$j.txt;
            rm $i$h"$k"EnrichedWithinGenes$j.bed;
          done;
      done;
  done;
**********ENDBASH**********

	Use Python GOATOOLS to analyze over-representation:
easy_install goatools;
easy_install fisher;
easy_install statsmodels;
wget http://geneontology.org/ontology/go-basic.obo;
wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo;

	Download a raw gene association file (GAF) from http://amigo.geneontology.org/amigo/search/annotation
	Choose species in 'Species' menu, remove 'not', 'contributes_to', and "colocalizes_with'
	in 'Annotation Qualifier' menu.  Call it PopRawGAF.txt. This GAF will include GO terms for 
	Biological Process, Molecular Function, and Cellular Component ontologies.
	Process the file into an "association" input file, called PopGAF.txt, for GOATOOLS script find_enrichment.py
**********BASH**********
t=Pop;
cut -d$'\t' -f3 "$t"RawGAF.txt | sort -u > GAFgenestmp.txt;
raw=$(cut -d$'\t' -f3,5 "$t"RawGAF.txt);
>tmp.txt;
while read -r line;
do
    f="";
    f=$(grep ^"$line"$'\t' <<<"$raw" | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" >> tmp.txt;
    fi;
done < GAFgenestmp.txt;
sed 's/;$//g' tmp.txt > "$t"GAF.txt;
rm GAFgenestmp.txt;
rm tmp.txt;
**********ENDBASH**********
	Takes 20 minutes or so on glitch.

	Run GOATOOLS find_enrichment.py to get the enrichment status, over (e), or under (p),
	representation (column 3 in output).
	The --no_propagate_counts selects only the least inclusive GO term, i.e no parent terms.
	This eliminates multiple significant results along a parent-child path, but has the
	undesirable consequence of making significant under-representation (p) a questionable
	result.

	For Populus only, need to add a 'g' to the end of each line in source and population input
	files.  These are the '...EnrichedGenesList..." files.
for p in 0.05 0.001;
  do echo $p;
  for e in Geo Pco;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          sed 's/$/g/g' Pop"$h""$e"EnrichedGenesList"$p".txt > tmp.txt
          mv tmp.txt Pop"$h""$e"EnrichedGenesList"$p".txt;
        done;
    done;
  done;
sed 's/$/g/g' PopBackgroundGenesList.txt > tmp.txt
mv tmp.txt PopBackgroundGenesList.txt;

	Move the relevant items (inc. *GAF.txt, go-basic.obo) to a nested folder, "statistical tests", then run
	the GOATOOLS script find_enrichment.py:
**********BASH**********
mv *EnrichedGenesList* "statistical tests";
mv PopBackgroundGenesList.txt "statistical tests";
cd "statistical tests";

t=Pop;
comm="--alpha 0.05 --pval 0.05 --obo go-basic.obo --no_propagate_counts --method holm";

for p in 0.05 0.001;
  do echo $p;
  for e in Geo Pco;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          find_enrichment.py $(echo "$comm" "$t""$h""$e"EnrichedGenesList"$p".txt "$t"BackgroundGenesList.txt "$t"GAF.txt) > +"$t""$h""$e""$p"GOAout.txt;
        done;
    done;
  done;
**********ENDBASH**********

	Examine first part of output files to determine if any GO terms are significantly over-
	or under-represented.
head -20 +*GOA*;

	FINAL RESULT:
There are no significantly over- or under-represented GO terms for Populus among genes that
are well-collected or poorly-collected by geographic or environmental information. True for:
	Rcut
	Rcut2
	Hclust--above also true here for maxbl=50, others only tested with maxbl=200
	Kmeans



	Make a pie chart showing GOslim category representation for genes enriched by M+ at the
	0.05 level. Do this in a folder called "pie charts".

	Extract the three GO categories (biological_process, cellular_component, molecular_function)
	from goslim_plant.obo.  Use PERL modules go-perl > go-filter-subset.pl.

for i in biological_process cellular_component molecular_function;
  do
    /Users/shrub/perl5/bin/go-filter-subset.pl -namespace "$i" goslim_plant.obo > goslim_plant_"$i".obo;
  done;

	For whatever reason (actually, it is because go-filter-subset.pl includes 'part-of's, not
	just 'is_a's as members of a category), two "cellular_component"s remain in goslim_plant_biological_process.obo.
	Ten "biological_process"s remain in goslim_plant_molecular_function.obo.
	Remove them manually.
	
	Make a gene association file (GAF) for M+ well-collected and poorly-collected genes.
	In files: Pop[h,l][Geo,Pco]EnrichedGenesList0.05.txt

**********BASHPARALLEL**********
myp() {
    line=$1;
    o=$(mktemp /tmp/tmp.XXXXXXXX);
    f="";
    f=$(grep ^"$line"$'\t' in.txt | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" > $o;
    fi;
}
export -f myp;

t=Pop;
rm /tmp/tmp.*; #clear tmp directory of tmp. files generated by this script
cut -d$'\t' -f3,5 "$t"RawGAF.txt > in.txt; #make a simplified GAF file of all genes to query
  for e in Geo Pco;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          cat $t$h$e"EnrichedGenesList0.05.txt" | parallel --env myp myp;

          cat /tmp/tmp.* | sort | sed 's/;$//g' > ./$t$h$e"0.05GAF.txt";
          rm /tmp/tmp.*;
        done;
    done;
rm in.txt;
**********ENDBASH**********

	Count the number of occurrences of each GO slim term in the gene association files
	derived from the EnrichedGenesList(s). These are the files Pop[h,l][Pco,Geo]0.05GAF.txt.
	Do this for each major GO category (bp, cc, mf). Output files are like Pop[h,l][Pco,Geo]piechart_[bp,cc,mf].txt.
	Goal here is to calculate the frequency of the function in the well/poorly collected genes and
	compare that to the frequency of the same functions when all annotated genes are considered.
**********BASH**********
  
t="Pop";
b=$(wc -l "$t"GAF.txt | awk '{print $1}'); #the total number of annotated genes
for j in biological_process cellular_component molecular_function;
  do echo $j;
    map_to_slim.py --slim_out=direct --association_file=$t"GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t"slim_"$j".txt; #get the goslim_plant terms associated with all annotated genes using goatools map_to_slim.py
    a=$(awk -F$'\t' '{print $2}' "$t"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for all collected genes.

    for e in Geo Pco;
      do echo "  $e";
        for h in "h" "l";
          do echo "    $h";
            map_to_slim.py --slim_out=direct --association_file=$t$h$e"0.05GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t$h$e"slim_"$j".txt; #get the goslim_plant terms associated with the genes using goatools map_to_slim.py

            go=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | tr ";" "\n" | sort -u); #list of unique GO terms for well/poorly collected genes
            o=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for well/poorly collected genes.
            n=$(wc -l "$t$h$e"slim_"$j".txt | awk '{print $1}'); #number of lines = number of genes for well/poorly collected genes
            n=$(( $n - 4 )); #subtract off header lines
            echo "GOterm ObsTotNumGenes ObsTotNumFunctionHits ObsGOcount obsfreq TotNumAnnotatedGenes ExpTotNumFuncHits ExpGOcount expfreq ratio absratio plratio diff absdiff" > Piechart"$t$h$e"_"$j".txt; #for each unique GO term, count how many times it occurs in slim.txt, divide by number of genes to get the proportion of genes representing that function in the M+ enriched/purified gene set.

            for i in $go;
              do f=$(grep $i "$t$h$e"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for well/poorly collected genes
                x=$(grep $i "$t"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for all genes
                g=$(echo "scale=4;$f/$o" | bc); #freq of function among all enriched functions
                y=$(echo "scale=4;$x/$a" | bc); #freq of function among all functions
                z=$(echo "$g - $y" | bc); #diff btw freq of function among all functions, and among enriched functions (negative means less common in enriched subset)
                zz=${z#-}; #absolute value of z
                rrx=$(echo "scale=4;$g/$y" | bc); #for log2 plotting
                #use below 4 lines for non-log2 axes
                if (( $(bc <<< "$g >= $y") )); 
                  then rr=$(echo "scale=4;$g/$y" | bc); #over represented results in positive fold enrichment
                  else rr=$(echo "scale=4;-$y/$g" | bc); #under represented results in negative fold enrichment
                fi; #ratio of frequency in enriched regions / frequency in all regions
                rrr=${rr#-}; #absolute value of rr
                 echo "$i $n $o $f $g $b $a $x $y $rr $rrr $rrx $z $zz" >> Piechart"$t$h$e"_"$j".txt;
              done;
              
            rm "$t$h$e"slim_"$j".txt;
          done;
      done;
  done;
**********ENDBASH**********

	Get verbal description of slimmed GO terms for eventual plotting. Exclude the cellular_component
	as it is just the location the protein is found. Filter out singletons, where there is a unique
	observation of a GO term in the set of GO terms from enriched regions:

**********BASH**********
t="Pop";
for e in Geo Pco;
  do echo "  $e";
    for h in "h" "l";
      do echo "    $h";
        > o.txt; #freq output
        > p.txt; #ratio output
        for j in biological_process molecular_function;
          do echo $j;
            #filter out singletons, GO terms with only 1 observation in the enriched regions
            sin=$(grep -v "GO:"[0-9]*" "[0-9]*" "[0-9]*" "1" " Piechart$t$h$e"_"$j.txt | tail -n +2);

            #extract columns holding freq and ratio metrics for reformatting
            f=$(echo "$sin" | cut -d' ' -f1); #get GO slim terms
            p=$(echo "$sin" | cut -d' ' -f13); #get diff in freq of function btw enriched functions and all functions
            pp=$(echo "$sin" | cut -d' ' -f14); #get abs of diff in freq of function btw enriched functions and all functions
            q=$(echo "$sin" | cut -d' ' -f10); #get ratio between enriched/all
            qq=$(echo "$sin" | cut -d' ' -f11); #get abs of ratio between enriched/all
            pl=$(echo "$sin" | cut -d' ' -f12); #get ratio between enriched/all to be plotted on log2 axis
            c=$(for i in $f;
              do sed -n -e '/id: '$i'/,$p' goslim_plant_"$j".obo | head -2 | tail -1 | cut -d' ' -f2-;
              done;);
            d=$(for i in $f; do echo "$j"; done;); #create a header column with bp, cc, mf category
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$p")  <(echo "$pp") >> o.txt; #write to temporary freq output file
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$q")  <(echo "$qq") <(echo "$pl")>> p.txt; #write to temporary ratio output file
         done;
          echo Gocat$'\t'Goterm$'\t'name$'\t'diff$'\t'absdiff > PlotRpieFREQ$t$h$e.txt; #make header for freq output file
          echo Gocat$'\t'Goterm$'\t'name$'\t'ratio$'\t'absratio$'\t'plratio > PlotRpieRATIO$t$h$e.txt; #make header for ratio output file
          sort -t$'\t' -nr -k5,5 o.txt >> PlotRpieFREQ$t$h$e.txt;
          sort -t$'\t' -nr -k5,5 p.txt >> PlotRpieRATIO$t$h$e.txt;
         rm o.txt p.txt;
      done;
  done;

**********ENDBASH**********
	
	Use R to construct charts, change hclust to kmeans as necessary in setwd:
**********RSCRIPT**********
#install.packages("ggplot2")
#library(ggplot2)
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Populus/Experiment2/+revision\ with\ hclust/+results/hclust/+PopHclustGO\ work/ToMaxBL50/pie\ charts")

t<-"Pop"
rvals<-c("RATIO", "FREQ")
cvals<-c("h", "l")
evals<-c("Geo","Pco")
for (rr in rvals)
{
  for (cc in cvals)
  {
    for (ee in evals)
    {
      pdf(file=paste("GOchart",rr,t,cc,ee,".pdf", sep=""))
      #par(mfrow=c(1,2))
  
      g <- read.table(paste("PlotRpie",rr,t,cc,ee,".txt", sep=""), header=TRUE, sep="\t")
      h <- g[g$Goterm!="GO:0008150",] #Remove overarching GO categories biological_process, molecular_function
      h <- h[h$Goterm!="GO:0003674",]
  
      par(mar = c(4,20,4,2) + 0.1)
      if(rr == "FREQ") {
        barplot(h$diff, main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$diff),max(h$diff)))
        } 
      else {
        #barplot(h$plratio, log="x", main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$plratio),max(h$plratio)))
        hh=h #transfer h table to a new variable, which will be modified
        hh$id=paste(hh$Goterm,hh$name,sep=" ") #assemble label term
        hh$id <- factor(hh$id, levels=hh$id) #fix order by making id a factor with levels

        ggout<-ggplot(hh,aes(id,plratio,width=0.8)) + 
        geom_bar(stat="identity",color="black",fill="dark grey",size=0.25) + 
        scale_y_continuous(trans='log2', breaks=c(0.125,0.25,0.5,1,2,4,8)) + 
        coord_flip() + 
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text=element_text(size=12, family="ArialMT"))
        
        print(ggout)
        }      
      dev.off()
    }
  }
}
**********ENDR**********





























***SORGHUM***
	Create new haplotypista input, excluding the outlying Sorghum population in India.  This
	causes problems no matter what kind of binning is done. This work is done in folder 
	/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/remakeSordata
	
	Remove population 23 from SorLatLongsp4MOD.txt and save as SorLatLongsp4MOD2.txt.
	Remove population 23 from [geo,env]DIVAtable.txt and save as [geo,env]DIVAtableMOD.txt.
	
	In the folder ../PCO analysis, calculate PCO values with the new set of 22 populations.
	The first four coordinates now account for 89.0517% of the variation, before, with 23
	it was 88.2995%. Make file pcoDIVAtableMOD.txt with the new values.
	

	Modify SorghumHaplotypistaInput.txt, the haplotypista input file, to contain the same names as, and be in the same
	order as SorLatLongsp4MOD.txt.  SorghumHaplotypistaInput.txt contains phased data in two rows per sample with
	missing data imputed.
	
a=$(tail -n +2 SorLatLongsp4MOD.txt | cut -d$'\t' -f1); #list of sample names with popid as prefix
b=$(tail -n +2 SorLatLongsp4MOD.txt | cut -d$'\t' -f1 | sed 's/^[0-9]*.//'); #list of sample names without popid as prefix
c=$(paste -d' ' <(echo "$a") <(echo "$b"));
cp SorghumHaplotypistaInput.txt shitmp.txt; #make a temporary copy to modify
while read -r line; #substitute the popid-prefixed sample name for the bare sample name
do 
    s=$(echo $line | cut -d' ' -f2);
    r=$(echo $line | cut -d' ' -f1);
    echo "$s";
    sed 's/'"$s"'/'"$r"'/g' shitmp.txt > sstmp.txt;
    mv sstmp.txt shitmp.txt;
done <<<"$c";
#sort the temp haplotypista input file the same as SorLatLongsp4MOD.txt
header=$(head -3 shitmp.txt);
tail +4 shitmp.txt | sort -t'.' -n -k1,1 -k2,2 > sstmp.txt;
#verify that sort order is the same
e=$(paste -d' ' <(echo "$a") <(echo "$a") | tr " " "\n"); #double list of sample names with popid as prefix from SorLatLongsp4MOD.txt
d=$(cut -d' ' -f1 sstmp.txt);
diff <(echo "$e") <(echo "$d"); <- hooray, they are the same!  
#sort worked so rename as SorghumHaplotypistaInputtmp.txt
echo "$header" > SorghumHaplotypistaInputtmp.txt;
cat sstmp.txt >> SorghumHaplotypistaInputtmp.txt;
rm shitmp.txt;


	Remove the lines in SorghumHaplotypistaInputtmp.txt that are associated with population 23.  The individuals in
	population 23 can be found in SorLatLongsp4MOD.txt and are:
		23.IS5667.D0D0HACXX.5.250046312
		23.IS6351.628NVAAXX.5.250006919
		23.IS6354.D0D0HACXX.5.250046315
		23.PI533842.C00TDABXX.1.250004174
		23.PI533856.C00TDABXX.1.250004166
		23.PI534021.C00TDABXX.2.250020253
	These are now located in the last 12 lines of the file, so just get rid of those using gnu head:
ghead -n -12 SorghumHaplotypistaInputtmp.txt > SorghumHaplotypistaInputMOD2.txt;

	Verify finally that the order of lines in SorLatLongsp4MOD2.txt and SorghumHaplotypistaInputMOD2.txt is the same
f=$(tail -n +2 SorLatLongsp4MOD2.txt | cut -d$'\t' -f1);
g=$(paste -d' ' <(echo "$f") <(echo "$f") | tr " " "\n"); #double list of sample names with popid as prefix from SorLatLongsp4MOD2.txt
h=$(cut -d' ' -f1 SorghumHaplotypistaInputMOD2.txt | tail -n +4);
diff <(echo "$g") <(echo "$h");  <- yes, they are the same!

	The new haplotypista input file is now complete!!!

	Make a new set of haplotypista output files that contain start and end points for each block.
	Should take ~18hrs. Uses ~120 GB memory, so be cautious, partition among nodes.
	
	on ceres, use 5 short nodes in screen for the commands below:
screen;
sshort;

	Takes ~24 hours.
	Rsync the results back:
rsync -avz --remove-source-files --progress pat.reeves@scinet-login.bioteam.net:"/home/pat.reeves/haplotypista/SorSNP*" .




	on glitch:
ssh compute-0-1
cd haplotypista
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog1.txt -b 1 1 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog2.txt -b 2 2 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog3.txt -b 3 6 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog4.txt -b 7 14 -m ? -p 2 &
ssh compute-0-2
cd haplotypista
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog5.txt -b 15 30 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog6.txt -b 31 60 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog7.txt -b 61 200 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog8.txt -b 201 250 -m ? -p 2 &
ssh compute-0-7
cd haplotypista
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog9.txt -b 251 300 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog10.txt -b 301 350 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog11.txt -b 351 400 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog12.txt -b 401 450 -m ? -p 2 &
ssh compute-0-6
cd haplotypista
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog13.txt -b 451 500 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog14.txt -b 501 575 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog15.txt -b 576 650 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog16.txt -b 651 725 -m ? -p 2 &
ssh compute-0-4
cd haplotypista
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog17.txt -b 726 800 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog18.txt -b 801 875 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog19.txt -b 876 950 -m ? -p 2 &
./haplotypista -i SorghumHaplotypistaInputMOD2.txt -o SorSNP -l SorSNPlog20.txt -b 951 1000 -m ? -p 2 &
	Takes ~24 hours
	
	Undouble the haplotypista output using the following bash script:
**********BASH**********
 MyUndub() {
 b=$1;
 #file root (f) and outfile root (o)
 f=SorSNP;
 o=SorM;
 
 # make doubled version of header to reflect diploid genotypes
    echo $f.b$b;
    z="";
    for i in {1..6}; #6 header lines
      do a=$(head -"$i" $f.b$b | tail -n -1 | tr ' ' '\n'); #extract ith line
        a=$(paste -d' ' <(echo "$a") <(echo "$a") | tr '\n' ' '); #double it side by side, then make it one line
        z=$z"$a
";
     done;
    
   # un-double lines containing genotypes, for each individual
    t4="";
    for i in {1..113}; #half of (total line number minus 6)
      do start=$(expr $i \* 2 + 5);
        t1=$(tail -n +$start $f.b$b | head -1 | tr ' ' '\n');
        t2=$(tail -n +$(expr $start + 1) $f.b$b | head -1 | tr ' ' '\n');
        t3=$(paste -d' ' <(echo "$t1") <(echo "$t2") | tr '\n' ' ');
        t4=$t4"$t3
";
      done;

    # assemble the un-doubled matrix with the doubled header lines
    { echo -n "$z"; head -c -1 <<<"$t4" | cut -d' ' -f2-; } > $o.b$b; #ghead (for mac) trims off white space at end, cut starts at col 2 to skip doubled sample name
 }
export -f MyUndub;
 
seq 1 1000 | parallel --env MyUndub MyUndub;
**********BASHEND**********
	Takes ~10 minutes on ceres short node.


	extract the NG:SS:NS values for each locus from the headers of the SorM.b* files.  Make a table
	called "SorLocusWeights". These values will serve to compare with enrichment values from m+1 of the same loci.

**********BASH**********
f="SorLocusWeights";
g="Sor";
header=$'genomeindex\tblocklength\tblocklengthindex\tlocus\tchr\tchrindex\tNGcount\tSScount\tNScount\tNGpr\tSSpr\tNSpr\tblockstart\tblockend';
> $f.txt;
for ((b=1;b<=1000;b++));
do echo "$g"M.b$b;
  chr=$(head -1 "$g"M.b$b | tr " " "\n" | awk 'NR%2'); #extract line 1 (chr), get odd numbered lines
  loc=$(head -3 "$g"M.b$b | tail -1 | tr " " "\n" | awk 'NR%2'); #extract line 3 (midpt), get odd numbered lines
  st=$(head -5 "$g"M.b$b | tail -1 | tr " " "\n" | awk 'NR%2'); #extract blockstart
  en=$(head -6 "$g"M.b$b | tail -1 | tr " " "\n" | awk 'NR%2'); #extract blockend
  w=$(wc -w <<< $loc); #get the number of loci
  blI=$(seq 1 $w); #compute an index value specific to this blocklength, use jot $w on macos.
  chrloc=$(paste -d'_' <(echo "$chr") <(echo "$loc"));#fuse chr and locus midpoint into locus name
  
  #calculate an index by chromosome, within blocklength
chrI=$(i=1;
       prevc=0;
       for c in $chr;
         do if [ $c -eq $prevc ];
            then i=$(($i + 1));
                echo $i;
            else i=1;
                echo $i;
            fi;
         prevc=$c;
         done;
       )
               
  aacats=$(head -4 "$g"M.b$b | tail -1 | tr " " "\n" | awk 'NR%2' | sed "s/:/"$'\t'"/g"); #extract line 4 (NG:SS:NS), get odd numbered lines
  aapr=$(echo "$aacats" | awk '{sum = $1 + $2 + $3; ng = ($1/sum); ss = ($2/sum); ns = ($3/sum); print ng "\t" ss "\t" ns;}'); #sum NGSSNS counts, and find the proportion of each category

  paste -d$'\t' <(echo "$blI") <(echo "$chrloc") <(echo "$chr") <(echo "$chrI") <(echo "$aacats") <(echo "$aapr") <(echo "$st") <(echo "$en") | sed "s:^:$b"$'\t'":g" >> $f.txt;
done;
awk 'BEGIN { OFS = "\t" } ; {print NR,$0}' $f.txt > tmp.txt; #add a sequential numeric index value (line number) to each row
#nl $f.txt | sed 's/ //g' > tmp.txt; #add a sequential numeric index value (line number) to each row
#cat tmp.txt | sed 's/\./_/g' > tmp.txt; #change . to _ in locus description so r doesn't think its a number
{ echo "$header" ; cat tmp.txt; } > $f.txt;
rm tmp.txt;
**********BASHEND**********
	takes a few minutes.

	make some test plots of the frequency of NG, NS, and SS categories in each haplotype block using
	this data (SorLocusWeights.txt) using R.  Here you use an r script called SorLocusWeightsHeatmap.r.
	This stuff is here: /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+results/+test\ heatmaps

	###REVISED REGION###
	remake the M+ .dat and .var input files with the new geo, env, pco binned variables, calculated 
	with R's 'hclust' function instead of 'cut' or 'cut2':

	Convert haplotypista data sets to M+:
	Run script "2. Convert haplotypista to M+.scpt"  #set so block midpoint positions end up in var file. 
	Use 1 script instance on MacMini (Mavericks).  The script takes
	as input the SorM files, which have diploid data in a single line per individual.  These input files
	are created by the "un-doubling" process.
	Takes ~1 hr.
	
	Add Geo, Env, PCOEnv data to end of M+ files for various optimizations:
	To get new bins using R's 'hclust' function, run script "5b. Calc generalized variance for env data.scpt" with the settings 
	CenterAndScale="bin", cuttype="hclust", and binonly="yes.  Rename the output files like "Sorhclustenvinput7.txt", "Sorkmeanspcoinput7.txt"...


	Then run script "7. Add target data to M+.scpt",
	which will add the new binned values to the .dat and .var files from above.  This was performed 
	with one processor on MacMini (Mavericks).
	Takes ~1.5 hrs on old Script Editor.
	

###

Run m+1 analyses on ceres. This parallelizes on locus, and post processes the m+1 output,
	leaving a SUM file with stats, and a RAW file which is an archive of all the m+1 output, concatenated.
	Structure of runs on Ceres:
	b1, 8dSorb1MOD.sh numslice=4, longmem:                 sbatch --array=1 slurmSorArrayb1longmem.sh
	b2-3, 8dSorb2-3MOD.sh numslice=2, longmem:             sbatch --array=2-3 slurmSorArrayb2-3longmem.sh
	b4-b10, 8dMOD.sh numslice=1, mem, individually:        sbatch --array=4-10 slurmSorArrayb4-10mem.sh
	b11-b25, 8dMOD.sh numslice=1, medium, individually:    sbatch --array=11-25 slurmSorArrayb11-25medium.sh
	b26-b1000, 8dMOD.sh numslice=1, short, in 75 groups of 13: sbatch --array=26-100 slurmSorArrayb26-100short.sh
	
	Split up b26-1000 using:
p=/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/script7outputHclust/;
dp=/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/m+deployment/hclust/;
f=26; #starting blocklength
for j in {1..13}; #counts the number of files to deposit per folder
do for i in {26..100}; #folder range
  do cp "$p"*.b"$f".envgeopco.dat "$dp$i"; #distribute files among folders
    cp "$p"*.b"$f".RgenTenv.var "$dp$i"; #distribute files among folders
    let f++;
  done;
done;

	Split up b1-25 using:
p=/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/script7outputHclust/;
dp=/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/m+deployment/hclust/;
f=1; #starting blocklength
for j in {1..1}; #counts the number of files to deposit per folder
do for i in {1..25};
  do cp "$p"*.b"$f".envgeopco.dat "$dp$i"; #distribute files among folders
    cp "$p"*.b"$f".RgenTenv.var "$dp$i"; #distribute files among folders
    let f++;
  done;
done;

	Compress before transferring to ceres:
tar -czvf sor.tar.gz {1..100};
scp sor.tar.gz pat.reeves@scinet-login.bioteam.net:"/home/pat.reeves/mplusrunsSor/sorb4-1000/kmeans";
tar -xzvf sor.tar.gz;

	Use a scheduled rsync to keep ceres from getting clogged with files:
	hclust and kmeans:
SSHPASS='qwerQWER12#$'; #enter the password as a shell variable
export SSHPASS;
while (true);
do date
  sshpass -e rsync -avz --remove-source-files --progress pat.reeves@scinet-login.bioteam.net:"/home/pat.reeves/mplusrunsSor/sorb4-1000/hclust/*/[R,S][A,U][W,M]*" "/Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/Sorghum/Experiment2/+revision with hclust/+results/hclust/";
  sshpass -e rsync -avz --remove-source-files --progress pat.reeves@scinet-login.bioteam.net:"/home/pat.reeves/mplusrunsSor/sorb4-1000/kmeans/*/[R,S][A,U][W,M]*" "/Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/Sorghum/Experiment2/+revision with hclust/+results/kmeans/";
  echo;
  echo "going to sleep for an hour starting "$(date);
  sleep 3600;
done;



	Run m+1 analyses on glitch. This parallelizes on locus, and post processes the m+1 output,
	leaving a SUM file with stats, and a RAW file which is an archive of all the m+1 output, concatenated.
	Here, b1 was run with numslice=4, procpernode=8. b2-3 were run
	with numslice=2, procpernode=12. b4-b1000 were run with numslice=1, procpernode=24.  Used
	script '8d. Manage M+ runs on clusterSorb1.sh', '...clusterSorb2-b3.sh', '...clusterSorb4-b1000.sh'.
	Script is saved at /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2.
	Following more problems with large analyses (b1,b2) not being repeatable, set up to run analyses on ceres.
	Build a script based on those developed for running At analyses on ceres, start with '8d.sh'.  Modify
	so that numslice=2 (for b2-b3 runs), and procpernode is irrelevant.  Call this script '8dSorghumSlice2RpcoCeres.sh'
	and '8dSorghumSlice2RgeoCeres.sh'.
	
	
	Verify that Rgeo and Rpco SUM files have the same number of lines.  The '8d' script that
	controls the M+ runs occasionally fails.  Verifying that output has the same number of lines is
	a partial check of data integrity. More checks of data integrity at the bottom of this file.
	
for i in {1..1000}; 
do f=$(wc -l SUM.b$i.* | awk '{print $1}' | head -2 | sort -u);
  if [[ $(echo "$f" | wc -l) == 2 ]]; then
    ls -l SUM.b$i.*;
    echo;
  fi;
done;

	Calculate mean M+ enrichment across core sizes, for each locus, for each blocklength. This uses the
	SUMfiles as input and produces the +STATS files. Do this on the cluster.
	
**********BASHPARALLEL-FASTER-RUNINSCREEN**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all core sizes, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum "\t" sum/NF "\t" sqrt(var) "\t" NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

main() {
        l=$1;
        ofile="$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt"; #temporary outfile for each parallel process
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";

        boo=$(grep -F $l $sumfileIN | awk '{ print $13 - $10 }' | tr "\n" " "); #extract all core sizes for this locus, calculate target enrichment
        st=$(mysd "$boo"); #calculate sum, mean, sd, and n using function mysd()
        echo "$b"$'\t'"$l"$'\t'"$st" > $ofile; #write stats to output file
}
export -f main;

#parameters
t="Sorghum";
elist="geo pco";
bmin=1; #min blocklength
bmax=1000;
nnode=8; #number of cluster nodes
blist=$(seq $bmin $bmax); #blocklengthrange
p=$(pwd)"/"; export p;

#parallelize on locus, within loop on blocklength
for e in $elist;
  do echo $e;
    export e;
    for b in $blist;
      do echo "b=$b";
        export b;
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";
        uloc=$(cut -d$'\t' -f2 "$sumfileIN" | tail -n +2 | uniq); #find unique loci from all loci in column 2
        numloc=$(echo "$uloc" | wc -l);
        
        #parallel step
        #sends $gnuN records to each node via parallel step #1, this is piped to parallel step #2, which starts 96 jobs per node from the instruction set it has received
        echo "processing $sumfileIN...";
        gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node
        #echo "$uloc" | parallel --sshloginfile ~/machines --jobs 24 --env main --env mysd --env e --env b --env p main;
        echo "$uloc" | parallel --sshloginfile ~/machines --jobs 1 --env main --env mysd --env e --env b --env p --pipe -N"$gnuN" /home/reevesp/bin/parallel -j96 --env main --env mysd --env e --env b --env p main;

        #consolidate results for each locus into one file
        ofile="+STATS.R"$e"Tgen.b"$b".txt"; #temporary outfile
        > "$p"$ofile;
        for l in $uloc;
          do echo "concatenating loc file $l to +STATS.R"$e"Tgen.b"$b".txt ...";
            cat "$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt" >> "$p"$ofile;
          done;
          #test whether all were concatenated
          #s1=$(wc -l "$p""+STATS.R"$e"Tgen.loc"*".b"$b".txt" | grep total | awk '{print $1}');
          s1=$(ls "$p" | grep -c "loc"); #count the number of files with "loc" in their name using grep. can't wc -l because too many arguments
          s2=$(wc -l "$p"$ofile |  awk '{print $1}');
         if [ $s1 == $s2 ]; 
          then find "$p" -maxdepth 1 -name "*loc*" -delete; #remove files when there are too many for rm
          #rm "$p""+STATS.R"$e"Tgen.loc"*".b"$b".txt"; #remove files if number of lines in concat file is the same as the sum of all input files
          else echo "s1="$s1", s2="$s2". Consolidated file "$ofile" no good. Quitting..." > err.txt;
            exit 1;
        fi;
      done; #$b
      
    #concatenate files for this Rgeo/Rpco
    statsfile="$p""+STATS."$t".R"$e"Tgen.txt";
    echo blocklength$'\t'locus$'\t'sum$'\t'mean$'\t'sd$'\t'n > $statsfile;
    >"$statsfile"TMP;
    for b in $blist;
      do echo "concatenating +STATS.R"$e"Tgen.b"$b".txt ...";
        cat "$p""+STATS.R"$e"Tgen.b"$b".txt" >> "$statsfile"TMP;
      done;
    cat "$statsfile"TMP >> "$statsfile";

    #test whether all individual +STATS files have been added to the final concatenated file using number of lines
    s1=$(wc -l +STATS.R"$e"Tgen.b*.txt | grep total | awk '{print $1}');
    s2=$(wc -l "$statsfile" | awk '{print $1}');
    s2=$(( $s2 - 1 ));
    if [ $s1 == $s2 ]; 
      then rm +STATS.R"$e"Tgen.b*.txt; #remove files if number of lines in concat file is the same as the sum of all input files
        rm "$statsfile"TMP;
      else echo "the number of lines ain't the same. something is 'crewed.";
    fi;

  done; #$e
**********BASHEND**********
	~Takes around 5 days.  Or takes ~10 hrs.  Not sure why the discrepancy.  Possibly you
	need to hide 'screen' so it doesn't have to print steps over network?





	add the proportion of NG,SS,NS sites from SorLocusWeights.txt to the +STATS summary file,
	generating the +v2STATS file.  Requires the SorLocusWeights.txt file as input.  

**********BASH**********
myAddaa() {
  weightsfile=$1;
  e=$2;
  key=$3;
  
  statsfile="+STATS.""$key"".R"$e"Tgen.txt";
  s=$(cut -d$'\t' -f2 $statsfile | tail -n +2 | md5sum); #get md5sum of locus ID column from the +STATS file
  w1=$(cut -d$'\t' -f4 $weightsfile | sed 's/_/./g'); #get the locus ID less blocklength index column, replace _ with . for md5 check
  w2=$(cut -d$'\t' -f3 $weightsfile); #get the blocklength index column
  w=$(paste -d'.' <(echo "$w1") <(echo "$w2") | tail -n +2 | md5sum); #get md5sum of locus ID column from the Weights file
  if [ "$s" == "$w" ]; #verify that files are in same order using a checksum of the unique locus ID column
    then boo=$(sed 's/blocklength/blocklength2/g' $weightsfile | sed 's/locus/locus2/g'); #modify the header in the Weights file so that there are no redundant column names after combining
      paste -d$'\t' $statsfile <(echo "$boo") > "+v2STATS.""$key"".R"$e"Tgen.txt"; #paste the statsfile and the weights file together
    else echo "md5sums do not match, $e, aborting.";
      kill -INT $$; #terminate the script, return to the shell
  fi;
}
export -f myAddaa;

weightsfile="SorLocusWeights.txt";
elist="geo pco";
key="Sorghum";
parallel --env myAddaa myAddaa ::: "$weightsfile" ::: $elist ::: "$key";

**********BASHEND**********
	Takes ~1 minute.
	
#Plot summed enrichment by genomic position across blocklengths as a heatmap. Use Rscript
#SorGenomicGeography.r. Uses the +v2STATS files as input.


	create a file '+SorSNPtoBlockMap.txt' that shows which SNPs fall in which blocks for
	b2-b200.
**********BASHPARALLEL**********
myp() {
  b=$1;
  #define parameters
  nchr=10; #number of chromosomes
  snps=$(cat "$path""$p""snps.tmp.txt");
  
    echo "blocklength=$b";
    bnames=b"$b"; #bnames will contain the haplotype block names associated with each snp for a given blocklength
    bla=$(grep "^[0-9]*"$'\t'"$b"$'\t' "$path""$p"LocusWeights.txt); #get all lines for blocklength b

    #repeat through chromosomes
    for i in $(seq 1 $nchr);
      do echo "  chr""$i";
        snpsizes=$(echo "$snps" | grep ^"$i"_ | sed 's/'$i'_//g'); #get snp sizes for this chromosome and blocklength 1
        
        declare -a cname=( $(echo "$bla" | grep $'\t'"$i_" | awk '{print $4}') ); #locus names for this chr & blocklength
        declare -a cstart=( $(echo "$bla" | grep $'\t'"$i_" | awk '{print $13}') ); #start points for each locus in $cname
        declare -a cend=( $(echo "$bla" | grep $'\t'"$i_" | awk '{print $14}') ); #end points for each locus in $cname
        clen=$(echo ${#cname[@]}); #total length of the arrays for chromosome i
       
        startpt=0;
        for s in $snpsizes;
          do found=0;
            #test whether snp position is found within haplotype block by seeing whether it is both >= start point and <= end point
            for ((c=$startpt;c<$clen;c++)); 
              do if (( $(bc <<< "$s >= ${cstart[$c]} && $s <= ${cend[$c]}") )); then
              #do if (( $(bc <<< "$s >= ${cstart[$c]}") )) && (( $(bc <<< "$s <= ${cend[$c]}") )); then
                  bnames="$bnames"$'\n'${cname[$c]}; #add the haplotype block name to the list that maps them to SNP name
                  found=1;
                  startpt=$c; #make search more efficient by ignoring earlier blocks in future searches
                fi;
                if [[ $found == 1 ]]; then break; fi;
              done;
            #deal with situation where SNP position is not found in a haplotype block (this will
            #sometimes happen at the end of chromosomes)
            if [[ $found == 0 ]]; then
              bnames="$bnames"$'\n'"--"; 
            fi;
          done; #s
      done; #i
      
      echo "$bnames" > "$path""$p""b$b.tmp";
}
export -f myp;

#define parameters
maxbl=1000; #maximum blocklength
p="Sor"; export p;
path=$(pwd)"/"; export path; #path to current folder
a=$(grep "^[0-9]*"$'\t'"1"$'\t' "$path""$p"LocusWeights.txt); #get all lines with blocklength 1
snps=$(echo "$a" | awk '{print $4}'); #get a list of all SNP positions
echo "$snps" > "$path""$p""snps.tmp.txt"; #put in a file accessible by all procs so it doesn't have to be passed by parallel or loaded into the environment
outf=$(head -1 "$path""$p"LocusWeights.txt | cut -d$'\t' -f1-6);
outf="$outf"$'\n'$(echo "$a" | cut -d$'\t' -f1-6);

#begin parallel
seq 2 $maxbl | parallel --sshloginfile ~/machines --jobs 24 --env myp --env path --env p myp;

#reassemble
echo "reassembling...";
comm=""; #create a list of tmp files to paste together
for b in $(seq 2 $maxbl);
  do comm="$comm""$p""b$b.tmp ";
  done;
paste -d$'\t' <(echo "$outf") $comm > "$path""+$p""SNPtoBlockMap.txt"; #consecutively paste results from each blocklength

#clean up
rm "$path""$p"snps.tmp.txt;
for b in $(seq 2 $maxbl);
  do rm "$path""$p""b$b.tmp"; #clean up
  done;
**********BASHEND**********
	Takes ~4 days.

	Using the file +*SNPtoBlockMap.txt, calculate the sum of M+ enrichment value across blocklengths,
	for each SNP, i.e. for each position in the genome.  These data will be extracted from the 
	+v2STATS* files in several steps.  The first makes a file containing a table with the max 
	enrichment value at each site for each block, *EnrichAcrossBlocks*. Step 1:

**********BASHPARALLEL**********
myp() {
  e=$1;
  b=$2;
  v2name="$p""+v2STATS."$v2p".R"$e"Tgen.txt";
        
        echo "e=$e  b=$b";
        bout="$p""rr.R"$e"Tgen.b""$b"".tmp"; #temporary outfile

        #cut out the locus id and sum columns for current blocklenth
        gg=$(grep ^$b$'\t' "$v2name" | cut -d$'\t' -f1-3);
        
        #cut column with locus names for current blocklength from SNPtoBlockMap
        if [ $b = 1 ]; then
          cc=$(cut -d$'\t' -f4 "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        else col=$(( $b + 5 )); #calculate the column to cut out for each blocklength
          cc=$(cut -d$'\t' -f$col "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        fi;

        #search for locus names in $gg, retrieve enrichment value
        rr="b$b";

#put another parallel here if necessary
        for l in $cc;
          do if [[ $l = "--." ]]; then
                r="--"; #test for empty enrichment value
              else
                r=$(grep -m 1 ^"$b"$'\t'"$l" <<<"$gg" | awk '{print $3}');
              fi;
            rr+=$'\n'$r; #add current locus enrichment value to list
          done;

         #write out the result
         echo "$rr" > "$bout";
         #truncate -s -1 "$bout"; #remove the hanging newline, use 'gtruncate' only on osx
}
export -f myp;

#parameters
v2p="Sorghum"; export v2p;
t="Sor";
minbl=1;
maxbl=50; #50, 1000
blist=$(seq $minbl $maxbl);
elist="geo pco";

p=$(pwd)"/"; export p;
tname="$p""+""$t""SNPtoBlockMap.txt"; export tname;

#begin parallel
parallel --sshloginfile ~/machines --jobs 24 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist;
  
#reconstitute
echo "reconstituting..."
for e in $elist;
  do echo "e=$e";
    comm=""; #create a list of tmp files to paste together
    for b in $blist;
      do comm+="$p""rr.R"$e"Tgen.b$b.tmp ";
    done;
    pos=$(cut -d$'\t' -f4 $tname | sed 's/'^.*_'//g' | sed 's/locus/pos/g');
    h=$(cut -d$'\t' -f1-6 $tname | tr $'\t' ' ');
    h2=$(paste -d' ' <(echo "$h") <(echo "$pos"));
    paste -d' ' <(echo "$h2") $comm > "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt";
    rm $comm; #clean up
  done;#e
  
**********BASHEND**********
	Takes ~24 hours.

	Process each line (SNP position) of *EnrichAcrossBlocks*.txt to get sum across all blocklengths
	of cumulative enrichment across all core sizes of max enrichment value from 10 reps of M+.
	Produces a file called SUM*EnrichAcrossBlocks... Step 2:

**********BASH**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all blocklengths, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       #returns space delimited string like "sum mean sd n"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum " " sum/NF " " sqrt(var) " " NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

myp() {
    l=$1

    p=$(echo "$l" | cut -d' ' -f8-); #cut out the portion of the line with max enrichment values

    #skip loci that have missing blocks
    if [[ $p != *"--"* ]]; then
      m=$(mysd "$p"); #calculate sum, mean, sd, n
      chr=$(echo "$l" | cut -d' ' -f5);
      pos=$(echo "$l" | cut -d' ' -f7);
      i=$(echo "$l" | cut -d' ' -f1); #get genome index
      #loc=$chr"."$pos"."$(echo "$l" | cut -d' ' -f6);
      loc=$chr"."$pos"."$i;
      echo "snp=$i";

      echo $i" "$chr" "$pos" "$loc" "$m > "$path"$i".tmp"; #add new values to growing string, write to tmp file
    fi;
}
export -f myp;

#parameters
elist="geo pco";
t="Sor";
path=$(pwd)"/"; export path;

#parallel
for e in $elist;
do echo $e;
  nnode=8;
  numloc=$(wc -l +"$t"EnrichAcrossBlocks.R"$e"Tgen.txt | awk '{print $1}');
  numloc=$(( $numloc - 1 ));
  gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node

  #tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --sshloginfile ~/machines --jobs 24 --env myp --env mysd --env path myp;
  tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --sshloginfile ~/machines --jobs 1 --env myp --env mysd --env path --pipe --round-robin /home/reevesp/bin/parallel -j96  --env myp --env mysd --env path myp;
  
  #assemble results
  echo "concatenating...";
  ofile="+"$t"SUMEnrichAcrossBlocks.R"$e"Tgen.txt";
  o="genomeindex chr pos loc sum mean sd n";
  echo "$o" > $ofile;
  
  find "$path" -name "*.tmp" -print0 | xargs -0 cat | sort -n -t' ' -k1 >> $ofile;
    
  #clean up
  echo "clean up...";
  find "$path" -name "*.tmp" -print0 | xargs -0 rm -f
done;
**********ENDBASH**********
	Takes ~35 minutes


	Use R to make some pretty plots of significant enrichment across the genome

**********BEGINR**********
#install.packages("fitdistrplus")
#install.packages("plotrix")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/+results/hclust/+SorHclustSummaryFiles/ToMaxBL50")
library(fitdistrplus)
library(plotrix)

t<-"Sor"
usezero="no" #this sets whether to include significant sites only when they are also > or < 0
             #to use the zero criterion, say "yes", to just use the ecdf around the mean, say "no"
pvals<-c(0.001,0.05)
for (pp in pvals)
{
	pdf(file=paste(t,"M+GenomeScan",pp,".pdf", sep=""))
	par(mfrow=c(2,1))

		i="geo"
		g <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",i,"Tgen.txt", sep=""), header=TRUE, sep=" ")
		mycdf <- ecdf(g$sum) #calculate the cdf of the empirical distribution
		p <- mycdf(g$sum) #get the probability of each observation (or lower)
		g$pecdf <- p #add new column to g, containing the probability from the empirical cdf
		g$index <- 1:nrow(g) #add a column with the sequential index
		h <- g[which(g$pecdf > 1-pp),] #get snps that are highly enriched
		if ( usezero == "yes" ) {
			hgeo <- h[which(h$sum>0),] #remove snps where sum<0, i.e. less than random
		}
		else {
			hgeo <- h # do not remove snps based on positive or negative
		}
		l <- g[which(g$pecdf < pp),] #get snps that are lowly enriched
		if ( usezero == "yes" ) {
			lgeo <- l[which(l$sum<0),] #remove snps where sum>0, i.e. more than random
		}
		else {
			lgeo <- l # do not remove snps based on positive or negative
		}
		highposgeo<-hgeo$index #get the sequential position of the highly enriched
		lowposgeo<-lgeo$index #get the sequential position of the deficient

		plot(g$index, g$sum, type="l", main=paste(t,i,pp,sep=""))
		color.scale.lines(g$index, g$sum, col=factor(g$chr))
		points(g$index[highposgeo], g$sum[highposgeo], col = "red", cex=0.7) #mark significant values
		points(g$index[lowposgeo], g$sum[lowposgeo], col = "blue", cex=0.7) #mark significant values

		i="pco"
		g <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",i,"Tgen.txt", sep=""), header=TRUE, sep=" ")
		mycdf <- ecdf(g$sum) #calculate the cdf of the empirical distribution
		p <- mycdf(g$sum) #get the probability of each observation (or lower)
		g$pecdf <- p #add new column to g, containing the probability from the empirical cdf
		g$index <- 1:nrow(g) #add a column with the sequential index
		h <- g[which(g$pecdf > 1-pp),] #get snps that are highly enriched
		if ( usezero == "yes" ) {
		hpco <- h[which(h$sum>0),] #remove snps where sum<0, i.e. less than random
		}
		else {
			hpco <- h
		}
		l <- g[which(g$pecdf < pp),] #get snps that are lowly enriched
		if ( usezero == "yes" ) {
		lpco <- l[which(l$sum<0),] #remove snps where sum>0, i.e. more than random
		}
		else {
			lpco <- l
		}
		highpospco<-hpco$index #get the sequential position of the highly enriched
		lowpospco<-lpco$index #get the sequential position of the deficient

		plot(g$index, g$sum, type="l", main=paste(t,i,pp,sep=""))
		color.scale.lines(g$index, g$sum, col=factor(g$chr))
		points(g$index[highpospco], g$sum[highpospco], col = "red", cex=0.7) #mark significant values
		points(g$index[lowpospco], g$sum[lowpospco], col = "blue", cex=0.7) #mark significant values

		#calculate intersection of well-collected significant snps
		isect <- hgeo[is.element(hgeo$pos, intersect(hgeo$pos,hpco$pos)),]
		write.table(isect, file=paste(t,"HighCommonToGeoPco",pp,".txt", sep=""), row.names=FALSE, quote=FALSE)

		#calculate intersection of poorly-collected significant snps
		isect <- lgeo[is.element(lgeo$pos, intersect(lgeo$pos,lpco$pos)),]
		write.table(isect, file=paste(t,"LowCommonToGeoPco",pp,".txt", sep=""), row.names=FALSE, quote=FALSE)
	
		#create BED formatted output for hgeo, hpco. verify first whether there exist any snps that are
		#both significantly enriched and better than random (or significantly suppressed and worse than random.)
		#if not, do not print an output table.
		if ( length(hgeo$pos) != 0 ) {
		hgeobed <- data.frame(chr=paste("Chr", hgeo$chr, sep=""), chromStart=hgeo$pos, chromEnd=hgeo$pos+1)	
		write.table(hgeobed, file=paste(t,"hGeo",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}
		if ( length(hpco$pos) != 0 ) {
		hpcobed <- data.frame(chr=paste("Chr", hpco$chr, sep=""), chromStart=hpco$pos, chromEnd=hpco$pos+1)	
		write.table(hpcobed, file=paste(t,"hPco",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}
		
		#create BED formatted output for lgeo, lpco
		if ( length(lgeo$pos) != 0 ) {
		lgeobed <- data.frame(chr=paste("Chr", lgeo$chr, sep=""), chromStart=lgeo$pos, chromEnd=lgeo$pos+1)	
		write.table(lgeobed, file=paste(t,"lGeo",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}
		if ( length(lpco$pos) != 0 ) {
		lpcobed <- data.frame(chr=paste("Chr", lpco$chr, sep=""), chromStart=lpco$pos, chromEnd=lpco$pos+1)	
		write.table(lpcobed, file=paste(t,"lPco",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}

	dev.off()
}
**********ENDR**********

	Determine the identity of each important Sorghum SNP

	Download Sorghum bicolor annotation v1.4
wget ftp://ftp.jgi-psf.org/pub/JGI_data/Sorghum_bicolor/v1.0/Sbi/annotation/Sbi1.4/Sbi1.4.gff3.gz
gunzip Sbi1.4.gff3.gz

	Convert to BED format using bedops, extract only genes from GFF file
cat Sbi1.4.gff3 | gff2bed | grep -w gene > Sb14genes.bed

	Create background data file, a bed file containing only those genes for which there are SNPs.
	Create files listing genes enriched by M+ analyses for downstream GO over-representation analysis.
	Uses bedops. For whatever reason, this script has to be run in two parts or it fails.
**********BASH**********
***Part 1***
i="Sor"; #taxon identifier
ag="Sb14genes.bed"; #bed file containing all annotated genes for given genome release

#Create a bed file for each SNP locus where we have a measurement (only needs to be done for geo or pco, since SNPs are the same):
chr=$(cut -d' ' -f2 "+"$i"SUMEnrichAcrossBlocks.RgeoTgen.txt" | tail -n +2 | sed 's/^/chromosome_/g');
pos=$(cut -d' ' -f3 "+"$i"SUMEnrichAcrossBlocks.RgeoTgen.txt" | tail -n +2);
pos2=$(awk '{$1 = $1 + 1; print}' <<<"$pos"); 
paste -d$'\t' <(echo "$chr") <(echo "$pos") <(echo "$pos2") > $i"AllSNPs.bed";

#Sort to make sure it will work for bedops --element-of
/usr/local/bin/sort-bed $i"AllSNPs.bed" > tmp.txt;
mv tmp.txt $i"AllSNPs.bed";

#Calculate the intersection of SorAllSNPs.bed and Sb14genes.bed, the list of genes for which there are SNPs in the study (the background):
/usr/local/bin/bedops --element-of 1 "$ag" $i"AllSNPs.bed" > $i"AllGenesWithSNPs.bed";
  
#Make background bed files into a simple gene list.
cut -d$'\t' -f4 $i"AllGenesWithSNPs.bed" | sort -u > $i"BackgroundGenesList.txt";
#clean up
rm $i"AllSNPs.bed";
rm $i"AllGenesWithSNPs.bed";

***Part 2***
#Compute lists of genes enriched by M+
for k in "Geo" "Pco";
  do echo $k;
    #Verify that bed files are sorted correctly using bedops:sort-bed, e.g.
    for j in "0.001" "0.05";
      do echo "  $j";
        for h in "h" "l";
          do echo "    $h";
            /usr/local/bin/sort-bed $i$h$k$j.bed | sed 's/^Chr/chromosome_/g' > tmp.txt;
            mv tmp.txt $i$h$k$j.bed;
    
            #Determine nearest gene to highly enriched snps using closest-features
            /usr/local/bin/closest-features --closest --dist $i$h$k$j.bed "$ag" > $i$h"$k"EnrichedGenes$j.bed;
    
            #Extract hits that are within genes, not just close
            grep "|0$" $i$h"$k"EnrichedGenes$j.bed > $i$h"$k"EnrichedWithinGenes$j.bed;
            rm $i$h"$k"EnrichedGenes$j.bed;
    
            #Get gene identifiers for downstream GO analysis
            cut -d$'\t' -f6 $i$h"$k"EnrichedWithinGenes$j.bed | sort -u > $i$h"$k"EnrichedGenesList$j.txt;
            rm $i$h"$k"EnrichedWithinGenes$j.bed;
          done;
      done;
  done;
**********ENDBASH**********

	Use Python GOATOOLS to analyze over-representation:
easy_install goatools;
easy_install fisher;
easy_install statsmodels;
wget http://geneontology.org/ontology/go-basic.obo;
wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo;

	Download a raw gene association file (GAF) from http://amigo.geneontology.org/amigo/search/annotation
	Choose species in 'Species' menu, remove 'not', 'contributes_to', and "colocalizes_with'
	in 'Annotation Qualifier' menu.  Call it SorRawGAF.txt. This GAF will include GO terms for 
	Biological Process, Molecular Function, and Cellular Component ontologies.
	Process the file into an "association" input file, called SorGAF.txt, for GOATOOLS script find_enrichment.py
**********BASH**********
t=Sor;
cut -d$'\t' -f3 "$t"RawGAF.txt | sort -u > GAFgenestmp.txt;
raw=$(cut -d$'\t' -f3,5 "$t"RawGAF.txt);
>tmp.txt;
while read -r line;
do
    f="";
    f=$(grep ^"$line"$'\t' <<<"$raw" | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" >> tmp.txt;
    fi;
done < GAFgenestmp.txt;
sed 's/;$//g' tmp.txt > "$t"GAF.txt;
rm GAFgenestmp.txt;
rm tmp.txt;
**********ENDBASH**********
	Takes a few hours.

	Run GOATOOLS find_enrichment.py to get the enrichment status, over (e), or under (p),
	representation (column 3 in output).
	The --no_propagate_counts selects only the least inclusive GO term, i.e no parent terms.
	This eliminates multiple significant results along a parent-child path, but has the
	undesirable consequence of making significant under-representation (p) a questionable
	result.

	Move the relevant items (inc. *GAF.txt, go-basic.obo) to a nested folder, "statistical tests", then run
	the GOATOOLS script find_enrichment.py:
**********BASH**********
mv *EnrichedGenesList* "statistical tests";
mv SorBackgroundGenesList.txt "statistical tests";
cd "statistical tests";

t=Sor;
comm="--alpha 0.05 --pval 0.05 --obo go-basic.obo --no_propagate_counts --method holm";
for p in 0.05 0.001;
  do echo $p;
  for e in Geo Pco;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          find_enrichment.py $(echo "$comm" "$t""$h""$e"EnrichedGenesList"$p".txt "$t"BackgroundGenesList.txt "$t"GAF.txt) > +"$t""$h""$e""$p"GOAout.txt;
        done;
    done;
  done;
**********ENDBASH**********

	Examine first part of output files to determine if any GO terms are significantly over-
	or under-represented.
head -20 +*GOA*;

	FINAL RESULT:
For MaxBL=50:
	There are no significantly over- or under-represented GO terms among Sorghum genes that are well-
	or poorly- collected using either geographic or environmental data.

For MaxBL=1000:
	Using Rcut:
	There is one significantly over-represented GO term among Sorghum genes that are poorly collected
	using environmental data. There are no significantly over- or under-represented GO terms among 
	Sorghum genes that are well- or poorly- collected using geographic data and Rcut.
	==> +SorlPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0000054	BP	e	ribosomal subunit export from nucleus	5/1160	8/23278	1.5e-05	n.a.	5	0.0377

	Using Rcut2:
	There are two significantly over-represented GO terms among Sorghum genes that are poorly collected
	using environmental data. There are no significantly over- or under-represented GO terms among 
	Sorghum genes that are well- or poorly- collected using geographic data and Rcut2.
	==> +SorlPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0009407	BP	e	toxin catabolic process	20/1365	46/23278	3.52e-11	n.a.	20	8.81e-08
	GO:0006749	BP	e	glutathione metabolic process	23/1365	83/23278	2.38e-10	n.a.	23	5.95e-07

	Using Hclust:
	There are no significantly over- or under-represented GO terms among Sorghum genes that are well-
	or poorly- collected using either geographic or environmental data.
 
	Using Kmeans:
	There is one significantly over-represented GO term among Sorghum genes that are poorly collected 
	using environmental data.
	==> +SorlPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0009407	BP	e	toxin catabolic process	11/1209	46/23278	1.76e-05	n.a.	11	0.0441
	There are no significantly over- or under-represented GO terms for Sorghum among genes that
	are well- or poorly-collected using geographic information and Kmeans.




	Make a pie chart showing GOslim category representation for genes enriched by M+ at the
	0.05 level.

	Extract the three GO categories (biological_process, cellular_component, molecular_function)
	from goslim_plant.obo.  Use PERL modules go-perl > go-filter-subset.pl.

for i in biological_process cellular_component molecular_function;
  do
    /Users/shrub/perl5/bin/go-filter-subset.pl -namespace "$i" goslim_plant.obo > goslim_plant_"$i".obo;
  done;

	For whatever reason (actually, it is because go-filter-subset.pl includes 'part-of's, not
	just 'is_a's as members of a category), two "cellular_component"s remain in goslim_plant_biological_process.obo.
	Ten "biological_process"s remain in goslim_plant_molecular_function.obo.
	Remove them manually.
	
	Make a gene association file (GAF) for M+ well-collected and poorly-collected genes.
	Input files: Sor[h,l][Geo,Pco]EnrichedGenesList0.05.txt

**********BASHPARALLEL**********
myp() {
    line=$1;
    o=$(mktemp /tmp/tmp.XXXXXXXX);
    f="";
    f=$(grep ^"$line"$'\t' in.txt | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" > $o;
    fi;
}
export -f myp;

t=Sor;
rm /tmp/tmp.*; #clear tmp directory of tmp. files generated by this script
cut -d$'\t' -f3,5 "$t"RawGAF.txt > in.txt; #make a simplified GAF file of all genes to query
  for e in Geo Pco;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          cat $t$h$e"EnrichedGenesList0.05.txt" | parallel --env myp myp;

          cat /tmp/tmp.* | sort | sed 's/;$//g' > ./$t$h$e"0.05GAF.txt";
          rm /tmp/tmp.*;
        done;
    done;
rm in.txt;
**********ENDBASH**********

	Count the number of occurrences of each GO slim term in the gene association files
	derived from the EnrichedGenesList(s). These are the files Sor[h,l][Pco,Geo]0.05GAF.txt.
	Do this for each major GO category (bp, cc, mf). Output files are like PiechartSor[h,l][Pco,Geo]_[bp,cc,mf].txt.
	Goal here is to calculate the frequency of the function in the well/poorly collected genes and
	compare that to the frequency of the same functions when all annotated genes are considered.
	Requires go-basic.obo and SorGAF.txt as input.
**********BASH**********
  
t="Sor";
b=$(wc -l "$t"GAF.txt | awk '{print $1}'); #the total number of annotated genes
for j in biological_process cellular_component molecular_function;
  do echo $j;
    map_to_slim.py --slim_out=direct --association_file=$t"GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t"slim_"$j".txt; #get the goslim_plant terms associated with all annotated genes using goatools map_to_slim.py
    a=$(awk -F$'\t' '{print $2}' "$t"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for all collected genes.

    for e in Geo Pco;
      do echo "  $e";
        for h in "h" "l";
          do echo "    $h";
            map_to_slim.py --slim_out=direct --association_file=$t$h$e"0.05GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t$h$e"slim_"$j".txt; #get the goslim_plant terms associated with the genes using goatools map_to_slim.py

            go=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | tr ";" "\n" | sort -u); #list of unique GO terms for well/poorly collected genes
            o=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for well/poorly collected genes.
            n=$(wc -l "$t$h$e"slim_"$j".txt | awk '{print $1}'); #number of lines = number of genes for well/poorly collected genes
            n=$(( $n - 4 )); #subtract off header lines
            echo "GOterm ObsTotNumGenes ObsTotNumFunctionHits ObsGOcount obsfreq TotNumAnnotatedGenes ExpTotNumFuncHits ExpGOcount expfreq ratio absratio plratio diff absdiff" > Piechart"$t$h$e"_"$j".txt; #for each unique GO term, count how many times it occurs in slim.txt, divide by number of genes to get the proportion of genes representing that function in the M+ enriched/purified gene set.

            for i in $go;
              do f=$(grep $i "$t$h$e"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for well/poorly collected genes
                x=$(grep $i "$t"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for all genes
                g=$(echo "scale=4;$f/$o" | bc); #freq of function among all enriched functions
                y=$(echo "scale=4;$x/$a" | bc); #freq of function among all functions
                z=$(echo "$g - $y" | bc); #diff btw freq of function among all functions, and among enriched functions (negative means less common in enriched subset)
                zz=${z#-}; #absolute value of z
                rrx=$(echo "scale=4;$g/$y" | bc); #for log2 plotting
                #use below 4 lines for non-log2 axes
                if (( $(bc <<< "$g >= $y") )); 
                  then rr=$(echo "scale=4;$g/$y" | bc); #over represented results in positive fold enrichment
                  else rr=$(echo "scale=4;-$y/$g" | bc); #under represented results in negative fold enrichment
                fi; #ratio of frequency in enriched regions / frequency in all regions
                rrr=${rr#-}; #absolute value of rr
                 echo "$i $n $o $f $g $b $a $x $y $rr $rrr $rrx $z $zz" >> Piechart"$t$h$e"_"$j".txt;
              done;
              
            rm "$t$h$e"slim_"$j".txt;
          done;
      done;
  done;
**********ENDBASH**********

	Get verbal description of slimmed GO terms for eventual plotting. Exclude the cellular_component
	as it is just the location the protein is found:
t="Sor";
for e in Geo Pco;
  do echo "  $e";
    for h in "h" "l";
      do echo "    $h";
        > o.txt; #freq output
        > p.txt; #ratio output
        for j in biological_process molecular_function;
          do echo $j;
            #filter out singletons, GO terms with only 1 observation in the enriched regions
            sin=$(grep -v "GO:"[0-9]*" "[0-9]*" "[0-9]*" "1" " Piechart$t$h$e"_"$j.txt | tail -n +2);

            #extract columns holding freq and ratio metrics for reformatting
            f=$(echo "$sin" | cut -d' ' -f1); #get GO slim terms
            p=$(echo "$sin" | cut -d' ' -f13); #get diff in freq of function btw enriched functions and all functions
            pp=$(echo "$sin" | cut -d' ' -f14); #get abs of diff in freq of function btw enriched functions and all functions
            q=$(echo "$sin" | cut -d' ' -f10); #get ratio between enriched/all
            qq=$(echo "$sin" | cut -d' ' -f11); #get abs of ratio between enriched/all
            pl=$(echo "$sin" | cut -d' ' -f12); #get ratio between enriched/all to be plotted on log2 axis
            c=$(for i in $f;
              do sed -n -e '/id: '$i'/,$p' goslim_plant_"$j".obo | head -2 | tail -1 | cut -d' ' -f2-;
              done;);
            d=$(for i in $f; do echo "$j"; done;); #create a header column with bp, cc, mf category
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$p")  <(echo "$pp") >> o.txt; #write to temporary freq output file
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$q")  <(echo "$qq") <(echo "$pl")>> p.txt; #write to temporary ratio output file
         done;
          echo Gocat$'\t'Goterm$'\t'name$'\t'diff$'\t'absdiff > PlotRpieFREQ$t$h$e.txt; #make header for freq output file
          echo Gocat$'\t'Goterm$'\t'name$'\t'ratio$'\t'absratio$'\t'plratio > PlotRpieRATIO$t$h$e.txt; #make header for ratio output file
          sort -t$'\t' -nr -k5,5 o.txt >> PlotRpieFREQ$t$h$e.txt;
          sort -t$'\t' -nr -k5,5 p.txt >> PlotRpieRATIO$t$h$e.txt;
         rm o.txt p.txt;
      done;
  done;
	
	Use R to construct charts:
**********RSCRIPT**********
#install.packages("ggplot2")
#library(ggplot2)
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/+results/hclust/+SorHclustGO\ work/ToMaxBL50/pie\ charts")

t<-"Sor"
rvals<-c("RATIO", "FREQ")
cvals<-c("h", "l")
evals<-c("Geo","Pco")
for (rr in rvals)
{
  for (cc in cvals)
  {
    for (ee in evals)
    {
      pdf(file=paste("GOchart",rr,t,cc,ee,".pdf", sep=""))
      #par(mfrow=c(1,2))
  
      g <- read.table(paste("PlotRpie",rr,t,cc,ee,".txt", sep=""), header=TRUE, sep="\t")
      h <- g[g$Goterm!="GO:0008150",] #Remove overarching GO categories biological_process, molecular_function
      h <- h[h$Goterm!="GO:0003674",]
  
      par(mar = c(4,20,4,2) + 0.1)
      if(rr == "FREQ") {
        barplot(h$diff, main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$diff),max(h$diff)))
        } 
      else {
        #barplot(h$plratio, log="x", main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$plratio),max(h$plratio)))
        hh=h #transfer h table to a new variable, which will be modified
        hh$id=paste(hh$Goterm,hh$name,sep=" ") #assemble label term
        hh$id <- factor(hh$id, levels=hh$id) #fix order by making id a factor with levels

        ggout<-ggplot(hh,aes(id,plratio,width=0.8)) + 
        geom_bar(stat="identity",color="black",fill="dark grey",size=0.25) + 
        scale_y_continuous(trans='log2', breaks=c(0.125,0.25,0.5,1,2,4,8)) + 
        coord_flip() + 
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text=element_text(size=12, family="ArialMT"))
        
        print(ggout)
        }      
      dev.off()
    }
  }
}
**********ENDR**********



















***ARABIDOPSIS***
	Modify the Arabidopsis lat/long file AtLatLongsp4MOD.txt to have the correct number of populations.
	In old file populations 15 and 33 have the same lat/long.  I missed this in previous runs.
	Do the change manually. New file is AtLatLongsp4MOD2.txt.  At /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ hclust/remakeAtdata/AtLatLongsp4MOD2.txt.

	Verify that there are no more identical lat/longs with different MINDEX:
sort -t$'\t' -n -k2,2 -k3,3 AtLatLongsp4MOD2.txt > sorted.tmp.txt;
cat sorted.tmp.txt | awk '{print $2" "$3" "$6}' | sort -t' ' -u -n -k3,3 | wc -l; # -> 41, correct (includes header)
rm sorted.tmp.txt;

	Re-sort AtSNP4aa.txt, the haplotypista input file, into the same order as AtLatLongsp4MOD.txt.
	The new haplotypista input file will be called AtSNP4aaMOD2.txt:
tail -n +2 AtLatLongsp4MOD2.txt > lltmp.txt; #remove the header
tail -n +4 AtSNP4aa.txt > snptmp.txt; #remove three header lines
sed 's/^X//g' snptmp.txt > snptmp1.txt; #remove leading X's on sample names
mv snptmp1.txt snptmp.txt;
	use fancy awk script from the internets to perform the re-sort
awk 'FNR == NR { lineno[$1] = NR; next} {print lineno[$1], $0;}' lltmp.txt snptmp.txt | sort -k 1,1n | cut -d' ' -f2- > snptmp1.txt;
	verify that re-order is correct
a=$(awk '{print $1}' AtLatLongsp4MOD2.txt | tail -n +2);
b=$(awk '{print $1}' snptmp1.txt);
diff <(echo "$a") <(echo "$b"); # <- all OK for first column
	verify that some columns with genotypes have also been reordered (snps 2-50, say):
c=$(cut -d' ' -f2-50 AtSNP4aa.txt | tail -n +4);
d=$(cut -d' ' -f2-50 snptmp1.txt);
diff <(echo "$c") <(echo "$d"); #expect differences for affected rows that should be the same as for the first column, below
a=$(cut -d' ' -f1 AtSNP4aa.txt | tail -n +4 | sed 's/^X//g');
b=$(cut -d' ' -f1 snptmp1.txt);
diff <(echo "$a") <(echo "$b"); # everything checks out OK
	verify that SNPs are retained with the same genotypes before and after sort for all samples:
samnames=$(cut -d$'\t' -f1 AtLatLongsp4MOD2.txt | tail -n +2);
for i in $samnames;
  do c=$(grep ^X"$i " AtSNP4aa.txt | sed 's/^X//g' | shasum);
    d=$(grep ^"$i " snptmp1.txt | shasum);
    echo "    $i";
    diff <(echo "$c") <(echo "$d"); #there should be no output
  done;

	now get snptmp1.txt built into a new haplotypista input file
a=$(head -3 AtSNP4aa.txt); #get header
b=$(sed 's/^/X/g' snptmp1.txt); #add back leading Xs
echo "$a" > snptmp2.txt;
echo "$b" >> snptmp2.txt;

	clean up
mv snptmp2.txt AtSNP4aaMOD2.txt;
rm lltmp.txt snptmp.txt snptmp1.txt;

	Need also to remake the *DIVAtable.txt files.  Do this by modifying the old ones to *DIVAtableMOD.txt.
	Basically, remove population 33 (same as 15) and renumber.
	
for e in env geo pco;
do a=$(sed '34d' $e"DIVAtable.txt" | cut -d$'\t' -f2-);
  b=$(seq 1 40);
  b=$(echo MINDEX; echo "$b");
  paste -d$'\t' <(echo "$b") <(echo "$a") > $e"DIVAtableMOD.txt";
done;




	On cluster perform the following set of haplotypista analyses:

./haplotypista -i AtSNP4aaMOD2.txt -o AtSNP4 -l AtSNP4log1.txt -b 1 1 -m ? -p 1 &
./haplotypista -i AtSNP4aaMOD2.txt -o AtSNP4 -l AtSNP4log2.txt -b 2 2 -m ? -p 1 &
./haplotypista -i AtSNP4aaMOD2.txt -o AtSNP4 -l AtSNP4log3.txt -b 3 3 -m ? -p 1 &
./haplotypista -i AtSNP4aaMOD2.txt -o AtSNP4 -l AtSNP4log14.txt -b 26 50 -m ? -p 1 &
./haplotypista -i AtSNP4aaMOD2.txt -o AtSNP4 -l AtSNP4log15.txt -b 51 100 -m ? -p 1 &
./haplotypista -i AtSNP4aaMOD2.txt -o AtSNP4 -l AtSNP4log16.txt -b 101 200 -m ? -p 1 &

rocks run host compute-0-1 \
'/home/reevesp/haplotypista/haplotypista -i /home/reevesp/haplotypista/AtSNP4aaMOD2.txt -o /home/reevesp/haplotypista/AtSNP4 -l /home/reevesp/haplotypista/AtSNP4log4.txt -b 4 4 -m ? -p 1 & \
/home/reevesp/haplotypista/haplotypista -i /home/reevesp/haplotypista/AtSNP4aaMOD2.txt -o /home/reevesp/haplotypista/AtSNP4 -l /home/reevesp/haplotypista/AtSNP4log5.txt -b 5 5 -m ? -p 1 & \
/home/reevesp/haplotypista/haplotypista -i /home/reevesp/haplotypista/AtSNP4aaMOD2.txt -o /home/reevesp/haplotypista/AtSNP4 -l /home/reevesp/haplotypista/AtSNP4log6.txt -b 6 6 -m ? -p 1 & '
rocks run host compute-0-2 \
'/home/reevesp/haplotypista/haplotypista -i /home/reevesp/haplotypista/AtSNP4aaMOD2.txt -o /home/reevesp/haplotypista/AtSNP4 -l /home/reevesp/haplotypista/AtSNP4log7.txt -b 7 7 -m ? -p 1 & \
/home/reevesp/haplotypista/haplotypista -i /home/reevesp/haplotypista/AtSNP4aaMOD2.txt -o /home/reevesp/haplotypista/AtSNP4 -l /home/reevesp/haplotypista/AtSNP4log8.txt -b 8 8 -m ? -p 1 & \
/home/reevesp/haplotypista/haplotypista -i /home/reevesp/haplotypista/AtSNP4aaMOD2.txt -o /home/reevesp/haplotypista/AtSNP4 -l /home/reevesp/haplotypista/AtSNP4log9.txt -b 9 9 -m ? -p 1 & '
rocks run host compute-0-4 \
'/home/reevesp/haplotypista/haplotypista -i /home/reevesp/haplotypista/AtSNP4aaMOD2.txt -o /home/reevesp/haplotypista/AtSNP4 -l /home/reevesp/haplotypista/AtSNP4log10.txt -b 10 10 -m ? -p 1 & \
/home/reevesp/haplotypista/haplotypista -i /home/reevesp/haplotypista/AtSNP4aaMOD2.txt -o /home/reevesp/haplotypista/AtSNP4 -l /home/reevesp/haplotypista/AtSNP4log11.txt -b 11 11 -m ? -p 1 & \
/home/reevesp/haplotypista/haplotypista -i /home/reevesp/haplotypista/AtSNP4aaMOD2.txt -o /home/reevesp/haplotypista/AtSNP4 -l /home/reevesp/haplotypista/AtSNP4log12.txt -b 12 12 -m ? -p 1 & '
rocks run host compute-0-5 \
'/home/reevesp/haplotypista/haplotypista -i /home/reevesp/haplotypista/AtSNP4aaMOD2.txt -o /home/reevesp/haplotypista/AtSNP4 -l /home/reevesp/haplotypista/AtSNP4log13.txt -b 13 25 -m ? -p 1 & '


	move AtSNPM.bX files to Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ hclust/haplotypista

scp reevesp@10.133.164.230:"/home/reevesp/haplotypista/AtSNP4*" /Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/Arabidopsis/Experiment2/+revision with hclust/haplotypista




	extract the NG:SS:NS values for each locus from the headers of the AtSNP4.b* files.  Make a table
	called "AtLocusWeights". These values will serve to compare with enrichment values from m+1 of the same loci.
	this is done with the AtSNP4.b* files that are the haplotypista output from the reanalysis above, located here:
	/Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/Arabidopsis/Experiment2/+revision with hclust/haplotypista
	Run on cluster.
	
**********BASH**********
f="AtLocusWeights";
g="At";
t="SNP4"; #"SNP4" or "M", depending on haplotypista output naming scheme
ploidy="1"; #"1", "2"
header=$'genomeindex\tblocklength\tblocklengthindex\tlocus\tchr\tchrindex\tNGcount\tSScount\tNScount\tNGpr\tSSpr\tNSpr\tblockstart\tblockend';
> $f.txt;
for ((b=1;b<=200;b++));
do echo "$g""$t".b$b;
  
  #diploid
  if [[ $ploidy == "2" ]]; then
    chr=$(head -1 "$g""$t".b$b | tr " " "\n" | awk 'NR%2'); #extract line 1 (chr), get odd numbered lines
    loc=$(head -3 "$g""$t".b$b | tail -1 | tr " " "\n" | awk 'NR%2'); #extract line 3 (midpt), get odd numbered lines
    st=$(head -5 "$g""$t".b$b | tail -1 | tr " " "\n" | awk 'NR%2'); #extract blockstart
    en=$(head -6 "$g""$t".b$b | tail -1 | tr " " "\n" | awk 'NR%2'); #extract blockend
  #haploid
  elif [[ $ploidy == "1" ]]; then
    chr=$(head -1 "$g""$t".b$b | tr " " "\n"); #extract line 1 (chr), get odd numbered lines
    loc=$(head -3 "$g""$t".b$b | tail -1 | tr " " "\n"); #extract line 3 (midpt), get odd numbered lines
    st=$(head -5 "$g""$t".b$b | tail -1 | tr " " "\n"); #extract blockstart
    en=$(head -6 "$g""$t".b$b | tail -1 | tr " " "\n"); #extract blockend
  fi;
  
  w=$(wc -w <<< $loc); #get the number of loci
  blI=$(seq 1 $w); #compute an index value specific to this blocklength
  chrloc=$(paste -d'_' <(echo "$chr") <(echo "$loc"));#fuse chr and locus midpoint into locus name
  
  #calculate an index by chromosome, within blocklength
chrI=$(i=1;
       prevc=0;
       for c in $chr;
         do if [ $c -eq $prevc ];
            then i=$(($i + 1));
                echo $i;
            else i=1;
                echo $i;
            fi;
         prevc=$c;
         done;
       );
               
  #diploid
  if [[ $ploidy == "2" ]]; then
    aacats=$(head -4 "$g""$t".b$b | tail -1 | tr " " "\n" | awk 'NR%2' | sed "s/:/"$'\t'"/g"); #diploid, extract line 4 (NG:SS:NS), get odd numbered lines
  #haploid
  elif [[ $ploidy == "1" ]]; then
    aacats=$(head -4 "$g""$t".b$b | tail -1 | tr " " "\n" | sed "s/:/"$'\t'"/g"); #haploid, extract line 4 (NG:SS:NS), get odd numbered lines
  fi;  
  
  aapr=$(echo "$aacats" | awk '{
sum = $1 + $2 + $3;
if (sum == 0) {
  print "na\tna\tna"; }
else {
  ng = ($1/sum);
  ss = ($2/sum);
  ns = ($3/sum);
  print ng "\t" ss "\t" ns; }
}'); #sum NGSSNS counts, and find the proportion of each category

  paste -d$'\t' <(echo "$blI") <(echo "$chrloc") <(echo "$chr") <(echo "$chrI") <(echo "$aacats") <(echo "$aapr") <(echo "$st") <(echo "$en") | sed "s:^:$b"$'\t'":g" >> $f.txt;
done;
awk 'BEGIN { OFS = "\t" } ; {print NR,$0}' $f.txt > tmp.txt; #add a sequential numeric index value (line number) to each row
#nl $f.txt | sed 's/ //g' > tmp.txt; #add a sequential numeric index value (line number) to each row
#cat tmp.txt | sed 's/\./_/g' > tmp.txt; #change . to _ in locus description so r doesn't think its a number
{ echo "$header" ; cat tmp.txt; } > $f.txt;
rm tmp.txt;
**********BASHEND**********
	On cluster, takes just a few minutes.
	
	create a file '+AtSNPtoBlockMap.txt' that shows which SNPs fall in which blocks for
	b2-b200.
**********BASHPARALLEL**********
myp() {
  b=$1;
  #define parameters
  nchr=5; #number of chromosomes
  snps=$(cat "$path""$p""snps.tmp.txt");
  
    echo "blocklength=$b";
    bnames=b"$b"; #bnames will contain the haplotype block names associated with each snp for a given blocklength
    bla=$(grep "^[0-9]*"$'\t'"$b"$'\t' "$path""$p"LocusWeights.txt); #get all lines for blocklength b

    #repeat through chromosomes
    for i in $(seq 1 $nchr);
      do echo "  chr""$i";
        snpsizes=$(echo "$snps" | grep ^"$i"_ | sed 's/'$i'_//g'); #get snp sizes for this chromosome and blocklength 1
        
        declare -a cname=( $(echo "$bla" | grep $'\t'"$i_" | awk '{print $4}') ); #locus names for this chr & blocklength
        declare -a cstart=( $(echo "$bla" | grep $'\t'"$i_" | awk '{print $13}') ); #start points for each locus in $cname
        declare -a cend=( $(echo "$bla" | grep $'\t'"$i_" | awk '{print $14}') ); #end points for each locus in $cname
        clen=$(echo ${#cname[@]}); #total length of the arrays for chromosome i
       
        startpt=0;
        for s in $snpsizes;
          do found=0;
            #test whether snp position is found within haplotype block by seeing whether it is both >= start point and <= end point
            for ((c=$startpt;c<$clen;c++)); 
              do if (( $(bc <<< "$s >= ${cstart[$c]}") )) && (( $(bc <<< "$s <= ${cend[$c]}") )); then
                  bnames="$bnames"$'\n'${cname[$c]}; #add the haplotype block name to the list that maps them to SNP name
                  found=1;
                  startpt=$c; #make search more efficient by ignoring earlier blocks in future searches
                fi;
                if [[ $found == 1 ]]; then break; fi;
              done;
            #deal with situation where SNP position is not found in a haplotype block (this will
            #sometimes happen at the end of chromosomes)
            if [[ $found == 0 ]]; then
              bnames="$bnames"$'\n'"--"; 
            fi;
          done; #s
      done; #i
      
      echo "$bnames" > "$path""$p""b$b.tmp";
}
export -f myp;

#define parameters
maxbl=200; #maximum blocklength
p="At"; export p;
path=$(pwd)"/"; export path; #path to current folder
a=$(grep "^[0-9]*"$'\t'"1"$'\t' "$path""$p"LocusWeights.txt); #get all lines with blocklength 1
snps=$(echo "$a" | awk '{print $4}'); #get a list of all SNP positions
echo "$snps" > "$path""$p""snps.tmp.txt"; #put in a file accessible by all procs so it doesn't have to be passed by parallel or loaded into the environment
outf=$(head -1 "$path""$p"LocusWeights.txt | cut -d$'\t' -f1-6);
outf="$outf"$'\n'$(echo "$a" | cut -d$'\t' -f1-6);

#begin parallel
seq 2 $maxbl | parallel --sshloginfile ~/machines --jobs 24 --env myp --env path --env p myp;

#reassemble
echo "reassembling...";
comm=""; #create a list of tmp files to paste together
for b in $(seq 2 $maxbl);
  do comm="$comm""$p""b$b.tmp ";
  done;
paste -d$'\t' <(echo "$outf") $comm > "$path""+$p""SNPtoBlockMap.txt"; #consecutively paste results from each blocklength

#clean up
rm "$path""$p"snps.tmp.txt;
for b in $(seq 2 $maxbl);
  do rm "$path""$p""b$b.tmp"; #clean up
  done;
**********BASHEND**********
	Takes ~9 hrs.


	make some test plots of the frequency of NG, NS, and SS categories in each haplotype block using
	this data (AtLocusWeights.txt) using R.  Here you use an r script called AtLocusWeightsHeatmap.r.
	This stuff is here: /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+results/+test\ heatmaps

	This doesn't work because there are many uncategorized SNPs in the Arabidopsis data set.
	These prevent R from constructing the heatmaps.  This may cause problems later...
	
		###REVISED REGION###
	remake the M+ .dat and .var input files with the new geo, env, pco binned variables, calculated 
	with R's 'hclust' function instead of 'cut' or 'cut2':

	Convert haplotypista data sets to M+:
	Run script "2. Convert haplotypista to M+.scpt"  #set so block midpoint positions end up in var file. 
	Use 1 script instance on MacMini (Mavericks).  The script takes
	as input the AtSNP4 files, which have haploid data in a single line per individual.
	Takes ~1 hr.

	
	Add Geo, Env, PCOEnv data to end of M+ files for various optimizations:
	To get new bins using R's 'hclust' function, run script "5b. Calc generalized variance for env data.scpt" with the settings 
	CenterAndScale="bin", cuttype="hclust", and binonly="yes.  Rename the output files like "envinput7.txt", "pcoinput7.txt",
	"geoinput7.txt". Then run script "7. Add target data to M+.scpt",
	which will add the new binned values to the .dat and .var files from above.  This was performed 
	with one processor on MacMini (Mavericks).
	Takes 1.5 hrs on old Script Editor.
	
###


Run m+1 analyses on ceres. This parallelizes on locus, and post processes the m+1 output,
	leaving a SUM file with stats, and a RAW file which is an archive of all the m+1 output, concatenated.
	Structure of runs on Ceres:
	b1, 8dAtb1MOD.sh numslice=2, longmem:                 sbatch --array=24 slurmAtArrayb1longmem.sh
	b2-b8, 8dMOD.sh numslice=1, medium, individually:     sbatch --array=17-23 slurmAtArrayb2-8medium.sh
	b9-b200, 8dMOD.sh numslice=1, short, in groups of 12: sbatch --array=1-16 slurmAtArrayb9-200short.sh
	
	Split up b9-200 using:
p=/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ hclust/script7outputHclust/;
dp=/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ hclust/m+deployment/hclust/;
f=9;
for j in {1..12};
do for i in {1..16};
  do cp "$p"*.b"$f".envgeopco.dat "$dp$i"; #distribute files among folders
    cp "$p"*.b"$f".RgenTenv.var "$dp$i"; #distribute files among folders
    let f++;
  done;
done;

	Split up b2-8 using:
p=/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ hclust/script7outputHclust/;
dp=/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ hclust/m+deployment/hclust/;
f=2;
for j in {1..1};
do for i in {17..23};
  do cp "$p"*.b"$f".envgeopco.dat "$dp$i"; #distribute files among folders
    cp "$p"*.b"$f".RgenTenv.var "$dp$i"; #distribute files among folders
    let f++;
  done;
done;

	Use a scheduled rsync to keep ceres from getting clogged with files:
	hclust and kmeans:
SSHPASS='qwerQWER12#$'; #enter the password as a shell variable
export SSHPASS;
while (true);
do date
  sshpass -e rsync -avz --remove-source-files --progress pat.reeves@scinet-login.bioteam.net:"/home/pat.reeves/mplusrunsAt/atb2-200/*/[R,S][A,U][W,M]*" "/Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/Arabidopsis/Experiment2/+revision with hclust/+results/hclust/";
  sshpass -e rsync -avz --remove-source-files --progress pat.reeves@scinet-login.bioteam.net:"/home/pat.reeves/mplusrunsAt/atb2-200/kmeans/*/[R,S][A,U][W,M]*" "/Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/Arabidopsis/Experiment2/+revision with hclust/+results/kmeans/";
  echo;
  echo "going to sleep for an hour starting "$(date);
  sleep 3600;
done;


	
	calculate mean M+ enrichment across core sizes, for each locus, for each blocklength. This uses the
	SUMfiles as input and produces the +STATS files.  Do this on the cluster.



**********BASHPARALLEL-FASTER-RUNINSCREEN**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all core sizes, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum "\t" sum/NF "\t" sqrt(var) "\t" NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

main() {
        l=$1;
        ofile="$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt"; #temporary outfile for each parallel process
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";

        boo=$(grep -F $l $sumfileIN | awk '{ print $13 - $10 }' | tr "\n" " "); #extract all core sizes for this locus, calculate target enrichment
        st=$(mysd "$boo"); #calculate sum, mean, sd, and n using function mysd()
        echo "$b"$'\t'"$l"$'\t'"$st" > $ofile; #write stats to output file
}
export -f main;

#parameters
t="Arabidopsis";
elist="geo pco";
bmin=1; #min blocklength
bmax=200;
nnode=8; #number of cluster nodes
blist=$(seq $bmin $bmax); #blocklengthrange
p=$(pwd)"/"; export p;

#parallelize on locus, within loop on blocklength
for e in $elist;
  do echo $e;
    export e;
    for b in $blist;
      do echo "b=$b";
        export b;
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";
        uloc=$(cut -d$'\t' -f2 "$sumfileIN" | tail -n +2 | uniq); #find unique loci from all loci in column 2
        numloc=$(echo "$uloc" | wc -l);
        
        #parallel step
        #sends $gnuN records to each node via parallel step #1, this is piped to parallel step #2, which starts 96 jobs per node from the instruction set it has received
        echo "processing $sumfileIN...";
        gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node
        #echo "$uloc" | parallel --sshloginfile ~/machines --jobs 24 --env main --env mysd --env e --env b --env p main;
        echo "$uloc" | parallel --sshloginfile ~/machines --jobs 1 --env main --env mysd --env e --env b --env p --pipe -N"$gnuN" /home/reevesp/bin/parallel -j96 --env main --env mysd --env e --env b --env p main;

        #consolidate results for each locus into one file
        ofile="+STATS.R"$e"Tgen.b"$b".txt"; #temporary outfile
        > "$p"$ofile;
        for l in $uloc;
          do echo "concatenating loc file $l to +STATS.R"$e"Tgen.b"$b".txt ...";
            cat "$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt" >> "$p"$ofile;
          done;
          #test whether all were concatenated
          #s1=$(wc -l "$p""+STATS.R"$e"Tgen.loc"*".b"$b".txt" | grep total | awk '{print $1}');
          s1=$(ls "$p" | grep -c "loc"); #count the number of files with "loc" in their name using grep. can't wc -l because too many arguments
          s2=$(wc -l "$p"$ofile |  awk '{print $1}');
         if [ $s1 == $s2 ]; 
          then find "$p" -maxdepth 1 -name "*loc*" -delete; #remove files when there are too many for rm
          #rm "$p""+STATS.R"$e"Tgen.loc"*".b"$b".txt"; #remove files if number of lines in concat file is the same as the sum of all input files
          else echo "s1="$s1", s2="$s2". Consolidated file "$ofile" no good. Quitting..." > err.txt;
            exit 1;
        fi;
      done; #$b
      
    #concatenate files for this Rgeo/Rpco
    statsfile="$p""+STATS."$t".R"$e"Tgen.txt";
    echo blocklength$'\t'locus$'\t'sum$'\t'mean$'\t'sd$'\t'n > $statsfile;
    >"$statsfile"TMP;
    for b in $blist;
      do echo "concatenating +STATS.R"$e"Tgen.b"$b".txt ...";
        cat "$p""+STATS.R"$e"Tgen.b"$b".txt" >> "$statsfile"TMP;
      done;
    cat "$statsfile"TMP >> "$statsfile";

    #test whether all individual +STATS files have been added to the final concatenated file using number of lines
    s1=$(wc -l +STATS.R"$e"Tgen.b*.txt | grep total | awk '{print $1}');
    s2=$(wc -l "$statsfile" | awk '{print $1}');
    s2=$(( $s2 - 1 ));
    if [ $s1 == $s2 ]; 
      then rm +STATS.R"$e"Tgen.b*.txt; #remove files if number of lines in concat file is the same as the sum of all input files
        rm "$statsfile"TMP;
      else echo "the number of lines ain't the same. something is 'crewed.";
    fi;

  done; #$e
**********BASHEND**********
	Takes ~7hrs

	add the proportion of NG,SS,NS sites from AtLocusWeights.txt to the +STATS summary file,
	generating the +v2STATS file.  Requires the AtLocusWeights.txt file as input.  

**********BASH**********
myAddaa() {
  weightsfile=$1;
  e=$2;
  key=$3;
  
  statsfile="+STATS.""$key"".R"$e"Tgen.txt";
  s=$(cut -d$'\t' -f2 $statsfile | tail -n +2 | md5sum); #get md5sum of locus ID column from the +STATS file
  w1=$(cut -d$'\t' -f4 $weightsfile | sed 's/_/./g'); #get the locus ID less blocklength index column, replace _ with . for md5 check
  w2=$(cut -d$'\t' -f3 $weightsfile); #get the blocklength index column
  w=$(paste -d'.' <(echo "$w1") <(echo "$w2") | tail -n +2 | md5sum); #get md5sum of locus ID column from the Weights file
  if [ "$s" == "$w" ]; #verify that files are in same order using a checksum of the unique locus ID column
    then boo=$(sed 's/blocklength/blocklength2/g' $weightsfile | sed 's/locus/locus2/g'); #modify the header in the Weights file so that there are no redundant column names after combining
      paste -d$'\t' $statsfile <(echo "$boo") > "+v2STATS.""$key"".R"$e"Tgen.txt"; #paste the statsfile and the weights file together
    else echo "md5sums do not match, $e, aborting.";
      kill -INT $$; #terminate the script, return to the shell
  fi;
}
export -f myAddaa;

weightsfile="AtLocusWeights.txt";
elist="geo pco";
key="Arabidopsis";
parallel --env myAddaa myAddaa ::: "$weightsfile" ::: $elist ::: "$key";

**********BASHEND**********

#Plot summed enrichment by genomic position across blocklengths as a heatmap. Use Rscript
#AtGenomicGeography.r. Uses the +v2STATS files as input.


	Using the file +*SNPtoBlockMap.txt, calculate the sum of M+ enrichment value across blocklengths,
	for each SNP, i.e. for each position in the genome.  These data will be extracted from the 
	+v2STATS* files in several steps.  The first makes a file containing a table with the max 
	enrichment value at each site for each block, *EnrichAcrossBlocks*. Step 1:

**********BASHPARALLEL**********
myp() {
  e=$1;
  b=$2;
  v2name="$p""+v2STATS."$v2p".R"$e"Tgen.txt";
        
        echo "e=$e  b=$b";
        bout="$p""rr.R"$e"Tgen.b""$b"".tmp"; #temporary outfile

        #cut out the locus id and sum columns for current blocklenth
        gg=$(grep ^$b$'\t' "$v2name" | cut -d$'\t' -f1-3);
        
        #cut column with locus names for current blocklength from SNPtoBlockMap
        if [ $b = 1 ]; then
          cc=$(cut -d$'\t' -f4 "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        else col=$(( $b + 5 )); #calculate the column to cut out for each blocklength
          cc=$(cut -d$'\t' -f$col "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        fi;

        #search for locus names in $gg, retrieve enrichment value
        rr="b$b";

#put another parallel here if necessary
        for l in $cc;
          do if [[ $l = "--." ]]; then
                r="--"; #test for empty enrichment value
              else
                r=$(grep -m 1 ^"$b"$'\t'"$l" <<<"$gg" | awk '{print $3}');
              fi;
            rr+=$'\n'$r; #add current locus enrichment value to list
          done;

         #write out the result
         echo "$rr" > "$bout";
         #truncate -s -1 "$bout"; #remove the hanging newline, use 'gtruncate' only on osx
}
export -f myp;

#parameters
v2p="Arabidopsis"; export v2p;
t="At";
minbl=1;
maxbl=50; # 50, 200
blist=$(seq $minbl $maxbl);
elist="geo pco";

p=$(pwd)"/"; export p;
tname="$p""+""$t""SNPtoBlockMap.txt"; export tname;

#begin parallel
parallel --sshloginfile ~/machines --jobs 24 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist;
  
#reconstitute
echo "reconstituting..."
for e in $elist;
  do echo "e=$e";
    comm=""; #create a list of tmp files to paste together
    for b in $blist;
      do comm+="$p""rr.R"$e"Tgen.b$b.tmp ";
    done;
    pos=$(cut -d$'\t' -f4 $tname | sed 's/'^.*_'//g' | sed 's/locus/pos/g');
    h=$(cut -d$'\t' -f1-6 $tname | tr $'\t' ' ');
    h2=$(paste -d' ' <(echo "$h") <(echo "$pos"));
    paste -d' ' <(echo "$h2") $comm > "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt";
    rm $comm; #clean up
  done;#e
  
**********BASHEND**********
	Takes ~8 hours.


	Process each line (SNP position) of *EnrichAcrossBlocks*.txt to get sum across all blocklengths
	of cumulative enrichment across all core sizes of max enrichment value from 10 reps of M+.
	Produces a file called SUM*EnrichAcrossBlocks... Step 2:

**********BASH**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all blocklengths, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       #returns space delimited string like "sum mean sd n"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum " " sum/NF " " sqrt(var) " " NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

myp() {
    l=$1

    p=$(echo "$l" | cut -d' ' -f8-); #cut out the portion of the line with max enrichment values

    #skip loci that have missing blocks
    if [[ $p != *"--"* ]]; then
      m=$(mysd "$p"); #calculate sum, mean, sd, n
      chr=$(echo "$l" | cut -d' ' -f5);
      pos=$(echo "$l" | cut -d' ' -f7);
      i=$(echo "$l" | cut -d' ' -f1); #get genome index
      #loc=$chr"."$pos"."$(echo "$l" | cut -d' ' -f6);
      loc=$chr"."$pos"."$i;
      echo "snp=$i";

      echo $i" "$chr" "$pos" "$loc" "$m > "$path"$i".tmp"; #add new values to growing string, write to tmp file
    fi;
}
export -f myp;

#parameters
elist="geo pco";
t="At";
path=$(pwd)"/"; export path;

#parallel
for e in $elist;
do echo $e;
  nnode=8;
  numloc=$(wc -l +"$t"EnrichAcrossBlocks.R"$e"Tgen.txt | awk '{print $1}');
  numloc=$(( $numloc - 1 ));
  gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node

  #tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --sshloginfile ~/machines --jobs 24 --env myp --env mysd --env path myp;
  tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --sshloginfile ~/machines --jobs 1 --env myp --env mysd --env path --pipe --round-robin /home/reevesp/bin/parallel -j96  --env myp --env mysd --env path myp;
  
  #assemble results
  echo "concatenating...";
  ofile="+"$t"SUMEnrichAcrossBlocks.R"$e"Tgen.txt";
  o="genomeindex chr pos loc sum mean sd n";
  echo "$o" > $ofile;
  
  find "$path" -name "*.tmp" -print0 | xargs -0 cat | sort -n -t' ' -k1 >> $ofile;
    
  #clean up
  echo "clean up...";
  find "$path" -name "*.tmp" -print0 | xargs -0 rm -f
done;
**********ENDBASH**********
	Takes ~10 minutes

	Use R to make some pretty plots of significant enrichment across the genome. Modify setwd for Hclust and Kmeans.

**********BEGINR**********
#install.packages("fitdistrplus")
#install.packages("plotrix")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ hclust/+results/hclust/+AtHclustSummary\ files/ToMaxBL50")
library(fitdistrplus)
library(plotrix)

t<-"At"
usezero="no" #this sets whether to include significant sites only when they are also > or < 0
             #to use the zero criterion, say "yes", to just use the ecdf around the mean, say "no"
pvals<-c(0.001,0.05)
for (pp in pvals)
{
	pdf(file=paste(t,"M+GenomeScan",pp,".pdf", sep=""))
	par(mfrow=c(2,1))

		i="geo"
		g <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",i,"Tgen.txt", sep=""), header=TRUE, sep=" ")
		mycdf <- ecdf(g$sum) #calculate the cdf of the empirical distribution
		p <- mycdf(g$sum) #get the probability of each observation (or lower)
		g$pecdf <- p #add new column to g, containing the probability from the empirical cdf
		g$index <- 1:nrow(g) #add a column with the sequential index
		h <- g[which(g$pecdf > 1-pp),] #get snps that are highly enriched
		if ( usezero == "yes" ) {
			hgeo <- h[which(h$sum>0),] #remove snps where sum<0, i.e. less than random
		}
		else {
			hgeo <- h # do not remove snps based on positive or negative
		}
		l <- g[which(g$pecdf < pp),] #get snps that are lowly enriched
		if ( usezero == "yes" ) {
			lgeo <- l[which(l$sum<0),] #remove snps where sum>0, i.e. more than random
		}
		else {
			lgeo <- l # do not remove snps based on positive or negative
		}
		highposgeo<-hgeo$index #get the sequential position of the highly enriched
		lowposgeo<-lgeo$index #get the sequential position of the deficient

		plot(g$index, g$sum, type="l", main=paste(t,i,pp,sep=""))
		color.scale.lines(g$index, g$sum, col=factor(g$chr))
		points(g$index[highposgeo], g$sum[highposgeo], col = "red", cex=0.7) #mark significant values
		points(g$index[lowposgeo], g$sum[lowposgeo], col = "blue", cex=0.7) #mark significant values

		i="pco"
		g <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",i,"Tgen.txt", sep=""), header=TRUE, sep=" ")
		mycdf <- ecdf(g$sum) #calculate the cdf of the empirical distribution
		p <- mycdf(g$sum) #get the probability of each observation (or lower)
		g$pecdf <- p #add new column to g, containing the probability from the empirical cdf
		g$index <- 1:nrow(g) #add a column with the sequential index
		h <- g[which(g$pecdf > 1-pp),] #get snps that are highly enriched
		if ( usezero == "yes" ) {
		hpco <- h[which(h$sum>0),] #remove snps where sum<0, i.e. less than random
		}
		else {
			hpco <- h
		}
		l <- g[which(g$pecdf < pp),] #get snps that are lowly enriched
		if ( usezero == "yes" ) {
		lpco <- l[which(l$sum<0),] #remove snps where sum>0, i.e. more than random
		}
		else {
			lpco <- l
		}
		highpospco<-hpco$index #get the sequential position of the highly enriched
		lowpospco<-lpco$index #get the sequential position of the deficient

		plot(g$index, g$sum, type="l", main=paste(t,i,pp,sep=""))
		color.scale.lines(g$index, g$sum, col=factor(g$chr))
		points(g$index[highpospco], g$sum[highpospco], col = "red", cex=0.7) #mark significant values
		points(g$index[lowpospco], g$sum[lowpospco], col = "blue", cex=0.7) #mark significant values

		#calculate intersection of well-collected significant snps
		isect <- hgeo[is.element(hgeo$pos, intersect(hgeo$pos,hpco$pos)),]
		write.table(isect, file=paste(t,"HighCommonToGeoPco",pp,".txt", sep=""), row.names=FALSE, quote=FALSE)

		#calculate intersection of poorly-collected significant snps
		isect <- lgeo[is.element(lgeo$pos, intersect(lgeo$pos,lpco$pos)),]
		write.table(isect, file=paste(t,"LowCommonToGeoPco",pp,".txt", sep=""), row.names=FALSE, quote=FALSE)
	
		#create BED formatted output for hgeo, hpco. verify first whether there exist any snps that are
		#both significantly enriched and better than random (or significantly suppressed and worse than random.)
		#if not, do not print an output table.
		if ( length(hgeo$pos) != 0 ) {
		hgeobed <- data.frame(chr=paste("Chr", hgeo$chr, sep=""), chromStart=hgeo$pos, chromEnd=hgeo$pos+1)	
		write.table(hgeobed, file=paste(t,"hGeo",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}
		if ( length(hpco$pos) != 0 ) {
		hpcobed <- data.frame(chr=paste("Chr", hpco$chr, sep=""), chromStart=hpco$pos, chromEnd=hpco$pos+1)	
		write.table(hpcobed, file=paste(t,"hPco",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}
		
		#create BED formatted output for lgeo, lpco
		if ( length(lgeo$pos) != 0 ) {
		lgeobed <- data.frame(chr=paste("Chr", lgeo$chr, sep=""), chromStart=lgeo$pos, chromEnd=lgeo$pos+1)	
		write.table(lgeobed, file=paste(t,"lGeo",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}
		if ( length(lpco$pos) != 0 ) {
		lpcobed <- data.frame(chr=paste("Chr", lpco$chr, sep=""), chromStart=lpco$pos, chromEnd=lpco$pos+1)	
		write.table(lpcobed, file=paste(t,"lPco",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}

	dev.off()
}
**********ENDR**********

	Determine the identity of each important SNP

	Download TAIR 10 annotation
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff 
	Convert to BED format using bedops, extract only genes from GFF file
cat TAIR10_GFF3_genes.gff | gff2bed | grep -w gene > TAIR10genes.bed

	Create background data file, a bed file containing only those genes for which there are SNPs.
	Create files listing genes enriched by M+ analyses for downstream GO over-representation analysis.
	Uses bedops. For whatever reason, this script has to be run in two parts or it fails.
**********BASH**********
***Part 1***
i="At"; #taxon identifier
ag="TAIR10genes.bed"; #bed file containing all annotated genes for given genome release

#Create a bed file for each SNP locus where we have a measurement (only needs to be done for geo or pco, since SNPs are the same):
chr=$(cut -d' ' -f2 "+"$i"SUMEnrichAcrossBlocks.RgeoTgen.txt" | tail -n +2 | sed 's/^/Chr/g');
pos=$(cut -d' ' -f3 "+"$i"SUMEnrichAcrossBlocks.RgeoTgen.txt" | tail -n +2);
pos2=$(awk '{$1 = $1 + 1; print}' <<<"$pos"); 
paste -d$'\t' <(echo "$chr") <(echo "$pos") <(echo "$pos2") > $i"AllSNPs.bed";

#Sort to make sure it will work for bedops --element-of
/usr/local/bin/sort-bed $i"AllSNPs.bed" > tmp.txt;
mv tmp.txt $i"AllSNPs.bed";

#Calculate the intersection of AtAllSNPs.bed and TAIR10genes.bed, the list of genes for which there are SNPs in the study (the background):
/usr/local/bin/bedops --element-of 1 "$ag" $i"AllSNPs.bed" > $i"AllGenesWithSNPs.bed";
  
#Make background bed files into a simple gene list.
cut -d$'\t' -f4 $i"AllGenesWithSNPs.bed" | sort -u > $i"BackgroundGenesList.txt";
#clean up
rm $i"AllSNPs.bed";
rm $i"AllGenesWithSNPs.bed";

***Part 2***
#Compute lists of genes enriched by M+
for k in "Geo" "Pco";
  do echo $k;
    #Verify that bed files are sorted correctly using bedops:sort-bed, e.g.
    for j in "0.001" "0.05";
      do echo "  $j";
        for h in "h" "l";
          do echo "    $h";
            /usr/local/bin/sort-bed $i$h$k$j.bed > tmp.txt;
            mv tmp.txt $i$h$k$j.bed;
    
            #Determine nearest gene to highly enriched snps using closest-features
            /usr/local/bin/closest-features --closest --dist $i$h$k$j.bed "$ag" > $i$h"$k"EnrichedGenes$j.bed;
    
            #Extract hits that are within genes, not just close
            grep "|0$" $i$h"$k"EnrichedGenes$j.bed > $i$h"$k"EnrichedWithinGenes$j.bed;
            rm $i$h"$k"EnrichedGenes$j.bed;
    
            #Get gene identifiers for downstream GO analysis
            cut -d$'\t' -f6 $i$h"$k"EnrichedWithinGenes$j.bed | sort -u > $i$h"$k"EnrichedGenesList$j.txt;
            rm $i$h"$k"EnrichedWithinGenes$j.bed;
          done;
      done;
  done;
**********ENDBASH**********

	Use Python GOATOOLS to analyze over-representation:
easy_install goatools;
easy_install fisher;
easy_install statsmodels;
wget http://geneontology.org/ontology/go-basic.obo;
wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo;

	Download a raw gene association file (GAF) from http://amigo.geneontology.org/amigo/search/annotation
	Choose species in 'Species' menu, remove 'not', 'contributes_to', and "colocalizes_with'
	in 'Annotation Qualifier' menu.  For Arabidopsis, this was done in three downloads, one for
	Cellular Component, one for Molecular Function, and one for Biological Process due to limits
	on the number of lines that could be downloaded at once. Use cat * and call it ATRawGAF.txt. This GAF will include GO terms for 
	Biological Process, Molecular Function, and Cellular Component ontologies.
	Process the file into an "association" input file, called AtGAF.txt, for GOATOOLS script find_enrichment.py


**********BASHPARALLEL**********
myp() {
    line=$1;
    o=$(mktemp);
    f="";
    f=$(grep ^"$line"$'\t' in.txt | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" > $o;
    fi;
}
export -f myp;

t=At;
cut -d$'\t' -f3 "$t"RawGAF.txt | sort -u > GAFgenestmp.txt;
cut -d$'\t' -f3,5 "$t"RawGAF.txt > in.txt;
rm /tmp/tmp.*;
cat GAFgenestmp.txt | parallel --env myp myp;

cat /tmp/tmp.* | sort | sed 's/;$//g' > ./"$t"GAF.txt;
rm GAFgenestmp.txt;
rm /tmp/tmp.*;
**********ENDBASH**********
	Takes a minute or so on glitch.
	
	One gene, FIP1[V], is missing because grep can't handle the brackets.  So, get it manually:
f=$(grep ^'FIP1\[V\]'$'\t' in.txt | cut -d$'\t' -f2 | tr '\n' ';' | sed 's/;$//g');
cat AtGAF.txt <(echo 'FIP1[V]'$'\t'"$f") | sort > tmp.txt;
mv tmp.txt AtGAF.txt;

	Run GOATOOLS find_enrichment.py to get the enrichment status, over (e), or under (p),
	representation (column 3 in output).
	The --no_propagate_counts selects only the least inclusive GO term, i.e no parent terms.
	This eliminates multiple significant results along a parent-child path, but has the
	undesirable consequence of making significant under-representation (p) a questionable
	result because higher level terms are "shorted" their descendant terms.

	Move the relevant items (inc. *GAF.txt, go-basic.obo) to a nested folder, "statistical tests", then run
	the GOATOOLS script find_enrichment.py:
**********BASH**********
mv *EnrichedGenesList* "statistical tests";
mv AtBackgroundGenesList.txt "statistical tests";
cd "statistical tests";

t=At;
comm="--alpha 0.05 --pval 0.05 --obo go-basic.obo --no_propagate_counts --method holm";

for p in 0.05 0.001;
  do echo $p;
  for e in Geo Pco;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          find_enrichment.py $(echo "$comm" "$t""$h""$e"EnrichedGenesList"$p".txt "$t"BackgroundGenesList.txt "$t"GAF.txt) > +"$t""$h""$e""$p"GOAout.txt;
        done;
    done;
  done;
**********ENDBASH**********

	Examine first part of output files to determine if any GO terms are significantly over-
	or under-represented.
head -20 +*GOA*;

	FINAL RESULT:

For MaxBL=50:
	There are four significantly over-represented GO terms for Arabidopsis among genes that are
	well-collected by environmental data and Hclust.  Three of these GO terms are close but no cigar over-represented
	among genes that are well-collected by geographic data. The three significant results for +AthPco0.05GOAout.txt
	are attributable to the same cluster of genes: AT4G03440 AT4G03450 AT4G03460 AT4G03470 AT4G03480 AT4G03490 AT4G03500.
	These can be found in /Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/Arabidopsis/Experiment2/+revision with hclust/+results/hclust/+AtHclustGO work/ToMaxBL50/gene identification/README.txt.
	Similar findings occurred for MaxBL=200 Hclust and Rcut.
	Using Hclust:
	==> +AthPco0.001GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0009117	BP	e	nucleotide metabolic process	2/33	4/24259	1.07e-05	n.a.	2	0.0183
	==> +AthPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:2000031	BP	e	regulation of salicylic acid mediated signaling pathway	7/1424	10/24259	2.43e-07	n.a.	7	0.000413
	GO:0071446	BP	e	cellular response to salicylic acid stimulus	7/1424	11/24259	6.34e-07	n.a.	7	0.00108
	GO:0031347	BP	e	regulation of defense response	7/1424	17/24259	2.73e-05	n.a.	7	0.0464
	
For MaxBL=200:
	Using Rcut:
	There are three significantly over-represented GO terms among Arabidopsis genes that are well-collected
	using geographic data. There are no significantly over- or under-represented GO terms among 
	Arabidopsis genes that are well- or poorly- collected using environmental data and Rcut.
	==> +AthGeo0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:2000031	BP	e	regulation of salicylic acid mediated signaling pathway	7/896	10/24203	1.01e-08	n.a.	7	1.73e-05
	GO:0071446	BP	e	cellular response to salicylic acid stimulus	7/896	11/24203	2.7e-08	n.a.	7	4.59e-05
	GO:0031347	BP	e	regulation of defense response	7/896	17/24203	1.31e-06	n.a.	7	0.00222

	Using Rcut2:
	There is one significantly under-represented GO term among Arabidopsis genes that are well-collected
	using environmental data.  There are two significantly over-represented GO terms among Arabidopsis genes
	that are poorly collected using environmental data and Rcut2.
	==> +AthPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0008150	BP	p	biological_process	419/2232	5740/24203	4.67e-09	n.a.	419	7.95e-06
	==> +AtlPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0006952	BP	e	defense response	27/848	235/24203	6.76e-08	n.a.	27	0.000115
	GO:0008150	BP	e	biological_process	254/848	5740/24203	2.23e-05	n.a.	254	0.0379

	Using Hclust:
	There are three significantly over-represented GO terms for Arabidopsis among genes that are
	well-collected by geographic data and Hclust.  The same three GO terms are significantly over-represented
	among genes that are well-collected by environmental data. These same findings occurred for Rcut.
	There is one significantly over-represented GO term among genes that are poorly collected
	by geographic data.
	==> +AthGeo0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:2000031	BP	e	regulation of salicylic acid mediated signaling pathway	7/1236	10/24203	9.33e-08	n.a.	7	0.000159
	GO:0071446	BP	e	cellular response to salicylic acid stimulus	7/1236	11/24203	2.45e-07	n.a.	7	0.000417
	GO:0031347	BP	e	regulation of defense response	7/1236	17/24203	1.1e-05	n.a.	7	0.0187
	==> +AthPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:2000031	BP	e	regulation of salicylic acid mediated signaling pathway	7/1304	10/24203	1.35e-07	n.a.	7	0.00023
	GO:0071446	BP	e	cellular response to salicylic acid stimulus	7/1304	11/24203	3.53e-07	n.a.	7	0.000601
	GO:0031347	BP	e	regulation of defense response	7/1304	17/24203	1.56e-05	n.a.	7	0.0266
	==> +AtlGeo0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0000272	BP	e	polysaccharide catabolic process	5/1576	6/24203	6.6e-06	n.a.	5	0.0112

	Using Kmeans:
	There are two significantly over-represented GO terms for Arabidopsis among genes in regions that are
	well-collected by environmental information:
	==> +AthPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:2000031	BP	e	regulation of salicylic acid mediated signaling pathway	6/1652	10/24203	1.66e-05	n.a.	6	0.0282
	GO:0009117	BP	e	nucleotide metabolic process	4/1652	4/24203	2.16e-05	n.a.	4	0.0368
	There are three significantly over-represented GO terms for Arabidopsis among genes in regions that are
	poorly-collected using geographic information:
	==> +AtlGeo0.001GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0000272	BP	e	polysaccharide catabolic process	4/44	6/24203	1.42e-10	n.a.	4	2.42e-07
	GO:0006032	BP	e	chitin catabolic process	4/44	9/24203	1.19e-09	n.a.	4	2.02e-06
	GO:0016998	BP	e	cell wall macromolecule catabolic process	4/44	13/24203	6.71e-09	n.a.	4	1.14e-05

	Above result with p<0.001 consistent with below result with p<0.05:
	==> +AtlGeo0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0000272	BP	e	polysaccharide catabolic process	5/1859	6/24203	1.49e-05	n.a.	5	0.0254




	Make a pie chart showing GOslim category representation for genes enriched by M+ at the
	0.05 level.

	Extract the three GO categories (biological_process, cellular_component, molecular_function)
	from goslim_plant.obo.  Use PERL modules go-perl > go-filter-subset.pl.

for i in biological_process cellular_component molecular_function;
  do
    /Users/shrub/perl5/bin/go-filter-subset.pl -namespace "$i" goslim_plant.obo > goslim_plant_"$i".obo;
  done;

	For whatever reason (actually, it is because go-filter-subset.pl includes 'part-of's, not
	just 'is_a's as members of a category), two "cellular_component"s remain in goslim_plant_biological_process.obo.
	Ten "biological_process"s remain in goslim_plant_molecular_function.obo.
	Remove them manually.
	
	Make a gene association file (GAF) for M+ well-collected and poorly-collected genes.
	Input files: At[h,l][Geo,Pco]EnrichedGenesList0.05.txt, AtRawGAF.txt

**********BASHPARALLEL**********
myp() {
    line=$1;
    o=$(mktemp /tmp/tmp.XXXXXXXX);
    f="";
    f=$(grep ^"$line"$'\t' in.txt | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" > $o;
    fi;
}
export -f myp;

t=At;
rm /tmp/tmp.*; #clear tmp directory of tmp. files generated by this script
cut -d$'\t' -f3,5 "$t"RawGAF.txt > in.txt; #make a simplified GAF file of all genes to query
  for e in Geo Pco;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          cat $t$h$e"EnrichedGenesList0.05.txt" | parallel --env myp myp;

          cat /tmp/tmp.* | sort | sed 's/;$//g' > ./$t$h$e"0.05GAF.txt";
          rm /tmp/tmp.*;
        done;
    done;
rm in.txt;
**********ENDBASH**********

	Count the number of occurrences of each GO slim term in the gene association files
	derived from the EnrichedGenesList(s). These are the files At[h,l][Pco,Geo]0.05GAF.txt.
	Do this for each major GO category (bp, cc, mf). Output files are like At[h,l][Pco,Geo]piechart_[bp,cc,mf].txt.
	Goal here is to calculate the frequency of the function in the well/poorly collected genes and
	compare that to the frequency of the same functions when all annotated genes are considered.
**********BASH**********
  
t="At";
b=$(wc -l "$t"GAF.txt | awk '{print $1}'); #the total number of annotated genes
for j in biological_process cellular_component molecular_function;
  do echo $j;
    map_to_slim.py --slim_out=direct --association_file=$t"GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t"slim_"$j".txt; #get the goslim_plant terms associated with all annotated genes using goatools map_to_slim.py
    a=$(awk -F$'\t' '{print $2}' "$t"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for all collected genes.

    for e in Geo Pco;
      do echo "  $e";
        for h in "h" "l";
          do echo "    $h";
            map_to_slim.py --slim_out=direct --association_file=$t$h$e"0.05GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t$h$e"slim_"$j".txt; #get the goslim_plant terms associated with the genes using goatools map_to_slim.py

            go=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | tr ";" "\n" | sort -u); #list of unique GO terms for well/poorly collected genes
            o=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for well/poorly collected genes.
            n=$(wc -l "$t$h$e"slim_"$j".txt | awk '{print $1}'); #number of lines = number of genes for well/poorly collected genes
            n=$(( $n - 4 )); #subtract off header lines
            echo "GOterm ObsTotNumGenes ObsTotNumFunctionHits ObsGOcount obsfreq TotNumAnnotatedGenes ExpTotNumFuncHits ExpGOcount expfreq ratio absratio plratio diff absdiff" > Piechart"$t$h$e"_"$j".txt; #for each unique GO term, count how many times it occurs in slim.txt, divide by number of genes to get the proportion of genes representing that function in the M+ enriched/purified gene set.

            for i in $go;
              do f=$(grep $i "$t$h$e"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for well/poorly collected genes
                x=$(grep $i "$t"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for all genes
                g=$(echo "scale=4;$f/$o" | bc); #freq of function among all enriched functions
                y=$(echo "scale=4;$x/$a" | bc); #freq of function among all functions
                z=$(echo "$g - $y" | bc); #diff btw freq of function among all functions, and among enriched functions (negative means less common in enriched subset)
                zz=${z#-}; #absolute value of z
                rrx=$(echo "scale=4;$g/$y" | bc); #for log2 plotting
                #use below 4 lines for non-log2 axes
                if (( $(bc <<< "$g >= $y") )); 
                  then rr=$(echo "scale=4;$g/$y" | bc); #over represented results in positive fold enrichment
                  else rr=$(echo "scale=4;-$y/$g" | bc); #under represented results in negative fold enrichment
                fi; #ratio of frequency in enriched regions / frequency in all regions
                rrr=${rr#-}; #absolute value of rr
                 echo "$i $n $o $f $g $b $a $x $y $rr $rrr $rrx $z $zz" >> Piechart"$t$h$e"_"$j".txt;
              done;
              
            rm "$t$h$e"slim_"$j".txt;
          done;
      done;
  done;
**********ENDBASH**********


	Get verbal description of slimmed GO terms for eventual plotting. Exclude the cellular_component
	as it is just the location the protein is found:
t="At";
for e in Geo Pco;
  do echo "  $e";
    for h in "h" "l";
      do echo "    $h";
        > o.txt; #freq output
        > p.txt; #ratio output
        for j in biological_process molecular_function;
          do echo $j;
            #filter out singletons, GO terms with only 1 observation in the enriched regions
            sin=$(grep -v "GO:"[0-9]*" "[0-9]*" "[0-9]*" "1" " Piechart$t$h$e"_"$j.txt | tail -n +2);

            #extract columns holding freq and ratio metrics for reformatting
            f=$(echo "$sin" | cut -d' ' -f1); #get GO slim terms
            p=$(echo "$sin" | cut -d' ' -f13); #get diff in freq of function btw enriched functions and all functions
            pp=$(echo "$sin" | cut -d' ' -f14); #get abs of diff in freq of function btw enriched functions and all functions
            q=$(echo "$sin" | cut -d' ' -f10); #get ratio between enriched/all
            qq=$(echo "$sin" | cut -d' ' -f11); #get abs of ratio between enriched/all
            pl=$(echo "$sin" | cut -d' ' -f12); #get ratio between enriched/all to be plotted on log2 axis
            c=$(for i in $f;
              do sed -n -e '/id: '$i'/,$p' goslim_plant_"$j".obo | head -2 | tail -1 | cut -d' ' -f2-;
              done;);
            d=$(for i in $f; do echo "$j"; done;); #create a header column with bp, cc, mf category
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$p")  <(echo "$pp") >> o.txt; #write to temporary freq output file
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$q")  <(echo "$qq") <(echo "$pl") >> p.txt; #write to temporary ratio output file
         done;
          echo Gocat$'\t'Goterm$'\t'name$'\t'diff$'\t'absdiff > PlotRpieFREQ$t$h$e.txt; #make header for freq output file
          echo Gocat$'\t'Goterm$'\t'name$'\t'ratio$'\t'absratio$'\t'plratio > PlotRpieRATIO$t$h$e.txt; #make header for ratio output file
          sort -t$'\t' -nr -k5,5 o.txt >> PlotRpieFREQ$t$h$e.txt;
          sort -t$'\t' -nr -k5,5 p.txt >> PlotRpieRATIO$t$h$e.txt;
         rm o.txt p.txt;
      done;
  done;
	
	Use R to construct charts, change hclust and kmeans as necessary in setwd:
**********RSCRIPT**********
#install.packages("ggplot2")
#library(ggplot2)
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ hclust/+results/hclust/+AtHclustGO\ work/ToMaxBL50/pie\ charts")

t<-"At"
rvals<-c("RATIO", "FREQ")
cvals<-c("h", "l")
evals<-c("Geo","Pco")
for (rr in rvals)
{
  for (cc in cvals)
  {
    for (ee in evals)
    {
      pdf(file=paste("GOchart",rr,t,cc,ee,".pdf", sep=""))
      #par(mfrow=c(1,2))
  
      g <- read.table(paste("PlotRpie",rr,t,cc,ee,".txt", sep=""), header=TRUE, sep="\t")
      h <- g[g$Goterm!="GO:0008150",] #Remove overarching GO categories biological_process, molecular_function
      h <- h[h$Goterm!="GO:0003674",]
  
      par(mar = c(4,20,4,2) + 0.1)
      if(rr == "FREQ") {
        barplot(h$diff, main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$diff),max(h$diff)))
        } 
      else {
        #barplot(h$plratio, log="x", main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$plratio),max(h$plratio)))
        hh=h #transfer h table to a new variable, which will be modified
        hh$id=paste(hh$Goterm,hh$name,sep=" ") #assemble label term
        hh$id <- factor(hh$id, levels=hh$id) #fix order by making id a factor with levels

        ggout<-ggplot(hh,aes(id,plratio,width=0.8)) + 
        geom_bar(stat="identity",color="black",fill="dark grey",size=0.25) + 
        scale_y_continuous(trans='log2', breaks=c(0.125,0.25,0.5,1,2,4,8)) + 
        coord_flip() + 
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text=element_text(size=12, family="ArialMT"))
        
        print(ggout)
        }      
      dev.off()
    }
  }
}
**********ENDR**********

	For Kmeans, because sometimes well-collected regions do not exist for pco because they are all worse than random.
**********RSCRIPT**********
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ Rcut/+AtRcutGO\ work/pie\ charts")

t<-"At"
cvals<-c("h", "l")
evals<-c("Pco","Geo")
for (cc in cvals)
{
  for (ee in evals)
  {
    pdf(file=paste("GOchart",t,cc,ee,".pdf", sep=""))
    #par(mfrow=c(1,2))

    if ( ( cc=="h" ) & ( ee=="Pco" ) ) {
      #do nothing
    }
    else {
      g <- read.table(paste("PlotRpie",t,cc,ee,".txt", sep=""), header=TRUE, sep="\t")
      h <- g[g$Goterm!="GO:0008150",] #Remove overarching GO categories biological_process, molecular_function
      h <- h[h$Goterm!="GO:0003674",]

      par(mar = c(4,20,4,2) + 0.1)
      barplot(h$diff, main=paste("GOchart",t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$diff),max(h$diff)))
    }
    dev.off()
  }
}
**********ENDR**********






*******POSTANALYSIS TRICKS********
	Calculate the number of MF and BP GO terms that are enriched in well- and poorly-conserved genomic regions.
	Start in the folder 'pie charts':

t="At";
for e in Geo Pco;
do for h in h l;
  do echo $h$e;
  awk -F$'\t' '$6 >= 1.5' PlotRpieRATIO"$t$h$e".txt | tail -n +2 | wc -l;
  done;
done;









**********TOOLS**********

	Check whether +STATS* files have the same number of lines.
	They should:
wc -l +STATS*;

	Check whether +STATS* files have the same loci in column 2:
g=$(cut -d$'\t' -f2 +STATS.*.RgeoTgen.txt);
p=$(cut -d$'\t' -f2 +STATS.*.RpcoTgen.txt);
diff <(echo "$g") <(echo "$p");


	Check whether SUM.RgeoTgen and SUM.RpcoTgen have the same number of lines.
	They should:

v="1";
#(for i in $v;
(for i in {1..1000}; 
do f=$(wc -l SUM.b$i.R*Tgen.txt| awk '{print $1}' | head -2 | sort -u);
  if (( $(echo "$f" | wc -l) == 2 )); then
    ls -l SUM.b$i.R*Tgen.txt;
    echo;
  fi;
done;) > xlinenum.txt;

#If same (or different), figure out if Rgeo and Rpco list different loci:
#(for i in $v;
(for i in {1..1000}; 
do echo $i;
  g=$(cut -d$'\t' -f2 SUM.b$i.RgeoTgen.txt | sort -u | md5sum);
  p=$(cut -d$'\t' -f2 SUM.b$i.RpcoTgen.txt | sort -u | md5sum);
  if [[ $g != $p ]]; then
    ls -l SUM.b$i.R*Tgen.txt;
    echo;
  fi;
done;) > xdiffloci.txt;

If there are some files that have the same number of lines, but differ in the loci listed
in those lines, figure out where the differences are:

for i in "6";  #cycle through list of offending files
do g=$(cut -d$'\t' -f2 SUM.b$i.RgeoTgen.txt);
  p=$(cut -d$'\t' -f2 SUM.b$i.RpcoTgen.txt);
  diff <(echo "$g") <(echo "$p") > diff$i.txt;
done;




If the same, compute +STATS* files.

If different, determine where the locus names differ between +STATS* files
sg=$(cut -d$'\t' -f2 +STATS.Sorghum.RgeoTgen.txt);
sp=$(cut -d$'\t' -f2 +STATS.Sorghum.RpcoTgen.txt);
diff <(echo "$sg") <(echo "$sp") > f.txt

Determine which blocklengths have a scrambled order relative to what's expected:

#(for i in $v;
(for i in {1..1000}; 
do g=$(cut -d$'\t' -f2 SUM.b$i.RgeoTgen.txt | md5sum);
  gsort=$(cut -d$'\t' -f2 SUM.b$i.RgeoTgen.txt | sort -t'.' -n -k1,1 -k2,2 | md5sum);
  p=$(cut -d$'\t' -f2 SUM.b$i.RpcoTgen.txt | md5sum);
  psort=$(cut -d$'\t' -f2 SUM.b$i.RpcoTgen.txt | sort -t'.' -n -k1,1 -k2,2 | md5sum);
  if [[ $g != $gsort ]]; then
    echo "SUM.b$i.RgeoTgen.txt";
  fi;
  if [[ $p != $psort ]]; then
    echo "SUM.b$i.RpcoTgen.txt";
  fi;
done; ) > xscrambled.txt

#Verify that all loci are included, done on Mac:
(nq=63; #(AtHclust=117,AtKmeans=89 ,AtRcut=27) number of env+geo+pco characters, to be cut off
        #(SorHclust=63,SorKmeans=53)
t="../../script7outputHclust/SorM";
#for i in $v;
for i in {1..1000};
do totloc=$(awk '{print $1}' "$t".b"$i".RgenTenv.var | tail -n +4 | ghead -n -$nq | sort -u);
  g=$(cut -d$'\t' -f2 SUM.b$i.RgeoTgen.txt | tail -n +2 | sort -u);
  p=$(cut -d$'\t' -f2 SUM.b$i.RpcoTgen.txt | tail -n +2 | sort -u);
  if [[ $g != $totloc ]]; then
    echo "problem: SUM.b$i.RgeoTgen.txt";
  fi;
  if [[ $p != $totloc ]]; then
    echo "problem: SUM.b$i.RpcoTgen.txt";
  fi;
done; ) > xallloci.txt

If all loci are included, but order is scrambled, the file can be sorted:

****
myp() {
      i=$1;
      echo $i.$e;
      sortkey1=$(cut -d$'\t' -f2 SUM.b"$i".R"$e"Tgen.txt | cut -d'.' -f3);
      paste -d$'\t' SUM.b"$i".R"$e"Tgen.txt <(echo "$sortkey1") | sort -t$'\t' -n -k28 -n -k3,3 > tmp"$i"."$e".txt 
      #Determine whether blocklengths are scrambled again
      g=$(cut -d$'\t' -f2 tmp"$i"."$e".txt | md5sum);
      gsort=$(cut -d$'\t' -f2 tmp"$i"."$e".txt | sort -t'.' -n -k1,1 -k2,2 | md5sum);
      if [[ $g != $gsort ]]; then
        echo "tmp"$i"."$e".txt is still screwed";
      fi;
}
export -f myp;

slist=$(grep SUM xscrambled.txt | cut -d'.' -f2 | sed 's/b//g' | sort -nu); #get list of scrambled output files
elist="geo pco";
for e in $elist;
  do export e; 
    echo "$slist" | parallel --env myp --env e myp;
  done;
****




#If no longer scrambled, mv back to original file name, eliminating last column (added for sorting):
****
mya() {
      i=$1;
      echo $i.$e; 
      rev tmp"$i"."$e".txt | cut -d$'\t' -f2- | rev > SUM.b"$i".R"$e"Tgen.txt;
      rm tmp"$i"."$e".txt;
}
export -f mya;

slist=$(grep SUM xscrambled.txt | cut -d'.' -f2 | sed 's/b//g' | sort -nu); #get list of scrambled output files
elist="geo pco";
for e in $elist;
  do export e; 
    echo "$slist" | parallel --env mya --env e mya;
  done;
****



#Determine whether all loci contain the same number of reported lines in the SUM file.
#Here, the grep -v " 22 " refers to the number of core sizes searched for Sorghum, 2-23.
#This shows loci where there are not 22 summary lines:

cut -d$'\t' -f2 SUM.b1.RgeoTgen.txt | sort -n -t'.' -k3,3 | uniq -c | grep -v " 21 "; # returns '1 locus' if all loci have the correct number of lines
cut -d$'\t' -f2 SUM.b1.RpcoTgen.txt | sort -n -t'.' -k3,3 | uniq -c | grep -v " 21 ";


#repeat needed analyses:
a=$(cat rgeobad.txt);
for i in $a; do z=$(echo "$i" | sed 's/-/ /g'); ./8dSorghumSlice1RgeoCeresLociSubsetBeginEnd.sh $z; done;
for i in $a; do z=$(echo "$i" | sed 's/-/ /g'); ./8dAtSlice1RpcoLociSubsetBeginEnd.sh $z; done;


#take old sum file, remove bad lines from above:
v="10.948500.375681 10.952995.375683 10.953073.375684 10.953103.375685 1.949028.1064 5.953135.229205"
grep -e 10.948500.375681 -e 10.952995.375683 -e 10.953073.375684 -e 10.953103.375685 -e 1.949028.1064 -e 5.953135.229205 SUM.b1.RgeoTgen.txtORIG > SUM.b1.RgeoTgen.txtORIG2; 


#concatenate
a=$(tail -n +2 newSUM.txt);
cat SUM.b1.RgeoTgen.txtORIG2 <(echo "$a") > SUM.b1.RgeoTgen.txt;
mv SUM.b1.RgeoTgen.txt ../compare;


#find rgeo missing loci by finding non-consecutive locus id number, will not find missing loci at end of file: 

l=$(cut -d$'\t' -f2 SUM.b1.RgeoTgen.txt | cut -d'.' -f3 | sort -nu); conseq "$l" > torepeatrgeo.txt;
l=$(cut -d$'\t' -f2 SUM.b1.RpcoTgen.txt | cut -d'.' -f3 | sort -nu); conseq "$l" > torepeatrpco.txt;
	
	ll=$(cut -d$'\t' -f2 tmp1.geo.txt | cut -d'.' -f3 | uniq)
	
	none









#compute +STATS* files

#additional tools
#determine loci with fewer than expected lines in the SUM file, there should be (#pops - 1)
#lines per locus:

cut -d$'\t' -f2 SUM.b1.RgeoTgen.txt | sort | uniq -c | grep -v "22 " #assumes Sorghum with #pops=23


