#!/bin/bash

#start with the following files in the current directory:
#  *.RgenTenv.var files for each blocklength
#  *.b*.envgeopco.dat file for each blocklength 
#  m+1 executable
#  

#  This program creates var files for the following kinds of M+ analyses:
#  Ref:  env pco geo
#  Tar:  gen
#  The genomic targets consist of individual loci with a genome referenced position.
#  The loci are haplotypes with block length from 1-200 consecutive SNPs.

#  M+1 (single processor version) analyses are run in parallel using GNU parallel.
#  Post-M+ analysis of data consists of:
#  1. concatenate M+ output for all loci within a given block link into a single, archive file (RAW file)
#  2. calculate mean, sd, and n for each locus, for all core sizes, for all M+ output
#     columns (these include things like "random reference, optimized reference, random target
#     optimized target, random reference allele counts, optimized...).  Add these to a SUM file.
#     There will be a RAW and SUM file for each block length, for each Reference/Target combination.
#  3. delete all individual M+ output files.
#  4. concatenate all SUM files into a single, unified output table for all blocklengths.


#FUNCTIONS
mycalc() {
    local in="$(echo "$@" | sed -e 's/\[/(/g' -e 's/\]/)/g')";
    awk 'BEGIN {print '"$in"'}' < /dev/null;
}
export -f mycalc;

mysd() {
	#calculate n, mean, and variance/sd using Welford's algorithm and AWK
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
	  print sum/NF "\t" sqrt(var) "\t" NF;
	}' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

locusbylocus() {
	loc=$1;
	ignore=$2;
	target=$3;
	r=$4;
	infile=$5;
	mcom=$6;
	np=$7;
	dird=$8;
	tmpfile=$9;
	b=${10}; #has to be curly-braced or it is treated as $1
	
	echo "Analyzing locus b$b $loc";
	
	s="$loc $ignore";
	repl="$loc $target";
	froot=loc"$loc".R"$r"Tgen;
	
	#save the modified .var file
	varfile=$(echo "$infile" | sed "s:RgenTenv:$froot:g");
	sed "s:$s:$repl:g" $tmpfile > $varfile;


	#set parameters for m+ run
	mpioutfile=$(echo "$varfile" | sed 's:\.var:\.out:g');
	mpidatfile=$(echo "$varfile" | sed 's:\.l.\+:.envgeopco.dat:g');

	#invoke m+ serial version
	$(echo "$dird"/m+1 "$varfile" "$mpidatfile" "$mcom" "$mpioutfile"); #execute bash command normally
	#$echo $dird/m+1 "$varfile" "$mpidatfile" "$mcom" "$mpioutfile" >> t.tmp #just write out the commands
	#condor_run $(echo ./m+1 "$varfile" "$mpidatfile" "$mcom" "$mpioutfile"); #feed bash command into condor queue, this sometimes causes hold when run on head node
	
	#invoke m+ openmp version--causes problems with gnu parallel, have to limit threads or it oversubscribes
	#$(echo ./m+omp "$varfile" "$mpidatfile" "$mcom" "$mpioutfile");

	#invoke m+ mpi version on mac--caches to disk with np=2, mediocre performance
	#$(echo mpirun -np "$np" ./m+mpi "$varfile" "$mpidatfile" "$mcom" "$mpioutfile");

	#invoke m+ mpi version on linux cluster
	#$(echo mpirun -np "$np" -machinefile $HOME/machines ./M+.LINUX.exe "$varfile" "$mpidatfile" "$mcom" "$mpioutfile");

	#POST-PROCESS
	echo "Post-processing Locus b$b $loc";
	#delete the .var file
	rm "$varfile";
	
	#make temporary summary file
	sumtmpfile="$mpioutfile".sumtmp;
	> "$sumtmpfile"; #initialize temporary summary file
	mincore=$(echo $mcom | cut -d' ' -f2);
	maxcore=$(echo $mcom | cut -d' ' -f3);
	step=$(echo $mcom | cut -d' ' -f4);
	
	#cut out each column for each core size, calculate mean, sd, n, write to sumfile
	for ((c=$mincore;c<=$maxcore;c=c+$step));
	do echo -n "$b"$'\t'"$loc"$'\t'"$c"$'\t' >> "$sumtmpfile"; #add blocklength, locus name, and core size to sumfile

	#use if you want the mean of target for all reps, not just the best reps
#		for ((a=2;a<=9;a++));
#		do boo=$(grep ^$c$'\t' $mpioutfile | cut -d$'\t' -f$a  | tr "\n" " "); #extract a column of data, space delimit it
#			rout=$(mysd "$boo"); #calculate mean, sd, n using function mysd()
#			echo -n "$rout" >> $sumtmpfile;
#			if [ a != 9 ]; then echo -n $'\t' >> $sumtmpfile; fi; #add a tab separating results for each column of .out file analyzed, except the last
#		done;
#		echo >> $sumtmpfile; #advance summary to next line, which will be a new core size
	
		#use if you want the mean for best reps only
		boo=$(grep ^$c$'\t' "$mpioutfile"); #collect all 10 reps for the current core size
	
		#find the rep with the best optimized reference diversity (column 3)
		bx=0;
		while read -r k; 
		  do ord=$(echo "$k" | cut -d$'\t' -f3); #get the M+ optimality criterion value for reference
			#determine the best reference value
			if (( $(bc <<< "$ord > $bx") )); then 
			  best="$k";
			  bx="$ord";
			fi;
		  done <<< "$boo"
	  
		#calculate mean values for columns 2-9 for all reps for current core size with the optimal reference diversity
		for ((a=2;a<=9;a++));
		do boo2=$(grep ^$c$'\t'$k "$mpioutfile" | cut -d$'\t' -f$a  | tr "\n" " "); #grep lines of a given core size, extract a column of data, space delimit it
		  rout=$(mysd "$boo2"); #calculate mean, sd, n using function mysd()
		  echo -n "$rout" >> "$sumtmpfile";
		  if [ a != 9 ]; then echo -n $'\t' >> "$sumtmpfile"; fi; #add a tab separating results for each column of .out file analyzed, except the last
		done;
		echo >> "$sumtmpfile"; #advance summary to next line, which will be a new core size
	done;

}
export -f locusbylocus; #necessary for gnu parallel to work


#MAIN

#define basic parameters
np=2; #num processors for mpi run, not used for single processor implementations
dird=$(pwd); #get the current directory to assign path to M+ executable
mincore=2;
maxcore=56; #41,23,56
step=1;
reps=10;
mcom="-m $mincore $maxcore $step $reps";
summarize="yes"; #"yes", "no"
numslice=1; #define the number of chunks to break up the locusnames into
procpernode=24; #number of processors per node for GNU parallel	

#define coding
referencea='2 1 0 1 5';
target='2 0 1 1 5';
ignore='1 0 1 1 5';

#develop parameter sets for haplotype block size analysis
#declare -a parmset1=("$ignore" "$referencea" "$ignore" "$ignore"); #here, ref=env
declare -a parmset2=("$ignore" "$ignore" "$referencea" "$ignore"); #ref=geo
declare -a parmset3=("$ignore" "$ignore" "$ignore" "$referencea"); #ref=pco

#combine
declare -a parsetlist=("${parmset2[@]}" "${parmset3[@]}");
#declare -a parsetlist=("${parmset1[@]}" "${parmset2[@]}" "${parmset3[@]}");

#make var files
p=0; #p is indexed up by $nparms each cycle, this allows you to traverse $parsetlist in chunks corresponding to each $parmset
parsets=2; #number of parameter sets
nparms=4; #number of parameters per parameter set
	for ((j=1;j<=$parsets;j++));
		#define reference, target, or ignore
		do gen="${parsetlist[$p]}";
		envr="${parsetlist[$(( $p+1 ))]}";
		geo="${parsetlist[$(( $p+2 ))]}";
		pco="${parsetlist[$(( $p+3 ))]}";
		
		#determine the var file name specification
		if [ "$referencea" == "$gen" ]; then r="gen";
		elif [ "$referencea" == "$envr" ]; then r="env";
		elif [ "$referencea" == "$geo" ]; then r="geo";
		elif [ "$referencea" == "$pco" ]; then r="pco";
		fi;
		
		#get list of files to modify
		infiles=$(find "$dird" -maxdepth 1 -name '*.RgenTenv.var');
		
		#iterate over input files containing genotypes with different block lengths
		for infile in $infiles; 
			do echo "$infile", Reference="$r";
			
			#make a filename for the .var file in play (for recovery), add a random suffix in case more than one script is running at a time
			tmpstr=$(mktemp -p . | cut -d"/" -f2); #make random file, save its name
			rm $tmpstr; #rm random file
			tmpfile=$infile"."$tmpstr; #assemble modified var file name
			
			#extract the blocklength in play
			b=$(echo $infile | rev | cut -d'/' -f1 | rev | cut -d'.' -f2 | sed 's:b::g');

			
			#convert existing codings in .RgenTenv.var files to something neutral so we can do global changes later
			#search and replace, all gen loci will be "ignore" after this, reference will vary following parsetlist
			sed $'s:\t: :g' $infile |\
			  sed 's:2 1 0 1 5:gengengengengen:g' |\
			   sed 's:2 0 1 1 5:envenvenvenvenv:g' |\
			    sed "s:\(^geo[0-9]\+ \)1 0 1 1 5:\1geogeogeogeogeo:g" |\
			     sed "s:\(^pco[0-9]\+ \)1 0 1 1 5:\1pcopcopcopcopco:g" |\
			       sed "s:gengengengengen:$gen:g" |\
			        sed "s:envenvenvenvenv:$envr:g" |\
			         sed "s:geogeogeogeogeo:$geo:g"  |\
			          sed "s:pcopcopcopcopco:$pco:g" > $tmpfile;
		

			#create array of locus names, and set up to slice the array if necessary
			locusstr=$(grep '\.' "$infile" | cut -d$'\t' -f1 | uniq); #read the locus names into a space delimited string
			read -r -a locusnames <<< $locusstr; #read space delimited string into an array
			len="${#locusnames[@]}"; #get length of array
			chunk=$(( $len / $numslice )); #get the length of each chunk
						
####PARALLELIZE ON LOCI WITHIN BLOCKLENGTH#####

			#make each single locus a target, analyze with m+, choose one of the following parallelization strategies
			
			#serial solution
				#for loc in $locusnames;
				#do locusbylocus "$loc" "$ignore" "$target" "$r" "$infile" "$mcom" "$np" "$dird" "$tmpfile" "$b";
				#done;
			
			#parallel in N-process batches--about same speed as GNU parallel
				#N=24
				#(
				#for loc in $locusnames; do 
				#   ((i=i%N)); ((i++==0)) && wait
				#   locusbylocus "$loc" "$ignore" "$target" "$r" "$infile" "$mcom" "$np" "$dird" &
				#done
				#);
						
			#parallel using gnu parallel
				#headnode only
				#parallel locusbylocus ::: $locusnames ::: "$ignore" ::: "$target" ::: "$r" ::: "$infile" ::: "$mcom" ::: "$np";
				#across all nodes
				#parallel --env locusbylocus --env mysd --sshloginfile ~/machines --jobs 24 locusbylocus ::: $locusnames ::: "$ignore" ::: "$target" ::: "$r" ::: "$infile" ::: "$mcom" ::: "$np" ::: "$dird" ::: "$tmpfile" ::: "$b";

			#parallel using slices of gigantic locusnames array and gnu parallel
			for ((k=0;k<=$numslice-1;k++));
				do
				start=$(( $k * $chunk ));
				if [ $k == $(($numslice-1)) ]; #make a slice of array, add the last element to the final slice
					then slice=("${locusnames[@]:$start:($len)}"); #last slice
					else slice=("${locusnames[@]:$start:$chunk}"); #normal slice
				fi;
	
				parallel --env locusbylocus --env mysd --sshloginfile ~/machines --jobs $procpernode locusbylocus ::: "${slice[@]}" ::: "$ignore" ::: "$target" ::: "$r" ::: "$infile" ::: "$mcom" ::: "$np" ::: "$dird" ::: "$tmpfile" ::: "$b";
	
				done;
			
			#clean up
			rm $tmpfile;
			
#####END PARALLEL#####
			
			
			#PROCESS DATA for current blocklength and reference data type (env,geo,pco)
			if [ $summarize = "yes" ];
			then
				#get list of all outfiles
				outfiles=$(find "$dird" -maxdepth 1 -name "*.b"$b".loc*.R"$r"Tgen.out"  | sort -t'.' -k5 -n); #obtain a list of files, ordered by position in genome, for the current blocklength and reference data type (env,geo,gen)
				#make the raw file archive, will contain all outfiles concatenated
				rawfile="RAW.b"$b".R"$r"Tgen.txt";
				> "$dird"/"$rawfile"; #initialize raw data archive file
				
				#make the summary file, one line per core size, with means and standard deviations calculated across reps, for all columns
				sumfile=$(echo "$rawfile" | sed 's:RAW:SUM:g');
				header=$'blocklength\tlocus\tcoresize\tranrefx\tranrefsd\tn\toptrefx\toptrefsd\tn\trantarx\trantarsd\tn\topttarx\topttarsd\tn\tranrefcntx\tranrefcntsd\tn\toptrefcntx\toptrefcntsd\tn\trantarcnts\trantarcntsd\tn\topttarcntx\topttarcntsd\tn';
				echo "$header" > "$dird"/"$sumfile";
	
				
				#MAKE RAW FILE ARCHIVE and SUMMARY FILE
				echo "Concatenating...";
				for outfile in $outfiles; 

				do echo " $outfile";
					locpos=$(echo "$outfile" | rev | cut -d'/' -f1 | rev | cut -d'.' -f3-5 | sed 's:loc::g'); #extract the position identifier string from the file name (format = chr.midpoint.index)
					
					#make raw file archive by concatenating all individual M+ .out files generated for this blocklength
					sed "s:^:$b\t$locpos\t:g" "$outfile" >> "$dird"/"$rawfile"; #add block length and locus ID to start of line
				
					#make summary file by concatenating all .sumtmp files produced in parallel
					cat "$outfile"".sumtmp" >> "$dird"/"$sumfile";
					
					#clean up
					rm "$outfile";
					rm "$outfile"".sumtmp";

				done; 
			
				#FINISH RAW FILE ARCHIVE
				#get the header from the tmp file, prepend blocklength and locus columns, write to a new file
				header=$(head -1 "$dird"/"$rawfile" | cut -d$'\t' -f3-);
				header=$'blocklength\tlocus\t'"$header";
				echo "$header" > "$dird"/tmp.txt;
			
				#remove all the existing header lines from RAW file, add to tmp.txt, write back to original rawfile.
				grep -v "^.*members$" "$dird"/"$rawfile" >> "$dird"/tmp.txt; #inverse grep to write lines that are not header
				mv "$dird"/tmp.txt "$dird"/"$rawfile";
			
				#settle up
				echo "Done!\n";
			fi;
			

		done; #loop on blocklength

		p=$(( $p+$nparms )); #update index for traversing $parsetlist

	done; #loop on reference data type (env,geo,pco), target is always gen


