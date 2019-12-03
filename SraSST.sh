#!/bin/bash
##	Release v11.3 10-04-17
##	Contact yacinebadis@sams.ac.uk

########################################################################  
#v11   test each retrieved read to see if paired or not to send to two global file, one for single reads, one for paired reads, only the paired reads will be sent through OTU calling
#v11  added a command to save all stdout to final logfile for future debugging
#v11.2 also added "--skip-technical to all fastq-dump commands to keep only Biological Reads"
#v11.2 added "n" to the sort command in the definition of ALNLGT lgt
#v11.3 maxtargetseq set to 1 in the final NR blast to have only the best hit retrieved and facilitate sorting which OTU are of interest. Users will run global remote blasts themselves using maxtargetseq 10 to get the diversity of best matches
#v11.3 clustering OTU with Minsize 2 --> Excluding singletons
#v11.3 took out many * wilcards in the cleanup stage so that cp and mv are more precise, allowing running several instances of the moulinette in the same working directory

#v10.1 add "sort -u" before each fastq-dump to avoid extractin same read several times
#
#
# v10.0 - Implementing final reciprocal blast on of retrieved reads VS usearch Otus - see line 324, reciprocal blast will be made on R1 and R2 reads, but also on merged reads produced by usearch before the filtration step
# also includes remote blast of final OTU list
#	I also took outh the relabel option of usearch for enhanced trancability
#v9.2 --> FDB MOD
# 17-01-2017 - 	ALNLGT variable definition changed;solving the issue of wrapped vs unwrapped read data by grepping "lenght=" of reads using fastq-dump instead of lenght of lines
#				changed filenames of executables, removing the version IDs find-replace, cach-mgr to ./cache-mgr to use the working directory
### v9 Retrieving GPS before fastqdump to generate per-read global result file
##   #v8: GPS- Extraction at each Run and not the end of run. -- v7: addition of 
##   #SraSST-GPS-Runextract module and Osearch OTU calling to SraSST-v6.6)



################## DESCRIPTION - USAGE - OUTPUT
##################
## Description: A Query $1 is run against  A list of SRA Runs $2. For each run, the query is vdb-blasted on the Run and Reads above the selected 
## alignment length threshholds are retrieved in fastq and converted to fasta in separate folders corresponding to each SRA Run.
## A Final result fastq file and count file is generated. Runs with no hits should not be saved. 

## Requirements: sratoolkit (tested with sratoolkit.2.8.0-centos_linux64); entrezdirect (tested with Sept 2016 package); usearch (tested with usearch v9.1.13_i86linux32)

##### Important  Note: Stream of data does not mean data are not in cache. Up to a certain extent, keeping cache data allows 
##### extra fast retrieval on RUNs that were already queried --> interesting when automated with Many queries following intial scan of a bif runlist file.
##### Data are stored in the ncbi/cache folder,  attention should be paid tokeep the cache folder under control. Per default, cache is emptied after each
##### run, this can be modified by switching cache option -c (clear) to -r (read) near line 264. 

### Usage:  bash ./SraSST-v7.sh $1 $2 $3 $4 $5 $6   
### Usage example for metabarcoding runs: bash ./SraSST-v7.sh query.fasta runlist.txt 98.5 100 0.8 1e-100
### Query can be a multifasta

###OUTPUT will contain In One global SraSST run folder:
### - an *input folder containing copies of original input files 
### - a log file summarizing results and parameters used for the run (*.log)
### - a global count table summarizing the number of read pairs retrieved for each positive run (*counts.tab)
### - a global fastq file pooling all retrieved reads (*All-retrieved-reads-paired.fastq)
### - a global table containg GPS coordinates of every positive sra run (If such metadata were provided by submitters to sra)
### - a folder containg the whole usearch outputs generated (*-usearch folder, the final otu file has the *otus.fa suffix)
### - a folder for each positive run containing (i) blast results (ii) selected ids (iii) a fastq file for each extracted read pair

##Enter Debug mode (deactivated here)
#set -x

QUERY=$1 # fasta file e.g. Ectrogella or Olpidiopsis
LIST=$2 # name of input list of sra RUN in txt format (e.g runselector output on ncbi, e.g ERR867907 and many others)
PERCID=$3 # blast identity threshold
MAXTARGET=$4 #max number of target sequences in blast parameter
READCOVER=$5 # percentage read cover for blast algntlgt filter expressed in decimal e.g 0.8
EVALUE=$6 #  Blast evalue parameter

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#####-----------------------------------------------------------------INPUT VALIDATION----------------------
#####     IN RECONSTRUCTION
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################


##### -----------------------------------------------------------------LAUNCH - SraSST -- warning on run duration
RUNCOUNT=$(cat "$LIST" | wc -l)
QUERYNB=$(cat "$QUERY" | grep ">" | wc -l)

echo Launching SraSST RUN with "$RUNCOUNT" sra runs and $QUERYNB fasta query sequence-s
echo In Our hands, a single metagenomic run can typically take 0 to 15 min to scan with 9 queries
echo THIS MIGHT BE LONG
echo ""
echo ""
echo ""
echo QUERY is "$QUERY"  fasta queryfile
echo LIST is "$LIST"  name of input list of sra RUN 
echo PERCID "$PERCID"  blast identity threshold
echo MAXTARGET is "$MAXTARGET" "max number of target sequences in blast parameter"
echo READCOVER is "$READCOVER"  "percentage read cover for blast algntlgt filtering"
echo EVALUE is "$EVALUE"  "Blast evalue parameter"
echo ""
echo ""
echo ""
while true;
do
	echo "YOUR Input order DOES COUNT !!! Check that your numerical variables are properly assigned and that no inversions were made by improper order of inputs."
	echo ""
	echo ""
	echo ""
	echo -n " Do you wish to proceed? Please confirm (y or n) :"
	read CONFIRM
	case $CONFIRM in
		y|Y|YES|yes|Yes) break ;;
n|N|no|NO|No)
echo Aborting
exit 1
;;
*) echo Please enter only y or n
esac
done
echo Continuing ...

echo ""
echo ""
echo ""
echo ""
##### -----------------------------------------------------------------
##### -----------------------------------------------------------------
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

#Generate uniq SraSST Runid based on date hour min sec - Thi ID will be USED in filename
ID=$(date | awk '{print $6$2$3$4}' | sed 's/://g' )

##CreateDirectoryStructure - GlobalResultFile - GlobalGPScoordinates
mkdir "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
	echo "creating result directory " "SraSST-"$ID"-"$QUERY"_vs_"$LIST""

touch  "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired.fastq"
touch  "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-nonpaired.fastq"

touch "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_GPS_coordinates.tab"
	echo "RUN_LATITUDESTART_LONGITUDESTART" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_GPS_coordinates.tab"

touch "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_counts.tab"
	echo ""QUERY_RUN_Readcount_LATITUDE_LONGITUDE"" | sed 's/_/\t/g' >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_counts.tab"

touch "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_GlobalResults.tab"
#queryId, srareadID, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore , RUN, READ, LONGITUDE, LATITUDE, SraSSTRUNID
	echo Query_READID_PERCID_ALNLGT_MISM_GAP_QSTART_QEND_SbSTART_SbEND_EVAL_BITSCORE_RUN_READ_LONGITUDE_LATITUDE | sed 's/_/\t/g' >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_GlobalResults.tab"

#Save Original Input for future reference
mkdir "SraSST-"$ID"-input"
	cp "$QUERY" "SraSST-"$ID"-input"
	cp "$LIST" "SraSST-"$ID"-input"
	mv "SraSST-"$ID"-input" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
#Create Logfile for future reference
touch "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
DATE=$(date)
		# Fill in run information in .log file
		echo $DATE >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
		echo "SraSST-"$ID"-"$QUERY"_vs_"$LIST"" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
		echo QUERY is "$QUERY"  fasta queryfile >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
		echo LIST is "$LIST"  name of input list of sra RUN  >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
		echo PERCID "$PERCID"  blast identity threshold >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
		echo MAXTARGET is "$MAXTARGET" "max number of target sequences in blast parameter" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
		echo READCOVER is "$READCOVER"  "percentage read cover for blast algntlgt filtering" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
		echo EVALUE is "$EVALUE"  "Blast evalue parameter" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
		echo "for each $ i sra run (e.g. ERR867907) blastn_vdb command will be" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
		echo "./blastn_vdb -db "$ i" -query $QUERY -outfmt 6 -out ""$QUERY"_vs_"$i".blastres.tab" -perc_identity $PERCID -max_target_seqs $MAXTARGET  -evalue $EVALUE" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"

		echo ""
		echo ""
		echo ""



##Actual start of input processing
for i in $(cat $LIST); 

do

	echo "Processing ""$i"
	mkdir "SraSST-"$ID"-"$QUERY"_vs_"$i""
	echo "creating ""SraSST-"$ID"-"$QUERY"_vs_"$i""
#############################################################################################################################
#############################################################################################################################
##################################################   BLAST   ################################################################
#############################################################################################################################
#############################################################################################################################
##RunBlast
echo "Running blastn_vdb"
echo ""##blastn_vdb -db "$i" -query "$QUERY" -outfmt 6 -out ""$QUERY"_vs_"$i".blastres.tab" -perc_identity "$PERCID" -max_target_seqs "$MAXTARGET" -evalue "$EVALUE"""
./blastn_vdb -db "$i" -query $QUERY -outfmt 6 -out ""$QUERY"_vs_"$i".blastres.tab" -perc_identity $PERCID -max_target_seqs $MAXTARGET  -evalue $EVALUE

#############################################################################################################################
#############################################################################################################################
##################################################   READCOVER FILTER  ######################################################
#############################################################################################################################
#############################################################################################################################
## checking read length to adapt blast align fiter
echo "estimating mean read length on first 10000 reads"
	ALNLGT=$(./fastq-dump --split-files --skip-technical -X 10000 -Z $i | grep "length=" | sed 's/ /\t/g'| cut -f3 | sed 's/length=//g' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'| awk '{print $1 * '"$READCOVER"'}' | cut -f1 -d '.')
	echo "User Setting is" "$READCOVER" ".....SraSST will ONLY select blast alignments over "$ALNLGT"nt" 	
	echo "~"
	echo "~"
	echo "~"
##Filter results based on blast align filter 
	#print blastres with aligntlengt filter criteria and extracts hit name (e.g SRA:ERR867907.21038.2) to deduce sra run spot (or read) id number (e.g 21038)
	cat "${QUERY}_vs_$i.blastres.tab" | awk ' $4 >= '"$ALNLGT"' ' | cut -f2 | tr : '\t' | tr . '\t' | cut -f3 | sort -u > ""$QUERY"_vs_"$i"_alnlgt_ids"

	

	RESULT=$(cat ""$QUERY"_vs_"$i"_alnlgt_ids" | wc -l)
	echo "$RESULT" "results to retriev

	if [[ "$RESULT" -gt 0 ]]; then 
		#############################################################################################################################
		#############################################################################################################################
		##################################################   GPS EXTRACTION  ########################################################
		#############################################################################################################################
		#############################################################################################################################
		LONGITUDE=$(./esearch -db sra -query $i | ./efetch -format native | sed 's/\///g'| sed 's/></\n/g'|uniq | grep -A2 Longitude | grep -A1 Start| grep VAL | sed 's/VALUE>//g'| sed 's/<VALUE//g')
		LATITUDE=$(./esearch -db sra -query $i | ./efetch -format native | sed 's/\///g'| sed 's/></\n/g'|uniq | grep -A2 Latitude | grep -A1 Start| grep VAL | sed 's/VALUE>//g'| sed 's/<VALUE//g')
		echo "$i"_"$LATITUDE"_"$LONGITUDE" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_GPS_coordinates.tab"
		echo "GPS coordinates stored in " "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_GPS_coordinates.tab"
	fi

#############################################################################################################################
#############################################################################################################################
##################################################   READ EXTRACTION   ######################################################
#############################################################################################################################
#############################################################################################################################
	for k in $(cat ""$QUERY"_vs_"$i"_alnlgt_ids"); do
	##Fetch retrieved individual fastq read pairs
		./fastq-dump --split-files --skip-technical -N $k -X $k -Z $i > "SraSST-"$ID"-"$QUERY"_vs_"$i""/""$i"."$k".fastq-dump"
	##Add the above to global result fastq file
		READNUMB=$(cat "SraSST-"$ID"-"$QUERY"_vs_"$i""/""$i"."$k".fastq-dump" | grep $i | grep @ | wc -l)
	if [[ "$READNUMB" -eq  2 ]] ; then
		cat "SraSST-"$ID"-"$QUERY"_vs_"$i""/""$i"."$k".fastq-dump" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired.fastq"
	else 
		cat "SraSST-"$ID"-"$QUERY"_vs_"$i""/""$i"."$k".fastq-dump" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-nonpaired.fastq"
	fi
	##Append GPS info to blast results to create the Global Read Detail Table
	cat "${QUERY}_vs_$i.blastres.tab" | grep "$i"."$k" |  awk ' {print $1"\t"$2"\t" $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t" "'$i'" "\t" "'$k'" "\t" "'$LONGITUDE'" "\t" "'$LATITUDE'"}' >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_GlobalResults.tab"
	##Disabled here ##### Convert to fasta
		#cat "SraSST-"$ID"-"$QUERY"_vs_"$i""/""$i"."$k".fastq-dump" | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > "SraSST-"$ID"-"$QUERY"_vs_"$i""/"""$i"."$k".fastq-dump.fasta""



	done

## Emptying cache after each run (could be modified to only keep runs matching query)
echo "Running cache-mgr"
cache-mgr -c

##Move blastres table and sorted ids to corresponding RUN folder
mv ""$QUERY"_vs_"$i"_alnlgt_ids" "SraSST-"$ID"-"$QUERY"_vs_"$i""

mv "${QUERY}_vs_$i.blastres.tab" "SraSST-"$ID"-"$QUERY"_vs_"$i""

##Check if file does not contain any fastq file (e.g. no reads were extracted), and delete if true
ls SraSST-"$ID"-"$QUERY"_vs_"$i" | grep *fastq-dump
if [[ $? -ne 0 ]]; then
	rm -r "SraSST-"$ID"-"$QUERY"_vs_"$i""
	echo  "$i"" SORRY!! No hits retrieved in this Run " 
	echo  "$i"" SORRY!! No hits retrieved in this Run " >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
	echo "" 
	echo "" 
	echo "" 

else                                        
	echo "$i" "HITs found for this Run!!!"
	echo "$i" "HITs found for this Run!!!" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
	echo "" 



##Feed estimate of results to global result file
y="$(ls "SraSST-"$ID"-"$QUERY"_vs_"$i"" | wc -l)"
	## minus 2 to uncount blastid and blastres files --> count indiv reads fastq files
	let "y -= 2"
	echo ""$QUERY"_"$i"_"$y"_"$LATITUDE"_"$LONGITUDE"" | sed 's/_/\t/g' >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_counts.tab"

##MoveIndividualRunFolder to main Search directory
mv "SraSST-"$ID"-"$QUERY"_vs_"$i"" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""

fi

done

#############################################################################################################################
#############################################################################################################################
##################################################   OTU SEARCH   ###########################################################
#############################################################################################################################
#############################################################################################################################

##Manual creation of global fastq file --> creating R1 and R2 file

echo "" 
echo "" 
echo "" 
echo "LAUNCHING USEARCH OTU CALLING (Default parameters accept singletons)"
echo "" 
echo "" 
echo "" 

cat "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired.fastq" | awk 'BEGIN {OFS = "\n"} {header1 = $0 ;
 getline seq1 ; getline qheader1 ; getline qseq1; getline header2; getline seq2; getline qheader2; getline qseq2 ;
  {print header1, seq1, qheader1, qseq1}}' > "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R1.fastq"

cat "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired.fastq" | awk 'BEGIN {OFS = "\n"} {header1 = $0 ;
 getline seq1 ; getline qheader1 ; getline qseq1; getline header2; getline seq2; getline qheader2; getline qseq2 ;
  {print header2, seq2, qheader2, qseq2}}' > "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R2.fastq"


 ##Merge read pairs
./usearch -fastq_mergepairs "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R1.fastq" -fastqout "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.merged.fq"
## if "Label Mismatch error" the mentionned reads dont have their respective R1 or R2 mate, you'll have to delete them before OTU calling"

##Standard usearch filtering parameters
./usearch -fastq_filter "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.merged.fq" -fastq_maxee 1.0  -fastaout "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.filtered.fa"

##Sort uniq reads
./usearch -fastx_uniques "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.filtered.fa" -relabel Uniq -sizeout -fastaout "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.uniques.fa"

##Cluster OTUs  - Option -minsize 1 accepts singleton OTUs - radius 1.5 groups together all reads over 97 id
./usearch -cluster_otus "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.uniques.fa" -otu_radius_pct 3 -minsize 2 -otus "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.otus.fa" -relabel SraSST-"$ID"-Otu

echo "Number of OTUs (including singletons) retrieved by USEARCH for run" "SraSST-"$ID"-"$QUERY"_vs_"$LIST"" 
OTUS_RETRIEVED=$(cat "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.otus.fa" | grep '>' | wc -l )
echo $OTUS_RETRIEVED

##Add number of dertemined OTUs to Logfile
echo "Number of OTUs (including singletons) retrieved by USEARCH for run" "SraSST-"$ID"-"$QUERY"_vs_"$LIST"" >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"
echo $OTUS_RETRIEVED >> "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log"

##Build custom blastdb based on retrieved OTUs
##########Build future query fasta files
cat "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-nonpaired.fastq" | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-nonpaired.fasta"
cat "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R1.fastq" | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R1.fasta"
cat "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R2.fastq" | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R2.fasta"
cat "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.merged.fq" | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_u_search.merged.fasta"

mkdir SraSSTblastdb # create global blastdb directory in the SraSST workspace, will include the db of differnt runs called by their IDs

makeblastdb -in "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.otus.fa" -dbtype nucl -title SraSST-"$ID"-otus -out SraSSTblastdb/SraSST-"$ID"-otus -parse_seqids

##########Launch blasts with each query file
blastn -db SraSSTblastdb/SraSST-"$ID"-otus -query "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-nonpaired.fasta" -max_target_seqs 1 -outfmt "6" -out "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-nonpaired.besthits.RblastVsOtu.tab"
blastn -db SraSSTblastdb/SraSST-"$ID"-otus -query "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R1.fasta" -max_target_seqs 1 -outfmt "6" -out "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R1.besthits.RblastVsOtu.tab"
blastn -db SraSSTblastdb/SraSST-"$ID"-otus -query "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R2.fasta" -max_target_seqs 1 -outfmt "6" -out "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R2.besthits.RblastVsOtu.tab"
blastn -db SraSSTblastdb/SraSST-"$ID"-otus -query  "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_u_search.merged.fasta" -outfmt "6" -max_target_seqs 1 -out "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_u_search.merged.besthits.RblastVsOtu.tab"

#Remote Blast of retrieved OTUs via nr - best hits in nr
echo ""
echo ""
echo ""
echo "Final Remote Blast of retrieved OTUs on NR - This might take a while ... up to 10 minutes"
blastn -db nr -query "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.otus.fa" -max_target_seqs 1 -outfmt "6" -out "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_u_search.otus.NRBlast.tab" -remote


##Usearch Cleanup
mkdir "SraSST-"$ID"-usearch"
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R1.fastq" "SraSST-"$ID"-usearch"
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R2.fastq" "SraSST-"$ID"-usearch"
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R1.fasta" "SraSST-"$ID"-usearch"
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R2.fasta" "SraSST-"$ID"-usearch"
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_u_search.merged.fasta" "SraSST-"$ID"-usearch"
mv SraSST-"$ID"-"$QUERY"_vs_"$LIST"_usearch.* "SraSST-"$ID"-usearch" 

##### FINAL CLEANUP
##Redirecting all final results in global result directory
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_u_search.otus.NRBlast.tab" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""

##  mv *RblastVsOtu* "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-nonpaired.besthits.RblastVsOtu.tab" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R1.besthits.RblastVsOtu.tab" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired_R2.besthits.RblastVsOtu.tab" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_u_search.merged.besthits.RblastVsOtu.tab" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""


mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-paired.fastq" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-nonpaired.fastq" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_All-retrieved-reads-nonpaired.fasta" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_counts.tab" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_GPS_coordinates.tab" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST".log" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
mv "SraSST-"$ID"-"$QUERY"_vs_"$LIST"_GlobalResults.tab" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
mv "SraSST-"$ID"-usearch" "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
###Warn user about cache size 
echo "BEWARE OF CACHE USAGE"
echo "Running cache-mgr"
cache-mgr -r
echo "command cache-mgr with option -c to clear you cache - BUT doing so might affect other running instances of SraSST if they share the same ncbi/cache folder"
echo ""
echo "Your Whole results are stored in " "SraSST-"$ID"-"$QUERY"_vs_"$LIST""
echo "with the following content:"
ls "SraSST-"$ID"-"$QUERY"_vs_"$LIST""

echo "You can now MAP your new OTUs to GPS coordinates via reciprocal blast on the retrived reads"
echo "You can also run the tool again using the retrieved unique reads (usearch folder - unique -) to recursively broaden your search"
echo "You're Welcome :-)"

exit 0
