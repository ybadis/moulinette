#!/bin/bash

##	Release v11.7 03-12-19
##	Contact yacinebadis@sams.ac.uk


################## DESCRIPTION - USAGE - OUTPUT
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

PATH="/opt/sratoolkit/2.9.6/bin:$PATH"
PATH="/opt/edirect:$PATH"
#PATH="/opt/usearch/9.1.13:$PATH"
PATH="/opt/usearch/9.2.64:$PATH"
#PATH="/opt/usearch/11.0.667:$PATH"
PATH="/opt/ncbi-blast/2.9.0/bin:$PATH"

NTHREADS=16

get_GPS_coordinates()
{
	local COORDINATE=$(esearch -db sra -query $1 | efetch -format native | xtract -pattern SAMPLE_ATTRIBUTES -block SAMPLE_ATTRIBUTE -if TAG -equals "$2" -element VALUE)
	echo $COORDINATE
}


usage() { 
	echo "Usage: $0 [-q <fasta>] [-l <list>] [-t <threshold>] [-n <ntarget>] [-p <percentage>] [-e <evalue>] [-o <dir>] [-m <nthreads>]" 1>&2
	echo "fasta: file e.g. Ectrogella or Olpidiopsis"
	echo "list: name of input list of sra RUN in txt format (e.g runselector output on ncbi, e.g ERR867907 and many others)"
	echo "threshold: blast identity threshold"
	echo "ntarget: max number of target sequences in blast parameter"
	echo "percentage: percentage read cover for blast algntlgt filter expressed in decimal e.g 0.8"
	echo "evalue: Blast evalue parameter"
	echo "dir: output directory"
	echo "Usage example for metabarcoding:"
        echo "./SraSST.sh -q query.fa -l runlist.txt -t 98.5 -n 100 -p 0.8 -e 1e-100 -o out -m 16"
	exit 1
}

while getopts ":h:q:l:t:n:p:e:o:m:" o; do
    case "${o}" in
        q) QUERY=${OPTARG} ;; # fasta file e.g. Ectrogella or Olpidiopsis
        l) LIST=${OPTARG} ;; # name of input list of sra RUN in txt format (e.g runselector output on ncbi, e.g ERR867907)
        t) PERCID=${OPTARG} ;; # blast identity threshold
        n) MAXTARGET=${OPTARG} ;; # max number of target sequences in blast parameter
        p) READCOVER=${OPTARG} ;; # percentage read cover for blast algntlgt filter expressed in decimal e.g 0.8
        e) EVALUE=${OPTARG} ;; #  Blast evalue parameter
        o) OUTPUTDIR=${OPTARG} ; mkdir -p $OUTPUTDIR;; #  Outputdir
        m) NTHREADS=${OPTARG} ;; #  Number of threads
	\?) echo "Invalid option: -$OPTARG" >&2;exit ;;
        *) usage ;;
    esac
done
# shift so that $@, $1, etc. refer to the non-option arguments
shift $((OPTIND-1))



#----------------------------------INPUT VALIDATION----------------------

RUNCOUNT=$(cat "$LIST" | wc -l)
QUERYNB=$(cat "$QUERY" | grep ">" | wc -l)

echo QUERY is "$QUERY"  fasta queryfile
echo LIST is "$LIST"  name of input list of sra RUN 
echo PERCID "$PERCID"  blast identity threshold
echo MAXTARGET is "$MAXTARGET" "max number of target sequences in blast parameter"
echo READCOVER is "$READCOVER"  "percentage read cover for blast algntlgt filtering"
echo EVALUE is "$EVALUE"  "Blast evalue parameter"
echo OUTPUT Folder is "$OUTPUTDIR"
echo ""
echo ""
echo Launching SraSST RUN with "$RUNCOUNT" sra runs and $QUERYNB fasta query sequences
echo In Our hands, a single metagenomic run can typically take 0 to 15 min to scan with 9 queries
echo THIS MIGHT BE LONG ...
echo ""
echo ""

while true;
do
	echo "Please check your input variables."
	echo ""
	echo ""
	echo -n "Do you wish to proceed? Please confirm (y or n) :"
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

start=`date +%s`  #Measuring time

#############################################################################################################################


## Generate uniq SraSST Runid based on date hour min sec - Thi ID will be USED in filename
ID=$(date | awk '{print $6$2$3$4}' | sed 's/://g' )

## CreateDirectoryStructure - GlobalResultFile - GlobalGPScoordinates
#QUERY_f="$(basename -- $QUERY)"
QUERY_f="query"
#LIST_f="$(basename -- $LIST)"
LIST_f="list"
RUN="SraSST-"$ID"-"$QUERY_f"_vs_"$LIST_f
echo $RUN
FOLDER=$OUTPUTDIR"/"$RUN
mkdir -p $FOLDER
echo Creating result directory: $FOLDER

PAIRED_READS_f=$FOLDER/"All-retrieved-reads-paired.fastq"; touch $PAIRED_READS_f
UNPAIRED_READS_f=$FOLDER/"All-retrieved-reads-nonpaired.fastq"; touch $UNPAIRED_READS_f
UNPAIRED_READS_fa=$FOLDER/"All-retrieved-reads-nonpaired.fasta"
PAIRED_READS_R1_f=$FOLDER/"All-retrieved-reads-paired_R1.fastq";
PAIRED_READS_R1_fa=$FOLDER/"All-retrieved-reads-paired_R1.fasta";
PAIRED_READS_R2_f=$FOLDER/"All-retrieved-reads-paired_R2.fastq";
PAIRED_READS_R2_fa=$FOLDER/"All-retrieved-reads-paired_R2.fasta";
UNPAIRED_RETRIEVED_READS=$FOLDER/"All-retrieved-reads-nonpaired.besthits.RblastVsOtu.tab"
PAIRED_RETRIEVED_READS_R1=$FOLDER/"All-retrieved-reads-paired_R1.besthits.RblastVsOtu.tab"
PAIRED_RETRIEVED_READS_R2=$FOLDER/"All-retrieved-reads-paired_R2.besthits.RblastVsOtu.tab"

USEARCH_merged_f=$FOLDER/"usearch.merged.fq";
USEARCH_merged_fa=$FOLDER/"usearch.merged.fa";
USEARCH_merged_RETRIEVED=$FOLDER/"usearch.merged.besthits.RblastVsOtu.tab";
USEARCH_filtered_f=$FOLDER/"usearch.filtered.fa";
USEARCH_uniques_f=$FOLDER/"usearch.uniques.fa";
USEARCH_otus_f=$FOLDER/"usearch.otus.fa";
USEARCH_RETRIEVED=$FOLDER/"usearch.otus.NRBlast.tab";
USEARCH_otus_GIDS=$FOLDER/"usearch.otus.NRBlast.gids";
USEARCH_otus_GIDS_NAME=$FOLDER/"usearch.otus.NRBlast.gids-Namelist.txt";

GPS_f=$FOLDER/"GPS_coordinates.tab"; touch $GPS_f
COUNTS_f=$FOLDER/"counts.tab"; touch $COUNTS_f
GLOBAL_f=$FOLDER/"GlobalResults.tab"; touch $GLOBAL_f
LOG_f=$FOLDER/"logs.log"; touch $LOG_f


echo -e "RUNID\tLATITUDE\tLONGITUDE" >> $GPS_f
echo -e "QUERY\tRUNID\tReadcount\tLATITUDE\tLONGITUDE" >> $COUNTS_f
echo -e "Query\tREADID\tPERCID\tALNLGT\tMISM\tGAP\tQSTART\tQEND\tSbSTART\tSbEND\tEVAL\tBITSCORE\tRUN\tREAD\tLongitude\tLatitude" >> $GLOBAL_f

# Save Original Input for future reference

INPUT_FOLDER=$FOLDER"/input"
mkdir -p $INPUT_FOLDER

cp $QUERY $INPUT_FOLDER"/"
cp $LIST $INPUT_FOLDER"/"


# Fill in run information in .log file
DATE=$(date)
echo $DATE >> $LOG_f
echo "SraSST-"$ID"-"$QUERY"_vs_"$LIST"" >> $LOG_f
echo QUERY is "$QUERY"  fasta queryfile >> $LOG_f
echo LIST is "$LIST"  name of input list of sra RUN  >> $LOG_f
echo PERCID "$PERCID"  blast identity threshold >> $LOG_f
echo MAXTARGET is "$MAXTARGET" "max number of target sequences in blast parameter" >> $LOG_f
echo READCOVER is "$READCOVER"  "percentage read cover for blast algntlgt filtering" >> $LOG_f
echo EVALUE is "$EVALUE"  "Blast evalue parameter" >> $LOG_f
echo "for each $ i sra run (e.g. ERR867907) blastn_vdb command will be" >> $LOG_f
echo "./blastn_vdb -db "$ i" -query $QUERY -outfmt 6 -out ""$QUERY"_vs_"$i".blastres.tab" -perc_identity $PERCID -max_target_seqs $MAXTARGET  -evalue $EVALUE" >> $LOG_f

echo ""
echo ""
echo ""
#
#
#
###Actual start of input processing
for i in $(cat $LIST); 
do
	echo "Processing ""$i"
	INT_FOLDER=$FOLDER/$i
	mkdir $INT_FOLDER
	echo "creating "$INT_FOLDER >> $LOG_f

	##########  BLAST ###############################

	echo "Running blastn_vdb"
	BLAST_f=$INT_FOLDER/blastres.tab
	echo ""##blastn_vdb -db "$i" -query "$QUERY" -outfmt 6 -out "$BLAST_f" -perc_identity "$PERCID" -max_target_seqs "$MAXTARGET" -evalue "$EVALUE""" >> $LOG_f
	blastn_vdb -db "$i" -num_threads $NTHREADS -query $QUERY -outfmt 6 -out "$BLAST_f" -perc_identity $PERCID -max_target_seqs $MAXTARGET  -evalue $EVALUE &> $LOG_f

	##########  READCOVER FILTER  ####################

	# checking read length to adapt blast align fiter
	echo "estimating mean read length on first 100 reads"
	ALNLGT=$(fastq-dump --split-files --skip-technical -X 100 -Z $i | grep "length=" | sed 's/ /\t/g'| cut -f3 | sed 's/length=//g' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'| awk '{print $1 * '"$READCOVER"'}' | cut -f1 -d '.')
	echo "User Setting is" "$READCOVER" ".....SraSST will ONLY select blast alignments over "$ALNLGT"nt" 	
	echo "~"
	echo "~"
	echo "~"
	# Filter results based on blast align filter 
#	#print blastres with aligntlengt filter criteria and extracts hit name
#	(e.g SRA:ERR867907.21038.2) to deduce sra run spot (or read) id number
#	(e.g 21038)
	ALNLGT_f=$INT_FOLDER/"alnlgt_ids"
	cat "$BLAST_f" | awk ' $4 >= '"$ALNLGT"' ' | cut -f2 | tr : '\t' | tr . '\t' | cut -f3 | sort -u > $ALNLGT_f

	RESULT=$(cat "$ALNLGT_f" | wc -l)
	echo "$RESULT" "results to retrieve"



#	############   GPS EXTRACTION  ####################
	if [[ "$RESULT" -gt 0 ]]; then
		LATITUDE=$(get_GPS_coordinates $i 'Latitude');
		LONGITUDE=$(get_GPS_coordinates $i 'Longitude');
		if [ -z $LATITUDE ];
		then
                	LATITUDE=$(get_GPS_coordinates $i 'latitude');
			LONGITUDE=$(get_GPS_coordinates $i 'longitude');
		fi
		
		if [ -z $LATITUDE ];
		then
			LATITUDE=$(get_GPS_coordinates $i 'Latitude Start');
			LONGITUDE=$(get_GPS_coordinates $i 'Longitude Start');
		fi
		if [ -z $LATITUDE ];
		then
			#'geographic location (latitude and longitude)'
			COORDINATES=$(get_GPS_coordinates $i 'lat_lon');
			DECIMAL=$(echo $COORDINATES | awk '{ print $1$2,$3$4 }'| GeoConvert)
			LATITUDE=$(echo $DECIMAL | cut -d' ' -f1 )
			LONGITUDE=$(echo $DECIMAL | cut -d' ' -f2 )
		else
			echo "No GPS coordinates found"
		fi
		
			
		echo "Latitude: "$LATITUDE
		echo "Longitude: "$LONGITUDE

		echo -e "$i"'\t'"$LATITUDE"'\t'"$LONGITUDE" >> $GPS_f
		echo "GPS coordinates stored in "$GPS_f
	fi

	##########   READ EXTRACTION   ##################### 
	for k in $(cat "$ALNLGT_f"); do
		FASTQ_f=$INT_FOLDER/$i.$k.fastq-dump
		
		#Fetch retrieved individual fastq read pairs
		fastq-dump --split-files --skip-technical -N $k -X $k -Z $i > $FASTQ_f

		#Add the above to global result fastq file
		READNUMB=$(cat "$FASTQ_f" | grep $i | grep @ | wc -l)
		if [[ "$READNUMB" -eq  2 ]] ; then
			cat $FASTQ_f >> $PAIRED_READS_f
		else 
			cat $FASTQ_f >> $UNPAIRED_READS_f
		fi
		##Append GPS info to blast results to create the Global Read Detail Table
   		cat "$BLAST_f" | grep "$i"."$k" |  awk ' {print $1"\t"$2"\t" $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t" "'$i'" "\t" "'$k'" "\t" "'$LONGITUDE'" "\t" "'$LATITUDE'"}' >> $GLOBAL_f
		##Disabled here ##### Convert to fasta
		cat $FASTQ_f | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $FASTQ_f.fasta
	done

	### Emptying cache after each run (could be modified to only keep runs matching query)
#	echo "Running cache-mgr"
	cache-mgr -c
	
	###Check if file does not contain any fastq file (e.g. no reads were extracted), and delete if true
	echo $INT_FOLDER
	ls "$INT_FOLDER" | grep fastq-dump
	if [[ $? -ne 0 ]]; then
		#rm -r "$INT_FOLDER"
		echo  "$i"" SORRY!! No hits retrieved in this Run " 
		echo  "$i"" SORRY!! No hits retrieved in this Run " >> $LOG_f
		echo "" 
	else                                        
		echo "$i" "HITs found for this Run!!!"
		echo "$i" "HITs found for this Run!!!" >> $LOG_f
		echo "" 
	
		##Feed estimate of results to global result file
		y="$(ls $INT_FOLDER/*fastq-dump | wc -l)"
		## minus 2 to uncount blastid and blastres files --> count indiv reads fastq files
		#let "y -= 2"
		echo ""$QUERY"_"$i"_"$y"_"$LATITUDE"_"$LONGITUDE"" | sed 's/_/\t/g' >> $COUNTS_f
	fi
done

#UPdate GPS data into 
#join --header -12 -a2 $GPS_f $COUNTS_f -o 2.1,2.2,2.3,1.2 -e >> /home/sa04mp/MoulinetteOutput/COUNT_GPS_join


###################################################   OTU SEARCH   ###########################################################
#
###Manual creation of global fastq file --> creating R1 and R2 file
#
echo "" 
echo "" 
echo "LAUNCHING USEARCH OTU CALLING (Default parameters accept singletons)"
echo "" 
echo "" 


cat $PAIRED_READS_f | awk 'BEGIN {OFS = "\n"} {header1 = $0 ;
 getline seq1 ; getline qheader1 ; getline qseq1; getline header2; getline seq2; getline qheader2; getline qseq2 ;
  {print header1, seq1, qheader1, qseq1}}' > $PAIRED_READS_R1_f

cat $PAIRED_READS_f | awk 'BEGIN {OFS = "\n"} {header1 = $0 ;
 getline seq1 ; getline qheader1 ; getline qseq1; getline header2; getline seq2; getline qheader2; getline qseq2 ;
  {print header2, seq2, qheader2, qseq2}}' > $PAIRED_READS_R2_f


##Merge read pairs
usearch -fastq_mergepairs $PAIRED_READS_R1_f -fastqout $USEARCH_merged_f

#if "Label Mismatch error" the mentionned reads dont have their respective R1
#or R2 mate, you'll have to delete them before OTU calling"

##Standard usearch filtering parameters
usearch -fastq_filter $USEARCH_merged_f -fastq_maxee 1.0 -fastaout $USEARCH_filtered_f

##Sort uniq reads
usearch -fastx_uniques $USEARCH_filtered_f -relabel Uniq -sizeout -fastaout $USEARCH_uniques_f

##Cluster OTUs  - Option -minsize 1 accepts singleton OTUs - radius 1.5 groups together all reads over 99 id
usearch -cluster_otus $USEARCH_uniques_f -otu_radius_pct 1 -minsize 2 -otus $USEARCH_otus_f -relabel SraSST-"$ID"-Otu

echo "Number of OTUs (including singletons) retrieved by USEARCH for run " $FOLDER  
OTUS_RETRIEVED=$(cat $USEARCH_otus_f | grep '>' | wc -l )
echo $OTUS_RETRIEVED

##Add number of dertemined OTUs to Logfile
echo "Number of OTUs (including singletons) retrieved by USEARCH for run "$RUN >> $LOG_f
echo $OTUS_RETRIEVED >> $LOG_f

##Build custom blastdb based on retrieved OTUs
###########Build future query fasta files
cat $UNPAIRED_READS_f | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $UNPAIRED_READS_fa
cat $PAIRED_READS_R1_f | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $PAIRED_READS_R1_fa
cat $PAIRED_READS_R2_f | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $PAIRED_READS_R2_fa
cat $USEARCH_merged_f | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $USEARCH_merged_fa



# create global blastdb directory in the SraSST workspace, will include the db of different runs called by their IDs
DBFOLDER=$FOLDER/SraSSTblastdb
mkdir $DBFOLDER
TITLE=SraSST-"$ID"-otus
BLASTDB_ID=$DBFOLDER/$TITLE

makeblastdb -in $USEARCH_otus_f -dbtype nucl -title $TITLE -out $BLASTDB_ID -parse_seqids


###########Launch blasts with each query file
blastn -num_threads $NTHREADS -db $BLASTDB_ID -query $UNPAIRED_READS_fa -max_target_seqs 1 -outfmt "6" -out $UNPAIRED_RETRIEVED_READS
blastn -num_threads $NTHREADS -db $BLASTDB_ID -query $PAIRED_READS_R1_fa -max_target_seqs 1 -outfmt "6" -out $PAIRED_RETRIEVED_READS_R1
blastn -num_threads $NTHREADS -db $BLASTDB_ID -query $PAIRED_READS_R2_fa -max_target_seqs 1 -outfmt "6" -out $PAIRED_RETRIEVED_READS_R2
blastn -num_threads $NTHREADS -db $BLASTDB_ID -query  $USEARCH_merged_fa -outfmt "6" -max_target_seqs 1 -out $USEARCH_merged_RETRIEVED



##Remote Blast of retrieved OTUs via nr - best hits in nr

echo "Final Remote Blast of retrieved OTUs on NR - This might take a while ... up to 10 minutes"
blastn -db nr -query $USEARCH_otus_f -max_target_seqs 1 -outfmt "6" -out $USEARCH_RETRIEVED -remote

BLASTDB=/media/data/BLASTDB
#blastn -db nr -num_threads $NTHREADS -query $USEARCH_otus_f -max_target_seqs 1 -outfmt "6" -out $USEARCH_RETRIEVED




##Extracting organism names from hits retrieved in NR
cat $USEARCH_RETRIEVED | cut -f2 > $USEARCH_otus_GIDS

for s in $(cat $USEARCH_otus_GIDS); do
 TAXLIN=$(esearch -db nucleotide -query $s | efetch -format xml |xtract -pattern Bioseq -element Textseq-id_accession Seqdesc_title OrgName_lineage Date-std_year | grep "$s")
 echo "$TAXLIN" >> $USEARCH_otus_GIDS_NAME 
 echo "processed "$s""
done

####Warn user about cache size 
end=`date +%s`
runtime=$((end-start))
echo "Moulinette runtime: "$runtime
echo "BEWARE OF CACHE USAGE"
echo "Running cache-mgr ..."
cache-mgr -r
echo "command cache-mgr with option -c to clear you cache - BUT doing so might affect other running instances of SraSST if they share the same ncbi/cache folder."
echo "Your results are stored in "$FOLDER
echo "with the content listed above"
echo "You can find the best NR hits matching your OTUs in " $USEARCH_RETRIEVED
echo "You can find the Organism dESCRIPTION OF THOSE NR hits matching your OTUs in " $USEARCH_otus_GIDS_NAME
echo "You can now MAP your new OTUs to GPS coordinates via reciprocal blast on the retrived reads"
echo "You can also run the tool again using the retrieved unique reads (usearch folder - unique -) to recursively broaden your search"
echo "You're Welcome :-)"
echo ""

exit 0
