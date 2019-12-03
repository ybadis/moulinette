# Changelog

All changes to this software will be documented in this file

## v11.5 

- updated GPS extraction, now accepting multiple formats using the dedicated xtract xml parser utility
- only 100 first reads used for computing read length

##  v11.3 

- maxtargetseq set to 1 in the final NR blast to have only the best hit retrieved and facilitate sorting which OTU are of interest. Users will run global remote blasts themselves using maxtargetseq 10 to get the diversity of best matches
- clustering OTU with Minsize 2 --> Excluding singletons
- took out many * wilcards in the cleanup stage so that cp and mv are more precise, allowing running several instances of the moulinette in the same working directory

## v11.2 

- added "--skip-technical to all fastq-dump commands to keep only Biological Reads"
- added "n" to the sort command in the definition of ALNLGT lgt

## v11

- test each retrieved read to see if paired or not to send to two global file, one for single reads, one for paired reads, only the paired reads will be sent through OTU calling
- added a command to save all stdout to final logfile for future debugging


## v10.1 

- add "sort -u" before each fastq-dump to avoid extractin same read several times

## v10.0 

- Implementing final reciprocal blast on of retrieved reads VS usearch Otus - see line 324, reciprocal blast will be made on R1 and R2 reads, but also on merged reads produced by usearch before the filtration step
- also includes remote blast of final OTU list
- I also took outh the relabel option of usearch for enhanced trancability

## v9.2 --> FDB MOD

- ALNLGT variable definition changed;solving the issue of wrapped vs unwrapped read data by grepping "lenght=" of reads using fastq-dump instead of lenght of lines
- changed filenames of executables, removing the version IDs find-replace, cach-mgr to ./cache-mgr to use the working directory

## v9 

-Retrieving GPS before fastqdump to generate per-read global result file
## v8
- GPS- Extraction at each Run and not the end of run. -- v7: addition of 
- SraSST-GPS-Runextract module and Osearch OTU calling to SraSST-v6.6)



