# MOULINETTE
MOULINETTE:  SRA Sequence Search Tool (Srasst) 
Pipeline to report OTU ocurrence of an organism

![alt text](https://github.com/ybadis/moulinette/raw/master/images/pipeline.png "Moulinette pipeline")

## How does it work?

A query in fasta (or Multi query fasta file) is run against a list of SRA runs. For each run, the query is blasted (blastn-vdb; SRA-toolkit 2.8.0, http://www.ncbi.nlm.nih.gov/Traces/sra) on the run. Based on blast results, reads matching any query above the selected identity, e-value, and alignment length thresholds are retrieved in fastq (fastq-dump, SRA-toolkit) in separate folders corresponding to each SRA run. A final result fastq file and count file is generated and each dataset is cleared after processing. 

This data reduction approach allows to screen for large amount of data with no storage concerns. In many instances, the MOULINETTE can retrieve the GPS coordinates of any given SRA Dataset, provided that those metadata are correctly assigned by submitters (e.g Ocean Sampling Day dataset). The pipeline accepts single and paired reads datasets, and those are stored in separate files in fastq format. The current version sends only paired-end reads (pooled reads from all screened datasets) to the USEARCH OTU calling pipeline, using a filtering parameter of 1.0 (expected error), an OTU radius of 3 (clustering based on 97% Identity), and excluding singletons for OTU calling. A remote blast is conducted to retrieve the best match of each defined OTU on the NCBI N/R database. The OTUs defined are finally gathered in a blast database, and all reads (paired and single end) are blasted against OTUs to assign them to their best matching OTU with 97% threshold, thus allowing the definition of biogeographical patterns based on GPS coordinate retrieval.


### (a) - Selection of Input SRA Runs

The pipeline can handle broad range of datasets, selected with or without a priori. Because MOULINETTE analyses one SRA run after the other, and keeps only reads satisfying the user-criteria, the input list of SRA Runs can contain irrelevant data (e.g A full SRA Bioproject containing both eukaryotic 18S and bacterial 16s barcoding data while the user is only interested in 18S data) as it will not yield any result and will not be integrated in any kind of normalization, however this will significantly increase the length of a MOULINETTE run. In this paper, it is the case for the selection of 19000 SRA Freshwater datasets: to our knowledge, and as of January 2017, this selection covered all possible 18S freshwater barcoding SRA runs, but also contained 16S data and whole metagenome sequencing. In such instances, including irrelevant data and increasing the length of the MOULINETTE run is quicker than manually sorting only relevant datasets. The selection of various runs from different studies is not common, and has advantages and limitations that are discussed below. Merging MOULINETTE Runs for large input lists: The current version of the pipeline (V11.3) does not handle multithreading and processes only one run at a time. Users can split a list of datasets in subsets of 1000-2000 SRA runs, and launch several MOULINETTE runs in parallel. Independent output result folders can be merged with the MOULINETTE-MergeRunsV4 script, available upon request.

 

### (b) - Selection of Input query sequences

The main current purpose of the MOULINETTE pipeline is to screen metabarcoding data, to detect which SRA raw dataset might contain reads matching an organism of interest. Although many other uses are envisaged (e.g SNP discovery, gene reconstruction in environmental whole DNA sequencing), the current design of the pipeline has been focussed on metabarcoding strategies. In that sense, queries should be limited to classic markers like 18S, 16S, COI and others, and we encourage users to trim queries to suit their needs. For example, the v4 and v9 regions of the 18s marker of an organism of interest can be considered as two separate queries, this would reduce potential false positive originating from reads matching non variable regions of the markers. Alternatively, retrieved reads can be sorted based on the position of the blast alignment (e.g. locating reads matching v4 or v9 regions for a particular full length sequence). Both approaches yield similar results.

### (c) - Stringency Parameters – E-value

We provided the user with the possibility to choose an e-value threshold for retrieval of meaningful matching reads. This ensures flexibility of the tool for future developments. However, selecting an appropriate e-value threshold is tedious and highly dependent on the queries used (length), and likely to be irrelevant with multi-fasta queries. Additionally, because each SRA dataset is considered as an independent blast database through vdb-blast, e-value scores could vary for the same read between two different SRA RUNs (e.g two SRA runs can have very different sizes, in some instances few thousand reads, in other instances millions of reads). We recommend users to set the e-value threshold to a relaxed stringency (e.g 1e-40) and rely on the two parameters discussed below for retrieval of meaningful information. For reference, in the SRA run ERR867676, containing 115662 reads, a read matching a given query at 100% identity typically has a blast e-value below 1e-100. 

### (d) - Stringency Parameters - %ID and READCOVER

The percentage of identity and the “READCOVER” parameter ($3 and $5 in the Usage section above, respectively) are the main parameters to take in account for accurate retrieval of meaningful hits. The READCOVER parameter is crucial: For a value of 0.8, the pipeline retains a read if the length of the aligned region represents more than 80% of the read length. While many reads could match a query at 100% identity on only few nucleotides (false positives), this filtering is performed upon processing of each read blast result, and ensures that only meaningful alignments retained despite the low stringency e-value mentioned above. Stringent blast parameters (e.g 100% Identity) will retrieve the biogeography of a target organism (i.e only one SRA-OTU is retrieved), while relaxed parameters (e.g 97% identity) will allow to explore the associated OTU diversity and Biogeography, although this can vary depending on paired and single read layouts as discussed below.

### (e) - Screening of Paired-end data Vs Single Read Data

The current design of the pipeline uses vdb-blast to harvest reads matching a given query. In the case of paired-end layout, Vdb-blast considers forward and reverse reads separately and provides two blast results for a read pair. For the 18S V4 region amplified by classical Eukaryotic primers, the most variable region is not centred, See figure below. 

![alt text](https://raw.githubusercontent.com/ybadis/moulinette/master/images/18Sv4region.png "18S V4 region amplified by classical Eukaryotic primers,")

 
In the example shown above, the forward read contains most of the variable region, and will be much more sensitive at discriminating species. The reverse read spans a less variable region and will tend to be less specific. In term of identity scores a barcode (final merged read) 98.5% identical to our target organism is in fact made of two reads with 97.3% and 99.6% identity. This notion is important because the MOULINETTE will harvest both F+R reads of a read pair, even if only one read satisfies the user defined stringency criteria. A user threshold of 97% can thus harvest a Reverse read at 97.1% identity that will give a barcode (merged read) below this threshold, for example approx. 95% Identity to our target organism. This can be useful when only few known sequences cluster close from a target organism, for example in the major but distantly related Olpidiopsis clades of this study, where they share an average 90-92% identity on their V4 region. Such phenomenon is not problematic for the cautious user, as the limited number of final OTUs allows for manual curation. The “De-Novo Targeted OTU” strategy of the MOULINETTE typically give a small number of OTUs (e.g 1 to 30 Vs over 1000 for classic metabarcoding) that are reasonably easy to check and process compared with regular OTU strategies. Note that this does not apply to single read layouts, as they represent the equivalent of the merged read discussed above. In conclusion, the %Id stringency parameter does not have the same meaning for paired vs single read layouts. Future implementations of the pipeline will address this potential limitation through the use of the recently released Magic-BLAST, that can generate blast results taking into account simultaneously the two reads of a pair.

### (f) - OTU Calling

Only paired-end reads are used to call OTUs via the Usearch pipeline, but both paired-end and single-end reads are mapped to OTUs. Future developments of the MOULINETTE Pipeline, will aim at integrating both read types prior to OTU calling. The clustering threshold used by Usearch can be modified in the script (otu radius parameter), and should be adapted to the species delimitation offered by a region of interest. In this study, we chose the 97% clustering threshold (out radius 3) mostly as this standard is use throughout many metabarcoding campaigns, and to be coherent with the OTUs generated by the independent Australian Marine Microbe Project and the P.umbilicalis metabarcoding generated at JGI.

### (g) - Read Mapping

The MOULINETTE pipeline uses regular blast based scores to assign reads to their best matching OTU. Although not using the otu_tab command of the Usearch pipeline, our approach is consistent with the guidelines and discussion provided in the Usearch website (http://www.drive5.com/usearch/manual/mapreadstootus.html). Guidelines of the Usearch website suggest to map merged unfiltered reads to account for a larger fraction of the reads (http://www.drive5.com/usearch/manual/cmd_otutab.html), however considering the low number of read per organism for rare obligate parasites, we chose to also map unmerged reads provided that their reciprocal best blast hit in N/R is an Olpidiopsis sequence. This concerns 135 reads out of the 1283 final read pool, we manually checked that no SRA Runs were retained in the final Biogeography dataset based on the sole presence of such unmerged reads.

### (h) - OTU Occurrence Vs Abundance

The MOULINETTE pipeline is aimed at reporting the occurrence of an organism, but not its abundance. Indeed, as OTU calling is performed on a pool of reads from a whole diversity of datasets (potentially different sequencing technologies, DNA extraction protocols, years) complex normalization strategies would be needed to compare abundance between unrelated SRA Runs. As a consequence, the current version of the pipeline and our definition of OTUs is not appropriate to infer diversity metrics.

### (i) - Cross talk effects

As pooled SRA datasets of a large MOULINETTE run could have been generated on different machines - e.g pooling different SRA bioprojects or different years -  cross-talk effects (contamination of read, assignment of reads to the wrong sample during demultiplexing etc..) are intrinsically less likely. However, the cross-talk effects of related SRA runs (e.g. from the same SRA Bioprojects) in a whole Moulinette run could still give rise to spurious low count reads as discussed in http://www.drive5.com/usearch/manual/otu_count_interpret.html

### (j) - Spurious low count OTUs and Limitations of Biogeography

Firstly, and like with any metabarcoding approach, detecting an organism in a location is a supposedly genuine information, while absence of detection is not an information per se. This is especially valid for organisms occupying specialized ecological niches, and appearing in “boom and bust” epidemiologic cycles, as would be the case of the obligate parasites studied in this paper. In classic metabarcoding experiments/campaigns (usually normalized and treated all in the same pipeline) some OTUs may be filtered out in the final analysis during denoizing, and filtering of low count OTUs. However the limited number of OTUs yielded with our approach allows for careful checking of results. In our study, SRA-OTU9 and SRA-OTU12 (Table AS4) might fall in this category of “doubtful” OTUs, because of their low counts and because they do not clearly cluster with known sequences or other OTUs. In the case of SRA-OTU12, with 7 reads mapped across 3 different datasets of two different locations, all contained in the Ocean Sampling Day dataset, making it a likely “doubtful” OTU. Manually checking SRA metadata associated with those runs did not allow to ascertain the existence of cross-talk effects, but SRA-OTU12 is still to be considered “doubtful”. Concerning SRA-OTU9, with 13 reads over 6 runs in 3 locations, the situation is clearer, as it gathers reads of independent SRA Bioprojects (PRJNA311248 and PRJEB8682), thus excluding any cross-talk artefact. As a consequence, we wish to stay careful and not consider them potential Olpidiopsis OTUs. However we kept both OTUs as they originate from the same MOULINETTE run, and because other low count OTUs are clearly genuine information (e.g SRA-OTU11 has only 11 reads in two SRA runs of the same location but is 99.7% identical to the TARA-OTU-239197 as can be seen on Fig. S7). Among SRA runs potentially containing Olpidiopsis sequences, most of them contain 1 to 30 reads. Our Biogeography dataset also keeps datasets that have “Only” one read assigned to our target organisms. Because of 
* the expected scarcity/rarity of this group of obligate parasites,
* the above-mentioned inability to address comparative abundance between unnormalized datasets and 
* because all called OTUs are based on at least two filtered merged reads in the Usearch pipeline, we chose to keep datasets containing only one read matching a given OTU.

## MOULINETTE pipeline: Usage and script

### Installation Requirements

* SRAtoolkit (sratoolkit.2.8.0-centos_linux64)
* entrezdirect tools (September2016 package)
* usearch (usearch v9.1.13_i86linux32, free academic license)

### Usage

```
bash ./MOULINETTE.sh QUERY LIST PERCID MAXTARGET READCOVER EVALUE   
```

* QUERY = fasta query file 
* LIST = name of input list of SRA Runs in txt format (e.g runselector output on ncbi, ERR867907) 
* PERCID = blast identity threshold
* MAXTARGET = maximum number of target sequences in vdb-blast parameter
* READCOVER = percentage read cover for blast Alignment length filtering expressed in decimal (0.8)
* EVALUE = VDB-Blast e-value parameter

### Results

Output is contained in one global Result folder:

- an *input folder containing copies of original input text files 
- a folder containing the whole usearch outputs generated (*-usearch folder, the final otu file has the *otus.fa suffix)
- a folder for each positive run containing (i) blast results (ii) selected read ids (iii) a fastq file for each extracted read 
- a log file summarizing results and parameters used for the run (*.log)
- a global count table summarizing the number of read pairs retrieved for each positive run (*counts.tab)
- a global fastq file pooling all retrieved reads (separation of paired and single reads -> *All-retrieved-reads-paired.fastq and *All-retrieved-reads-nonpaired.fastq)
- Individual Blast results for each reads (Reciprocal Blast on defined OTUs for paired_R1, paired_R2, and single reads)
- a global table containing GPS coordinates of every positive SRA run (when available)

## Authors

* **Yacine Badis**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
