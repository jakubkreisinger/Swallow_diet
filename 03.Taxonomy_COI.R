#These scripts decribe, how was the COI reference database buid


#R packages
library(phyloseq)
library(ShortRead)
library(plyr)
library(taxonomizr)
library(dada2)

#Load data
setwd("/media/kreising/DATA/data/DIET_MICRO/DIET_DATA/PHYLOSEQ/")
load("Hirundo_all.R")

#Prepare fasta file with ASVs reference sequences
FASTA<-refseq(Hirundo_all)
writeFasta(FASTA,file = "/media/kreising/DATA/data/DIET_MICRO/DIET_DATA/BLAST/ref.fasta")

#blast hist for ASVs sequences
#assumes that blast is installed and nt database is downloaded
setwd("/media/kreising/DATA/data/DIET_MICRO/DIET_DATA/BLAST/")
system('~/software/bin/blastn -db /media/kreising/DATA/DB/ncbi/nt -negative_gilist /media/kreising/DATA/DB/ncbi/uncultured.gi.txt -query ref.fasta -max_target_seqs 200 -num_threads 8 -out outfile -outfmt "6 qseqid sseqid pident length mismatch gapope qstart qend sstart send evalue bitscore staxids  sseq"')

#filter and dereplicate blast hits
BLAST<-read.delim("/media/kreising/DATA/data/DIET_MICRO/DIET_DATA/BLAST/outfile",header = F,stringsAsFactors = F)
BLAST[,13]<-gsub("-","",BLAST[,13])
BLAST.sub<-BLAST[nchar(BLAST[,13])>300,] #only hits longer than 300bp retained
BLAST.sub<-BLAST.sub[BLAST.sub$V3>90,] #only hits with > 90% similarity to given ASVs

BLAST.sub.uniq<-BLAST.sub[duplicated(BLAST.sub[,2])==F,] #dereplicate blast hits

#downloading names, nodes and accession2taxid data from NCBI
prepareDatabase('accessionTaxa.sql')

#extracting and dereplicating tax ids from blast results
TAXA<-BLAST.sub.uniq[,12]
TAXA<-as.numeric(sapply(strsplit(TAXA, ";"), function(x) x[1], simplify=TRUE))
TAXA<-unique(TAXA)

#get hierachical taxonomy for each unique tax id
CLASS<-getTaxonomy(TAXA,'/media/kreising/DATA/accessionTaxa.sql')
TAXID<-rownames(CLASS)
CLASS.df<-data.frame(taxid=as.numeric(TAXID),CLASS)

#merge taxonomy with blast result
names(BLAST.sub.uniq)[12]<-"taxid"
BLAST.sub.uniq.tax<-join(BLAST.sub.uniq,CLASS.df)

#collapse all taxomonic levels for individual blast hits into semicolon-separated string 
#(i.e. names compatible with RDP classifier)
CLASS.vect<-apply(BLAST.sub.uniq.tax[,14:dim(BLAST.sub.uniq.tax)[2]],1,function(x) paste(x,collapse = "; "))
BLAST.sub.uniq.tax$taxcactor<-CLASS.vect

#remove blast hits, where the the genus-level taxonomy is not resolved 
BLAST.sub.uniq.tax.sub<-BLAST.sub.uniq.tax[!is.na(BLAST.sub.uniq.tax$genus),]

#preprare and save reference fasta
FASTA<-DNAStringSet(BLAST.sub.uniq.tax.sub$V13)
names(FASTA)<-BLAST.sub.uniq.tax.sub$taxcactor
writeFasta(FASTA,"/media/kreising/DATA/data/DIET_MICRO/DIET_DATA/BLAST/REF.fasta")

#Taxonomic classification by RDP classifier
TO_CLASS<-refseq(Hirundo_all)
taxa<- assignTaxonomy(as.character(TO_CLASS),
                      "/media/kreising/DATA/data/DIET_MICRO/DIET_DATA/BLAST/REF.fasta",
                      multithread=TRUE,minBoot = 80,tryRC=T)

