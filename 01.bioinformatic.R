#This script decribe, how to build phyloseq-compotible database from original fastq files. 
#Nearly the same commads can be used for 16s rRNA and diet data
#The major difference regards taxonomic assignations of resulting ASVs
#While the Silva database was used in the case (16S rRNA data) in this case, custom reference
#database was build in the case of COI data (see script Taxonomy_COI.R)

##########################################
#DEMULTIPLEXING###########################
##########################################

####!!!THESE COMMANDS MUST BE RUN IN TERMINAL!!!######

#Define files for demultiplexing ([1] primers in fasta format and [2] matrix of their combinations) 
ADAPT_PATH=$(echo /media/kreising/DATA/data/DIET_MICRO/MICRO_DATA/External_data/primers.fasta)
MATRIX_PATH=$(echo /media/kreising/DATA/data/DIET_MICRO/MICRO_DATA/External_data/matrix.txt)

#######Demultiplex reads based on primer sequences########
cd /media/kreising/DATA/data/DIET_MICRO/MICRO_DATA/00.RAW_DATA
mkdir ../demultiplexed
rm ../demultiplexed/*

#List of fastq files for demultiplexing   
SAMPLES=$(ls -1| grep "R[12].fastq.gz"|sed -r 's/R[12].fastq.gz//'|sort| uniq)
#demultiplexing based on primer sequences
for i in $SAMPLES; do skewer -x $ADAPT_PATH -M $MATRIX_PATH -b -m head -k 35 -d 0 -t 8 "$i"R1.fastq.gz "$i"R2.fastq.gz -o ../demultiplexed/"$i"; done >> ./log

cd ../demultiplexed/
gzip *fastq
cd ..
mkdir 02A.DEMULTI.16s
cd ./demultiplexed

#######Trimm primers from demultiplexed files########
SAMPLES_16S=$(ls -1| grep "assigned-F_"|sed -r 's/-pair[12].fastq.gz//'|sort| uniq)
for i in $SAMPLES_16S; do skewer -x NNNNCCTACGGGNGGCWGCAG -y GACTACHVGGGTATCTAATCC  -m head -k 35 -d 0 -t 8 "$i"-pair1.fastq.gz "$i"-pair2.fastq.gz -o ../02A.DEMULTI.16s/"$i"; done >> ../02A.DEMULTI.16s/log.trimmed.head
cd ../02A.DEMULTI.16s
gzip *fastq

###########################################
#QUALITY FILTERING#########################
###########################################

####!!!THESE COMMANDS MUST BE RUN IN R!!!######

library(dada2)
library(ggplot2)
library(ShortRead)
library(ape)
library(vegan)
library(ggplot2)
library(plyr)
library(DECIPHER)
library(ape)
library(phytools)

setwd("/media/kreising/DATA/data/DIET_MICRO/MICRO_DATA/02A.DEMULTI.16s/")

#List of forward and reverse reads
LIST<-list.files()
F_reads<-LIST[grep("pair1.fastq.gz",LIST)]
R_reads<-LIST[grep("pair2.fastq.gz",LIST)]

#graphical representation of quality profiles
plotQualityProfile(F_reads[100:105],aggregate = TRUE)+ggtitle("Forward reads")
plotQualityProfile(R_reads[100:105],aggregate = TRUE)+ggtitle("Rewerse reads")

sample.names<-gsub("_trus-trimmed-pair1.fastq.gz","",F_reads)
sample.names<-gsub("-assigned-","",sample.names)
filtFs <- paste0(sample.names, "_READ1_filt.fastq.gz")
filtRs <- paste0(sample.names, "_READ2_filt.fastq.gz")

#Quality filtering
for(x in 1:length(F_reads)) {
  print(sample.names[x])
  fastqPairedFilter(c(F_reads[x], R_reads[x]), c(filtFs[x], filtRs[x]),
                    maxN=0, maxEE=1, minQ=2,truncQ=2,
                    compress=TRUE, verbose=TRUE,
                    minLen = c(260,200),truncLen = c(260,200))
}


###############################################
#DADA DENOISING################################
###############################################

#These commands denoise quality-filtered fastq files and build abundance matrix,
#(samples in rows, ASVs in columns)


#List of quality filtered fastq files
fns <- list.files()
fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) 

fnFs <- fastqs[grepl("_READ1_filt.fastq.gz", fastqs)] 
fnRs <- fastqs[grepl("_READ2_filt.fastq.gz", fastqs)] 
sample.names <- gsub("_READ1_filt.fastq.gz","",fnFs)

#fastq dereplication
derepFs <- derepFastq(fnFs,n = 1e+05, verbose=F)
derepRs <- derepFastq(fnRs,n = 1e+05, verbose=F)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#deoising
dadaFs <- dada(derepFs, selfConsist = TRUE,MAX_CONSIST=20)
dadaRs <- dada(derepRs, selfConsist = TRUE,MAX_CONSIST=20)

#merge denoised forward and reverse ASVs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 10,maxMismatch=1,justConcatenate=F)

#abundance matrix
seqtab <- makeSequenceTable(mergers)

#########################################
#Chimeric sequences######################
#########################################

#extraxt ASVs fasta from abundance matrix
FASTA<-DNAStringSet(colnames(seqtab))
names(FASTA)<-colnames(seqtab)
writeFasta(FASTA,"haplo.fasta")

#elimination of chimeric sequences by uchime (Terminal command)
system("usearch8.0.1517_i86linux32 -uchime_ref haplo.fasta -db ~/DB/gold.fasta -nonchimeras haplo.uchime.fasta -strand plus")

############################################
#TAXONOMY###################################
############################################
DB="~/DB/DADA2/silva_nr99_v138_train_set.fa.gz"
FASTA<-readDNAStringSet("haplo.uchime.fasta")
taxa <- assignTaxonomy(as.character(FASTA), DB, multithread=8,minBoot = 80)

############################################
#Create phyloseq object#####################
############################################

#OTU TABLE
seqtab<-otu_table(seqtab,taxa_are_rows = F)

#HAPLO
HAPLO<-readDNAStringSet("haplo.fasta")

#TAXO
TAXO<-tax_table(taxa)

PHYLOSEQ<-merge_phyloseq(seqtab,TAXO,HAPLO)

#############################################
#Add sample metadata#########################
#############################################

MET<-read.delim("/media/kreising/DATA/data/DIET_MICRO/MICRO_DATA/External_data/metadata.txt",
                header = T,
                stringsAsFactors = F)

SN<-paste(MET$Adaptor_F_R_sekvence,"_F_",MET$primerF,"R",sep="")

MET<-sample_data(MET)
sample_names(MET)<-SN

PHYLOSEQ_SD<-merge_phyloseq(PHYLOSEQ,MET)

##############################################
#CONSISTENCY AMONG TECHNICAL DUPLICATES#######
##############################################

#this custom script takes technical duplicates and check which haplotypes are present in both samples
#haplotypes, that are NOT present in both duplicates are eliminated
dupl.concensus<-function(PHYLOS,NAMES){
  
  # exclude nonduplicated samples
  IDS<-as.character(data.frame(sample_data(PHYLOS))[,NAMES])
  IDS.dupl<-IDS[duplicated(IDS)]
  
  PHYLOSEQ<-prune_samples(IDS%in%IDS.dupl, PHYLOS)
  if(length(IDS.dupl)*2<length(IDS)) {NONUPLICATED<-prune_samples(!IDS%in%IDS.dupl, PHYLOS)
  print(paste("Following names are nonduplicated",sample_names(NONUPLICATED)))}
  
  CATS<-as.character(data.frame(sample_data(PHYLOSEQ))[,NAMES])
  CATS2<-levels(factor(CATS))
  OTU_TAB<-otu_table(PHYLOSEQ)
  rownames(OTU_TAB)<-CATS
  
  # i<-5
  for (i in 1:length(CATS2))
  {
    # print(CATS2[i])
    FILTER.act<-colSums(OTU_TAB[rownames(OTU_TAB)==CATS2[i],]>0)>1
    OTU_TAB[rownames(OTU_TAB)==CATS2[i],]
    OTU_TAB[rownames(OTU_TAB)==CATS2[i],]<-t(apply(OTU_TAB[rownames(OTU_TAB)==CATS2[i],],1,function(x) x*FILTER.act))
  }
  
  rownames(OTU_TAB)<-sample_names(PHYLOSEQ)
  otu_table(PHYLOSEQ)<-OTU_TAB
  PHYLOSEQ.clean<-prune_taxa(taxa_sums(PHYLOSEQ)>0,PHYLOSEQ)
  
  PHYLOSEQ.clean
}

#This script merge technical duplicates (specifired in "NAMES" argument)
merge.duplicates<-function(PHYLOSEQ,NAMES){
  CATS<-as.character(data.frame(sample_data(PHYLOSEQ))[,NAMES])
  sample_data(PHYLOSEQ)$duplic.id<-CATS
  SAMDAT<-sample_data(PHYLOSEQ)
  SAMDAT.sub<-subset(SAMDAT,duplicated(CATS)==F)
  FASTA<-refseq(PHYLOSEQ)
  rownames(SAMDAT.sub)<-SAMDAT.sub$duplic.id
  PHYLOSEQ.merge<-merge_samples(PHYLOSEQ,"duplic.id")
  sample_data(PHYLOSEQ.merge)<-SAMDAT.sub
  PHYLOSEQ.merge<-merge_phyloseq(PHYLOSEQ.merge,FASTA)
  PHYLOSEQ.merge
}

CONSIST<-dupl.concensus(PHYLOSEQ_SD,"ID_D")
CONSIST.merged<-merge.duplicates(PHYLOSEQ=CONSIST,NAMES="ID_D") #this is the final phyloseq database
PHYLOSEQ_16S_vlastovky.taxo<-CONSIST.merged

#remove chloroplast, mitochonria unassigned
FILTER<-as.logical(tax_table(PHYLOSEQ_16S_vlastovky.taxo)[,4]!="Chloroplast")
PHYLOSEQ_16S_vlastovky.taxo<-prune_taxa(FILTER,PHYLOSEQ_16S_vlastovky.taxo)
FILTER<-as.logical(tax_table(PHYLOSEQ_16S_vlastovky.taxo)[,5]!="Mitochondria")
PHYLOSEQ_16S_vlastovky.taxo<-prune_taxa(FILTER,PHYLOSEQ_16S_vlastovky.taxo)
FILTER<-as.logical(!is.na(tax_table(PHYLOSEQ_16S_vlastovky.taxo)[,1]))
PHYLOSEQ_16S_vlastovky.taxo<-prune_taxa(FILTER,PHYLOSEQ_16S_vlastovky.taxo)

#Alignment
#accounts for secondary structure of 16S reads
DNA<-DNAStringSet(refseq(PHYLOSEQ_16S_vlastovky.taxo))
RNA<-RNAStringSet(DNA)
ALIGN<-AlignSeqs(RNA,processors = 4)
ALIGN<-DNAStringSet(ALIGN)
writeFasta(ALIGN,file = "/media/kreising/DATA/data/DIET_MICRO/MICRO_DATA/External_data/ASVS_align.fasta")

#phylogeny with FastTree2
system("FastTree -gtr -nt /media/kreising/DATA/data/DIET_MICRO/MICRO_DATA/External_data/ASVS_align.fasta > /media/kreising/DATA/data/DIET_MICRO/MICRO_DATA/External_data/ASVS_align.tree")

#add phylogenty to phyloseq objects
TREE<-read.tree("/media/kreising/DATA/data/DIET_MICRO/MICRO_DATA/External_data/ASVS_align.tree")
TREE<-midpoint.root(TREE)
phy_tree(PHYLOSEQ_16S_vlastovky.taxo)<-TREE
HIRUNDO_16S.taxo.phylo<-PHYLOSEQ_16S_vlastovky.taxo

#rename ASVs names
taxa_names(HIRUNDO_16S.taxo.phylo)<-paste("ASV_", 1:ntaxa(HIRUNDO_16S.taxo.phylo),sep="")

save(HIRUNDO_16S.taxo.phylo, file = "/media/kreising/DATA/data/DIET_MICRO/MICRO_DATA/PHYLOSEQ/HIRUNDO_16S.taxo.phylo.R")
