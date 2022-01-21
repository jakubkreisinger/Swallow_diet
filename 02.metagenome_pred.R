library(phyloseq)
library(ShortRead)

#################################################
#1. Prepare filies for picrust2 and run pipeline#
#################################################

#load phyloseq object with 16S rRNA data
load("/media/kreising/DATA/data/VLASTOVKY/PHYLOSEQ_rerun/HIRUNDO_16S.taxo.phylo.R")

#save ASVs sequences as fasta
REFS<-refseq(HIRUNDO_16S.taxo.phylo)
writeFasta(REFS,file = "/media/kreising/DATA/data/VLASTOVKY/PICRUST2/REF.fasta")

#save ASVs abundance matrix as a tab delim file
OTU_TAB<-t(otu_table(HIRUNDO_16S.taxo.phylo))
OTU_TAB<-data.frame("#OTU ID"=taxa_names(HIRUNDO_16S.taxo.phylo),OTU_TAB)
write.table(OTU_TAB,file =  "/media/kreising/DATA/data/VLASTOVKY/PICRUST2/OTU_tab.txt",row.names = F,quote = F,sep="\t")

#run picrust2 pipeline a add descriptions to resulting feature abundance matrix
#!!!must be run in terminal!!!#### 
picrust2_pipeline.py -s REF.fasta -i OTU_tab.txt -o picrust2_out_pipeline -p 6
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
# -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
# add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
# -o pathways_out/path_abun_unstrat_descrip.tsv.gz

##########################################
#Convert EC predictions to phyloseq#######
##########################################

EC<-read.delim("/media/kreising/DATA/data/VLASTOVKY/PICRUST2/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv")

#convert descriptions to the tax_table
EC.tax<-EC[,1:2]
EC.tax<-tax_table(EC.tax)
EC.names<-EC[,1]
taxa_names(EC.tax)<-EC.names

#convert the rest of the table to the otu_table()
EC.otu<-EC[,3:dim(EC)[2]]
EC.otu<-round(EC.otu,0)
EC.otu<-otu_table(EC.otu,taxa_are_rows = T)
taxa_names(EC.otu)<-EC.names

#merge these objects with original metadata and save
HIRUNDO_picrust<-merge_phyloseq(EC.otu,EC.tax,sample_data(HIRUNDO_16S.taxo.phylo))
save(HIRUNDO_picrust,file = "/media/kreising/DATA/data/VLASTOVKY/PHYLOSEQ_rerun/HIRUNDO_picrust.R")

##########################################
#Convert pathways predictions to phyloseq#
##########################################

EC<-read.delim("/media/kreising/DATA/data/VLASTOVKY/PICRUST2/picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv")

#convert descriptions to the tax_table
EC.tax<-data.frame(Path=EC[,1:2])
EC.tax<-tax_table(EC.tax)
EC.names<-EC[,1]
taxa_names(EC.tax)<-EC.names

#convert the rest of the table to the otu_table()
EC.otu<-EC[,3:dim(EC)[2]]
EC.otu<-round(EC.otu,0)
EC.otu<-otu_table(EC.otu,taxa_are_rows = T)
taxa_names(EC.otu)<-EC.names

#merge these objects with original metadata and save
HIRUNDO_picrust_path<-merge_phyloseq(EC.otu,EC.tax,sample_data(HIRUNDO_16S.taxo.phylo))
save(HIRUNDO_picrust_path,file = "/media/kreising/DATA/data/VLASTOVKY/PHYLOSEQ_rerun/HIRUNDO_picrust_path.R")

