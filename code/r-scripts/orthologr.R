# # Install Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# # Install package dependencies
# BiocManager::install(c(
#   "Biostrings",
#   "GenomicRanges",
#   "GenomicFeatures",
#   "Rsamtools",
#   "rtracklayer"
# ))
# 
# # install CRAN dependencies
# install.packages(c("doParallel", "foreach", "ape", "Rdpack", "benchmarkme", "devtools"))
# 
# # install BLAST dependency metablastr from GitHub
# devtools::install_github("drostlab/metablastr")
# 
# # install DIAMOND dependency rdiamond from GitHub
# devtools::install_github("drostlab/rdiamond")
# 
# # install orthologr from GitHub
# devtools::install_github("drostlab/orthologr")
# # install metablastr for tblastn
# devtools::install_github("drostlab/metablastr", build_vignettes = TRUE, 
#                          dependencies = TRUE)

library(tidyverse)
library(orthologr)
library(biomartr)
library(R.utils)
# library(metablastr)
# library(readxl)
metadata.dir.path<-"./data/metadata" # directory with metadata
rtables.dir.path<-"./output/rtables"
rdafiles.dir.path<-"./output/rdafiles"
proteomes.dir.path<-"./data/proteomes"
gff.dir.path<-"./data/gff"
cds.dir.path<-"./data/cds"
protein_homologs.dir.path<-"./data/protein_homologs"
genomes.dir.path<-"./data/genomes"
downloaded.organims.ensembl<-c("Homo sapiens",
                       "Mus musculus",
                       "Heterocephalus glaber",
                       "Rattus norvegicus",
                       "Loxodonta africana",
                       "Equus caballus",
                       "Bos taurus",
                       "Rhinolophus ferrumequinum",
                       "Myotis lucifugus",
                       "Cavia porcellus",
                       "Sus scrofa",
                       "Canis lupus familiaris",
                       "Pan troglodytes",
                       "Tursiops truncatus",
                       "Dasypus novemcinctus",
                       "Monodelphis domestica",
                       "Castor canadensis",
                       "Sarcophilus harrisii")
downloaded.organims.refseq<-c("Fukomys damarensis",
                              "Myotis davidii",
                              "Pteropus alecto")
# download the proteomes of three different species at the same time
#### Database: ENSEMBL
# file_paths.ensembl <- getProteomeSet(db = "ensembl", 
#                              organisms = downloaded.organims.ensembl,
#                              path = proteomes.dir.path)
#                                            
# file_paths.refseq <- getProteomeSet(db = "refseq", organisms = downloaded.organims.refseq,
#                                     path = proteomes.dir.path)

# Download genomes
# TODO: getCDSSet
genomes.file_path<-biomartr::getCDSSet(db="ensembl",
                                       organisms = downloaded.organims.ensembl,
                                         path = cds.dir.path)
# genome.file_path.dmr<-biomartr::getGenome(db="refseq",
#                                          organism="Fukomys damarensis",
#                                          reference = FALSE,
#                                          path = genomes.dir.path)
# gff.file_path.dmr<-biomartr::getGFF(db="refseq",
#                                          organism="Fukomys damarensis",
#                                          reference = FALSE,
#                                          path = gff.dir.path)
R.utils::gunzip(file_path.dmr.new,remove=FALSE)

# Get bowhead whale proteome
# whale.proteome.url<-"http://www.bowhead-whale.org/static/bowhead_whale_proteins.zip"
# whale.proteome.destfile<-file.path(proteomes.dir.path,"bowhead_whale_proteins.zip")
# download.file(url=whale.proteome.url, destfile = whale.proteome.destfile)
# unzip(whale.proteome.destfile,exdir = proteomes.dir.path)
# Get bowhead whale homologs 
# whale.homologs.url<-"http://www.bowhead-whale.org/static/homologs.zip"
# whale.homologs.destfile<-file.path(protein_homologs.dir.path,"bowhead_whale_homologs.zip")
# download.file(url=whale.homologs.url, destfile = whale.homologs.destfile)
# unzip(whale.homologs.destfile,exdir = protein_homologs.dir.path)

# Import proteomes and store them in a list of dataframes
proteome.list<-list.files(proteomes.dir.path)
proteome.list<-proteome.list[grepl("\\.fa",proteome.list)]
proteomes<-list()
for (proteome_file_name in proteome.list){
  proteome.obj<-biomartr::read_proteome(file = file.path(proteomes.dir.path,proteome_file_name))
  proteome_file_name.tidy<-gsub("\\.faa$|\\.fa$|\\.fasta$|\\.pep\\.all|_protein_refseq|_proteins","",proteome_file_name)
  proteomes[[proteome_file_name.tidy]]<-proteome.obj
  rm (proteome.obj)
}

proteome.names<-c("btaurus"="Bos_taurus.ARS-UCD1.3",
                  "bowhead_whale"="bowhead_whale",
                  "Castor_canadensis.C.can_genome_v1.0"="Castor_canadensis.C.can_genome_v1.0",
                  "cporcellus"="Cavia_porcellus.Cavpor3.0",
                  "ecaballus"="Equus_caballus.EquCab3.0",
                  "Fukomys_damarensis"="Fukomys_damarensis",
                  "hgfemale"="Heterocephalus_glaber_female.Naked_mole-rat_maternal",
                  "hsapiens"="Homo_sapiens.GRCh38",
                  "lafricana"="Loxodonta_africana.loxAfr3",
                  "mmusculus"="Mus_musculus.GRCm39",
                  "Myotis_davidii"="Myotis_davidii",
                  "mlucifugus"="Myotis_lucifugus.Myoluc2.0",
                  "Pteropus_alecto"="Pteropus_alecto",
                  "rnorvegicus"="Rattus_norvegicus.mRatBN7.2",
                  "rferrumequinum"="Rhinolophus_ferrumequinum.mRhiFer1_v1.p")
# Read the table with whale orthologs
# whale.orthologs<-read_xlsx(file.path(protein_homologs.dir.path,"bowhead_whale_homologs.xlsx"))
# colnames(whale.orthologs)<-tolower(colnames(whale.orthologs))
# colnames(whale.orthologs)[grepl("homologs in other species",colnames(whale.orthologs))]<-"homologs"
# write.table(whale.orthologs,file=file.path(rtables.dir.path,"whale-orthologs.tsv"),
#             sep="\t",row.names = F)

whale.orthologs<-read.table(file=file.path(rtables.dir.path,"whale-orthologs.tsv"),
            sep="\t",header = T)


# Get gene names for whale
# Read the rds file with ortholog data to get gene names
all.orthologs.ensembl.data.list<-readRDS(file.path(rdafiles.dir.path,"20250324_14_12_55-all-ensembl-orthologs.rds"))
# Extract whale contig names based on gene names of their orthologs that we need.
# Use the ensembl data (gene ids, gene names, transcript ids, etc) and select
# only gene names and mouse ensembl ids (the common column between all datasets).
# Store in a list.
all.ortholog.gene_names<-list()
# Get canonical transcript ids for blast
all.ortholog.transcript_ids<-list()
for(df.name in names(all.orthologs.ensembl.data.list)){
  ortholog.df<-all.orthologs.ensembl.data.list[[df.name]]
  if(df.name=="mmusculus"){
    ortholog.df<-ortholog.df%>%
      mutate(mmusculus_gene_id=ensembl_gene_id)
  }
  colnames(ortholog.df)[grepl("gene_name",
                              colnames(ortholog.df))]<-paste(df.name,"gene_name",sep="_")
  colnames(ortholog.df)[grepl("ensembl_transcript_id",
                              colnames(ortholog.df))]<-paste(df.name,"ensembl_transcript_id",sep="_")
  
  ortholog.df.gene_names<-ortholog.df[,c("mmusculus_gene_id",
                              paste(df.name,"gene_name",sep="_"))]%>%
    distinct(.keep_all = T)
  ortholog.df.transcript_ids<-ortholog.df%>%
    drop_na(transcript_is_canonical)%>%
    dplyr::select(all_of(c("mmusculus_gene_id",paste(df.name,"ensembl_transcript_id",sep="_"))))%>%
    distinct(.keep_all = T)
  all.ortholog.gene_names[[df.name]]<-ortholog.df.gene_names
  all.ortholog.transcript_ids[[df.name]]<-ortholog.df.transcript_ids
  rm(ortholog.df.gene_names)
  rm(ortholog.df.transcript_ids)
}
# Merge lists into one table by the common column
all.ortholog.gene_names<-all.ortholog.gene_names %>% 
  purrr::reduce(full_join, 
                by = "mmusculus_gene_id")
all.ortholog.transcript_ids<-all.ortholog.transcript_ids %>% 
  purrr::reduce(full_join, 
                by = "mmusculus_gene_id")

# select only some species (no bats)
selected.species<-c("mmusculus","hsapiens",
                    "hgfemale","rnorvegicus",
                    "ecaballus","btaurus","cporcellus")
search.gene_names<-all.ortholog.gene_names[,paste(selected.species,"gene_name",sep="_")]%>%
  unlist()%>%
  unique()
search.transcript_ids<-all.ortholog.transcript_ids[,paste(selected.species,"ensembl_transcript_id",sep="_")]%>%
  unlist()%>%
  unique()

# Remove NA
search.gene_names<-search.gene_names[!is.na(search.gene_names)]
search.transcript_ids<-search.transcript_ids[!is.na(search.transcript_ids)]
# Remove empty character string ("")
search.gene_names<-search.gene_names[nzchar(search.gene_names)] 
search.transcript_ids<-search.transcript_ids[nzchar(search.transcript_ids)] 
# Grep to retain only contigs we need for whale
selected.contigs.whale <- filter(whale.orthologs, grepl(paste(search.gene_names, collapse='|'), 
                                              homologs))




# Extract the sequences that we need for BLAST ####
# Import ensembl data for mouse
# bm.all.coding.mouse<-read.table(file.path(rtables.dir.path,"20250319_15_19_51-mouse-bm-all-coding.tsv"),
#                                 header = T,sep = "\t")
# head(bm.all.coding.mouse)
# Import table with ortholog matches
# ortholog.matches<-read.table(file.path(rtables.dir.path,"20250324_14_12_56-ensembl-ortholog-matches.tsv"),
#                              header = T,sep = "\t")
all.transcript_ids.selected_cols<-all.ortholog.transcript_ids[,paste(selected.species,"ensembl_transcript_id",sep="_")]
all.selected.transcript_seqs<-Biostrings::AAStringSetList()
# Extract sequences whose ids are found in the ortholog dataframe
for(species.name in colnames(all.transcript_ids.selected_cols)){
  species.name.tidy<-gsub("_ensembl_transcript_id","",species.name)
  selected.proteome<-proteome.names[species.name.tidy]
  species.transcript_ids<-all.transcript_ids.selected_cols%>%
    dplyr::select(all_of(species.name))%>%
    distinct()%>%
    pull
  species.transcript_ids<-species.transcript_ids[!is.na(species.transcript_ids)]
  species.transcript_ids<-species.transcript_ids[nzchar(species.transcript_ids)] 
  
  species.transcript_seqs<-proteomes[[selected.proteome]][grepl(paste(species.transcript_ids,collapse="|"),
                                                      names(proteomes[[selected.proteome]]))]
  all.selected.transcript_seqs[[species.name.tidy]]<-species.transcript_seqs
  seqs.path<- file.path("./output/fasta",paste(species.name.tidy,
                                             "selected_aa.fasta",sep="_"))
  Biostrings::writeXStringSet(species.transcript_seqs,filepath =seqs.path)
  print(paste(species.name.tidy, "sequences were saved as",seqs.path))
  rm(species.name)
  rm(species.name.tidy)
  rm(species.transcript_seqs)
  rm(species.transcript_ids)
}
Biostrings::writeXStringSet(unlist(all.selected.transcript_seqs),
                            filepath =file.path("./output/fasta/all_selected_aa.fasta"))
# save.image("//wsl.localhost/Ubuntu-22.04/home/rakhimov/projects/nfkappab_evo/output/rdafiles/orthologr-workspace.RData")
load("//wsl.localhost/Ubuntu-22.04/home/rakhimov/projects/nfkappab_evo/output/rdafiles/orthologr-workspace.RData")

# BLAST
tblastn.results<-
  metablastr::blast_protein_to_nucleotide(query = "./output/fasta/all_selected_aa.fasta",
                                          subject =file.path(genomes.dir.path,"Fukomys_damarensis_genomic_refseq.fna"),
                                          task = "tblastn" ,
                                          out.format = "csv",
                                          cores = 4)
