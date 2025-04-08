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
library(Biostrings)
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
all.orthologs.ensembl.date_time<-"20250402_21_01_45" #using unfiltered
all.orthologs.ensembl.date_time<- "20250402_16_24_18"
all.orthologs.interpro.date_time<-"20250402_21_01_53"  #using unfiltered
# all.orthologs.interpro.date_time<-"20250402_16_27_03"
all.orthologs.ensembl.filepath<-file.path(rdafiles.dir.path,
                                          paste(all.orthologs.ensembl.date_time,
                                                "all-ensembl-organism-dfs.rds",sep="-"))
all.orthologs.interpro.filepath<-file.path(rdafiles.dir.path,
                                           paste(all.orthologs.interpro.date_time,
                                                 "all-ensembl-interpro_domains-dfs.rds",sep = "-"))

downloaded.organims.ensembl<-c("Homo sapiens",
                       "Mus musculus",
                       "GCA_944319715.1", #NMR Female
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
# Download the proteomes of different species at the same time ####
#### Database: ENSEMBL
# file_paths.ensembl <- getProteomeSet(db = "ensembl", 
#                              organisms = downloaded.organims.ensembl,
#                              path = proteomes.dir.path)
#                                            
# file_paths.refseq <- getProteomeSet(db = "refseq", organisms = downloaded.organims.refseq,
#                                     path = proteomes.dir.path)

# Download genomes ####
# TODO: getCDSSet
# genomes.file_path<-biomartr::getCDSSet(db="ensembl",
#                                        organisms = downloaded.organims.ensembl,
#                                          path = cds.dir.path)
# genome.file_path.dmr<-biomartr::getGenome(db="refseq",
#                                          organism="Fukomys damarensis",
#                                          reference = FALSE,
#                                          path = genomes.dir.path)
# gff.file_path.dmr<-biomartr::getGFF(db="refseq",
#                                          organism="Fukomys damarensis",
#                                          reference = FALSE,
#                                          path = gff.dir.path)
# R.utils::gunzip(file_paths.ensembl,remove=FALSE)

# Get bowhead whale proteome ####
# whale.proteome.url<-"http://www.bowhead-whale.org/static/bowhead_whale_proteins.zip"
# whale.proteome.destfile<-file.path(proteomes.dir.path,"bowhead_whale_proteins.zip")
# download.file(url=whale.proteome.url, destfile = whale.proteome.destfile)
# unzip(whale.proteome.destfile,exdir = proteomes.dir.path)
# Get bowhead whale homologs 
# whale.homologs.url<-"http://www.bowhead-whale.org/static/homologs.zip"
# whale.homologs.destfile<-file.path(protein_homologs.dir.path,"bowhead_whale_homologs.zip")
# download.file(url=whale.homologs.url, destfile = whale.homologs.destfile)
# unzip(whale.homologs.destfile,exdir = protein_homologs.dir.path)

# Import proteomes and store them in a list of dataframes ####
proteome.list<-list.files(proteomes.dir.path)
proteome.list<-proteome.list[grepl("\\.fa",proteome.list)]
proteomes<-list()
for (proteome_file_name in proteome.list){
  proteome.obj<-biomartr::read_proteome(file = file.path(proteomes.dir.path,proteome_file_name))
  proteome_file_name.tidy<-gsub("\\.faa$|\\.fa$|\\.fasta$|\\.pep\\.all","",proteome_file_name)
  proteome_file_name.tidy<-gsub("_protein_refseq|_proteins","",proteome_file_name.tidy)
  proteome_file_name.tidy<-gsub("\\..*","",proteome_file_name.tidy)
  proteomes[[proteome_file_name.tidy]]<-proteome.obj
  rm (proteome.obj)
}

# Convert proteome names to ensembl format (first letter of all strings before 
# species + species). For example, Homo sapiens -> hsapiens
# Canis lupus familiaris -> clfamiliaris
abbreviate_species <- function(species) {
  species<-tolower(species)
  parts <- strsplit(species, "_")[[1]]  # Split by underscore
  abbr <- paste0(substr(parts[-length(parts)], 1, 1), collapse = "")  # First letter of all but the last
  paste0(abbr, parts[length(parts)])  # Append the last word
}
proteome.names<-sapply(names(proteomes),abbreviate_species)

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
all.orthologs.ensembl.data.list<-readRDS(all.orthologs.ensembl.filepath)
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
    # If it's mouse dataset, create an additional column called mmusculus_gene_id 
    # (same as ensembl_gene_id)
    ortholog.df<-ortholog.df%>%
      mutate(mmusculus_gene_id=ensembl_gene_id)
  }
  # Find column called gene_name and attach the name of the organism
  colnames(ortholog.df)[grepl("gene_name",
                              colnames(ortholog.df))]<-paste(df.name,"gene_name",sep="_")
  # Find column called ensembl_transcript_id and attach the name of the organism
  colnames(ortholog.df)[grepl("ensembl_transcript_id",
                              colnames(ortholog.df))]<-paste(df.name,"ensembl_transcript_id",sep="_")
  # This is a df with mapping between gene_id and gene_name of mouse and other organism
  ortholog.df.gene_names<-ortholog.df[,c("mmusculus_gene_id",
                              paste(df.name,"gene_name",sep="_"))]%>%
    distinct(.keep_all = T)
  # Keep canonical transcripts and retain only gene_id and transcript_id
  ortholog.df.transcript_ids<-ortholog.df%>%
    drop_na(transcript_is_canonical)%>%
    dplyr::select(all_of(c("mmusculus_gene_id",paste(df.name,"ensembl_transcript_id",sep="_"))))%>%
    distinct(.keep_all = T)
  # Save gene_names and transcript_ids in the lists
  all.ortholog.gene_names[[df.name]]<-ortholog.df.gene_names
  all.ortholog.transcript_ids[[df.name]]<-ortholog.df.transcript_ids
  rm(ortholog.df.gene_names)
  rm(ortholog.df.transcript_ids)
  rm(ortholog.df)
}
# Merge lists into one table by the common column
all.ortholog.gene_names<-all.ortholog.gene_names %>% 
  purrr::reduce(full_join, 
                by = "mmusculus_gene_id")
all.ortholog.transcript_ids<-all.ortholog.transcript_ids %>% 
  purrr::reduce(full_join, 
                by = "mmusculus_gene_id")

# select only some species (no bats)
selected.species.for_whale<-c("mmusculus","hsapiens",
                    "hgfemale","rnorvegicus",
                    "ecaballus","btaurus","cporcellus")
gene_names.for_whale<-all.ortholog.gene_names[,paste(selected.species.for_whale,"gene_name",sep="_")]%>%
  unlist()%>%
  unique()
transcript_ids.for_whale<-all.ortholog.transcript_ids[,paste(selected.species.for_whale,"ensembl_transcript_id",sep="_")]%>%
  unlist()%>%
  unique()

# Remove NA
gene_names.for_whale<-gene_names.for_whale[!is.na(gene_names.for_whale)]
transcript_ids.for_whale<-transcript_ids.for_whale[!is.na(transcript_ids.for_whale)]
# Remove empty character string ("")
gene_names.for_whale<-gene_names.for_whale[nzchar(gene_names.for_whale)] 
transcript_ids.for_whale<-transcript_ids.for_whale[nzchar(transcript_ids.for_whale)] 
# Grep to retain only contigs we need for whale
selected.contigs.for_whale <- filter(whale.orthologs, grepl(paste(gene_names.for_whale, collapse='|'), 
                                              homologs))




# Extract transcript sequences that we need for BLAST ####
# Import ensembl data for mouse
# bm.all.coding.mouse<-read.table(file.path(rtables.dir.path,"20250319_15_19_51-mouse-bm-all-coding.tsv"),
#                                 header = T,sep = "\t")
# head(bm.all.coding.mouse)
# Import table with ortholog matches
# ortholog.matches<-read.table(file.path(rtables.dir.path,"20250324_14_12_56-ensembl-ortholog-matches.tsv"),
#                              header = T,sep = "\t")
all.transcript_ids.selected_cols<-all.ortholog.transcript_ids[,paste(selected.species.for_whale,"ensembl_transcript_id",sep="_")]
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






# Extract InterPro protein domains for BLAST ####
all.orthologs.interpro.data.list<-readRDS(all.orthologs.interpro.filepath)
species.for_blast<-intersect(paste(proteome.names,"ensembl_transcript_id",sep="_"),
                             colnames(all.ortholog.transcript_ids))
all.transcript_ids.selected_cols<-all.ortholog.transcript_ids[,species.for_blast]
all.selected.transcript_seqs<-Biostrings::AAStringSetList()
all.selected.domain_seqs<-Biostrings::AAStringSetList()
all.selected.domain_seqs<-list()
# Extract sequences whose ids are found in the ortholog dataframe
for(species.name in colnames(all.transcript_ids.selected_cols)){
  # start from the name that corresponds to ensembl dataset (mmusculus or hsapiens)
  species.name.short<-gsub("_ensembl_transcript_id","",species.name)
  # get the name of the proteome that corresponds to selected species name
  # (which has ensembl format)
  selected.proteome.name<-names(proteome.names[which(proteome.names==species.name.short)])
  # Get the proteome to extract sequences
  selected.proteome<-proteomes[[selected.proteome.name]]
  # Get the dataframe with Interpro domain coordinates of the organism
  selected.domains.coord<-all.orthologs.interpro.data.list[[species.name.short]]
  print(paste("Number of domain sequences in",species.name.short,"data:", 
              nrow(selected.domains.coord)))
  
  # Get all unique transcript ids
  selected.transcript_ids<-all.transcript_ids.selected_cols%>%
    dplyr::select(all_of(species.name))%>%
    distinct()%>%
    pull
  # Remove NA
  selected.transcript_ids<-selected.transcript_ids[!is.na(selected.transcript_ids)]
  # Remove empty strings
  selected.transcript_ids<-selected.transcript_ids[nzchar(selected.transcript_ids)] 
  selected.transcript_ids<-sort(selected.transcript_ids)
  # Get transcript sequences
  selected.transcript_seqs<-selected.proteome[grepl(paste0(selected.transcript_ids,collapse="|"),
                               names(selected.proteome))]
  # Save full protein sequences for each species
  all.selected.transcript_seqs[[species.name.short]]<-selected.transcript_seqs
  
  # For extracting domain sequences:
  # First, create an empty AAStringset called selected.aa.set.
  # Then, iterate through the transcripts in the dataset of interpro domains 
  # (row by row) and add the full protein sequence of the transcript
  # into the selected.aa.set. We will use the selected.aa.set to extract domains. 
  # Therefore, selected.aa.set will have repetitions (one transcript
  # can have many domains).
  selected.aa.set <- Biostrings::AAStringSet()
  for (protein in selected.domains.coord$ensembl_transcript_id) {
    selected.aa.set <- c(selected.aa.set, selected.transcript_seqs[grepl(protein,names(selected.transcript_seqs))])
  }
  # This is a quick way to extract domains: provide a vector of 
  # sequences, a vector of start coordinates, and a vector of end coordinates
  selected.domains.seqs<-Biostrings::subseq(selected.aa.set,
                          start=selected.domains.coord$interpro_start,
                          end=selected.domains.coord$interpro_end)
  # Extract Interpro domain ids from sequence names, the transcript sequence names,
  # their Ensembl transcript IDs, and domain(!) sequences.
  # selected.domains.coord$seqs<-selected.domains.seqs
  # selected.domains.coord$seq_name<-names(selected.domains.seqs)
  selected.domains.seqs.df <- data.frame(
    interpro_id=selected.domains.coord$interpro,
    seq_name = names(selected.domains.seqs),
    ensembl_transcript_id = selected.domains.coord$ensembl_transcript_id,
    seqs=selected.domains.seqs,
    start_end=paste(selected.domains.coord$interpro_start,
                    selected.domains.coord$interpro_end,sep = "_"),
    stringsAsFactors = FALSE
  )
  selected.domains.seqs.df<-selected.domains.seqs.df%>%
    dplyr::distinct(seqs,interpro_id,.keep_all = T)
  
  
  # Split the dataframe by interproid. We create a list of dataframes where
  # each dataframe contains domain information and sequences for that interpro id.
  selected.domains.seqs.list<-split(selected.domains.seqs.df,selected.domains.seqs.df$interpro_id)
  
  # Set names to each sequence according to the transcript id, interpro id, 
  # and positions of domain
  # selected.domains.seqs.list <- lapply(selected.domains.seqs.list, function(sub_df) {
  #   setNames(sub_df$seqs, sub_df$ensembl_transcript_id)
  # })
  
  selected.domains.seqs.list <- lapply(selected.domains.seqs.list, function(sub_df) {
    # In case there are sequences with duplicated ensembl ids and interpro ids, 
    # we add a '-' and the number of sequence (1, 2, 3, etc.)
    unique_names <- make.unique(paste(sub_df$ensembl_transcript_id, 
                                      sub_df$interpro_id, 
                                      sub_df$start_end,
                                      sep = "_"),
                                sep = "-")
    setNames(sub_df$seqs, unique_names)
  })
  
  
  # Convert each list entry into AAstringset
  selected.domains.seqs.list <- lapply(selected.domains.seqs.list, Biostrings::AAStringSet)
  # Convert the whole list into AAStringsetlist
  selected.domains.seqs.list <- AAStringSetList(selected.domains.seqs.list)
  

  # Verify that all transcript domains are separated correctly: transcript
  # names should match and sequences should match
  all_correct <- all(sapply(names(selected.domains.seqs.list), function(g) {
    names_match <- identical(names(selected.domains.seqs.list[[g]]), 
                             selected.domains.seqs.df$ensembl_transcript_id[selected.domains.seqs.df$interpro_id == g])
    sequences_match <- setequal(selected.domains.seqs.list[[g]], 
                                selected.domains.seqs.df$seqs[selected.domains.seqs.df$interpro_id == g])  # Ensure it's a valid comparison
    
    return(names_match && sequences_match)  
  }))
  
  # stopifnot(all_correct) #TODO: verify

  print(paste("Number of unique domains in",species.name.short,"data:", 
              length(unique(selected.domains.coord$interpro))))
  print(paste("Number of domain sequences in",species.name.short,"data:", 
              length(unlist(selected.domains.seqs.list))))
  # all.selected.domain_seqs[[species.name.short]]<-selected.domains.seqs #remove
  # Save the domain sequences in the list of all domains. First, we group by 
  # organism; then, by domain
  all.selected.domain_seqs[[species.name.short]]<-selected.domains.seqs.list
  
  organism.fasta.dir<-file.path("./output/fasta",paste(species.name.short,
                                                       "interpro_domains",sep="_"))
  ifelse(!dir.exists(file.path(organism.fasta.dir)),
         dir.create(file.path(organism.fasta.dir)),
         "Directory Exists")
  
  organism.seqs.ungrouped.path<- file.path("./output/fasta",paste(species.name.short,
                                               "domains_aa.fasta",sep="_"))
  
  
  # Biostrings::writeXStringSet(selected.domains.seqs,filepath =organism.seqs.ungrouped.path)
  Biostrings::writeXStringSet(unlist(selected.domains.seqs.list,use.names = F),
                              filepath =organism.seqs.ungrouped.path)
  # Write each domain set into a separate file (by interpro id) in the
  # organism folder
  for (interpro in names(selected.domains.seqs.list)) {
    seqs.path<- file.path(organism.fasta.dir,
                          paste(species.name.short,interpro,"selected_aa.fasta",sep="_"))
    writeXStringSet(selected.domains.seqs.list[[interpro]], filepath = seqs.path)
  }
  
  print(paste(species.name.short, "sequences were saved as",organism.seqs.ungrouped.path))
  rm(species.name)
  rm(species.name.short)
  rm(selected.transcript_seqs)
  rm(selected.transcript_ids)
  rm(selected.domains.coord)
  rm(selected.domains.seqs.df)
  rm(selected.aa.set)
  rm(selected.domains.seqs.list)
  rm(selected.domains.seqs)
  rm(selected.proteome)
  rm(selected.proteome.name)
  rm(protein)
}

species_groups<-all.selected.domain_seqs
all_domains <- unique(unlist(lapply(species_groups, names)))

# Initialize the inverted structure
inverted_groups <- list()

for (domain in all_domains) {
  species_list <- list()
  for (species_name in names(species_groups)) {
    # Get the AAStringSetList for the current species
    species <- species_groups[[species_name]]
    if (domain %in% names(species)) {
      # Extract sequences for the domain in this species
      species_list[[species_name]] <- species[[domain]]
    }
  }
  # Convert to AAStringSetList and assign to the domain
  inverted_groups[[domain]] <- AAStringSetList(species_list)
}


domain.fasta.dir<-file.path("./output/fasta/interpro")

ifelse(!dir.exists(domain.fasta.dir),
       dir.create(domain.fasta.dir),
       "Directory Exists")


for (interpro in names(inverted_groups)) {
  seqs.path<- file.path(domain.fasta.dir,
                        paste(interpro,"fasta",sep="."))
  writeXStringSet(unlist(inverted_groups[[interpro]],use.names = T), filepath = seqs.path)
}
