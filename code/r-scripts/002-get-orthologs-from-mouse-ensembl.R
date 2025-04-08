# BiocManager::install(c("vitkl/orthologsBioMART"), dependencies=T)
library(tidyverse)
library(biomaRt)
library(orthologsBioMART)
library(httr)
library(jsonlite)
source("code/r-scripts/get_ensembl_data_feature_and_structure.R")
source("code/r-scripts/check_interpro_hierarchy.R")
ensembl.df.date_time<-"20250402_16_09_42"
# interpro.df.date_time<-"20250402_16_09_43" #filtered
interpro.df.date_time<-"20250402_20_44_33" # unfiltered
# 
coding.or.noncoding<-"coding"
ensembl.df.file_name<-paste0(paste(ensembl.df.date_time,
                                   "mouse-bm-feature-structure-all",
                                   coding.or.noncoding,sep = "-"),".tsv")
interpro.df.file_name<-paste0(paste(interpro.df.date_time,
                                   "mouse-bm-interpro-all-clean",sep = "-"),".tsv")
rtables.dir.path<-"output/rtables"
rdafiles.dir.path<-"output/rdafiles"
ensembl.df <- read.table(file.path(rtables.dir.path,ensembl.df.file_name), 
                         header=T,
                         sep="\t")
interpro.df <- read.table(file.path(rtables.dir.path,interpro.df.file_name), 
                         header=T,
                         sep="\t")
colnames(ensembl.df)[grep("^X[0-9]_utr",colnames(ensembl.df))]<-
  gsub("^X","",colnames(ensembl.df)[grep("^X[0-9]_utr",colnames(ensembl.df))])

ensembl.gene_ids<-unique(ensembl.df$ensembl_gene_id)

#' ## Prepare the Ensembl dataset to search with biomaRt
#' I used the mouse dataset.
downloaded.organims.ensembl<-c("Homo sapiens",
  "Heterocephalus glaber female",
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
  "Sarcophilus harrisii")
abbreviate_species <- function(species) {
  species<-tolower(species)
  parts <- strsplit(species, " ")[[1]]  # Split by space
  abbr <- paste0(substr(parts[-length(parts)], 1, 1), collapse = "")  # First letter of all but the last
  paste0(abbr, parts[length(parts)],"_gene_ensembl")  # Append the last word
}
downloaded.species.dataset_names<-sapply(downloaded.organims.ensembl,abbreviate_species)

my.attr.feature<-c("entrezgene_id",
                   "chromosome_name",
                   "start_position","end_position","strand",
                   "ensembl_gene_id","ensembl_gene_id_version",
                   "ensembl_transcript_id","ensembl_transcript_id_version",
                   "external_gene_name","transcript_start",
                   "transcript_end","transcription_start_site",
                   "transcript_length","transcript_count","transcript_is_canonical")
# InterPro start and end coordinates refer to amino acid sequence
my.attr.feature.interpro<-c("ensembl_gene_id","ensembl_gene_id_version",
                            "ensembl_transcript_id","ensembl_transcript_id_version",
                            "transcript_is_canonical","external_gene_name",
                            "interpro","interpro_short_description",
                            "interpro_description","interpro_start","interpro_end")
my.attr.structure<-c("ensembl_gene_id","ensembl_gene_id_version",
                     "ensembl_transcript_id","ensembl_transcript_id_version",
                     "exon_chrom_start","exon_chrom_end",
                     "rank","cdna_coding_start",
                     "cdna_coding_end",
                     "cds_start","cds_end",
                     "genomic_coding_start","genomic_coding_end",
                     "5_utr_start","5_utr_end",
                     "3_utr_start","3_utr_end","cds_length")
all.organism.dfs<-list()
all.ortholog.mappings<-list()
all.domains.dfs<-list()
all.failed.interpro_ids<-list()
for(organism.name in names(downloaded.species.dataset_names)){
  organism.dataset.name<-unname(downloaded.species.dataset_names[organism.name])
  organism.dataset.name.short<-gsub("_gene_ensembl","",organism.dataset.name)
  
  # ensembl.organism<-useEnsembl(biomart="ensembl",
  #                              version = 112,
  #                              dataset = organism.dataset.name)
  # marts.from.to<-new.env()
  # marts.from.to$ensembl_from<-ensembl.mmusculus
  # marts.from.to$ensembl_to<-ensembl.organism
  # Get ortholog data (gene ids, coordinates, etc)
  organism.df<-findOrthologs(datasets_FROM_TO =#marts.from.to,
                             loadBIOMARTdatasets(from = "mmusculus_gene_ensembl",
                                                 to = paste(organism.dataset.name.short,"gene_ensembl",sep="_")),
                             from_filters = "ensembl_gene_id", from_values = ensembl.gene_ids,
                             to_attributes = my.attr.feature[!my.attr.feature%in%c('entrezgene_id')],
                             to_homolog_attribute = paste(organism.dataset.name.short,"homolog_ensembl_gene",sep = "_"),
                             from_gene_id_name = "mmusculus_gene_id",
                             to_gene_id_name = paste(organism.dataset.name.short,"gene_id",sep="_"))
  organism.df.interpro<-findOrthologs(datasets_FROM_TO =#marts.from.to,
                             loadBIOMARTdatasets(from = "mmusculus_gene_ensembl",
                                                 to = paste(organism.dataset.name.short,"gene_ensembl",sep="_")),
                             from_filters = "ensembl_gene_id", from_values = ensembl.gene_ids,
                             to_attributes = my.attr.feature.interpro,
                             to_homolog_attribute = paste(organism.dataset.name.short,"homolog_ensembl_gene",sep = "_"),
                             from_gene_id_name = "mmusculus_gene_id",
                             to_gene_id_name = paste(organism.dataset.name.short,"gene_id",sep="_"))
  join_df_vector<-c("ensembl_gene_id", 
                    "ensembl_gene_id_version", 
                    "ensembl_transcript_id", 
                    "ensembl_transcript_id_version")
  # Extract columns with ortholog ensembl ids (mapping between mus musculus and 
  # the other species). Example columns are mmusculus_gene_id and hsapiens_gene_id
  ortholog.mapping<-organism.df[c("mmusculus_gene_id",paste(organism.dataset.name.short,"gene_id",sep="_"))]
  ortholog.mapping<-ortholog.mapping%>%distinct(.keep_all = T)
  # Add the mapping between orthologs to the list of orthologs
  all.ortholog.mappings[[organism.dataset.name.short]]<-ortholog.mapping
  # Rename the column hsapiens_gene_id to generic ensembl_gene_id
  organism.df<-organism.df%>%
    dplyr::rename(ensembl_gene_id=paste(organism.dataset.name.short,"gene_id",sep="_"))
  organism.df.interpro<-organism.df.interpro%>%
    dplyr::rename(ensembl_gene_id=paste(organism.dataset.name.short,"gene_id",sep="_"))
  organism.df.feature.structure.all<-organism.df
  # organism.df.feature.structure.all<-get_ensembl_data_feature_and_structure(feature_attributes_vec = my.attr.feature[!my.attr.feature%in%c('entrezgene_id')],
  #                                                         structure_attributes_vec = my.attr.structure,
  #                                                         feature_values_vec =unique(organism.df$ensembl_gene_id),
  #                                                         filter_feature_by = "ensembl_gene_id",
  #                                                         mart_obj = ensembl.organism,filter_structure_by = "ensembl_gene_id",
  #                                                         join_dfs_by = join_df_vector)
  # Get InterPro domains
  # Take the ensembl data that mapped to mouse ensembl genes and contains Interpro IDs,
  # extract only canonical transcripts and keep unique Interpro IDs
  organism.df.interpro.interpro_id<-organism.df.interpro%>%
    filter(transcript_is_canonical==1)%>%
    filter(interpro!="")%>%
    distinct(interpro)%>%
    pull
  # TODO: Run the custom function to remove Interpro families and homologous superfamilies
  # interpro_ids.to.remove<-check_interpro_hierarchy(organism.df.interpro.interpro_id)
  # interpro_ids.to.remove.tryagain<-check_interpro_hierarchy(interpro_ids.to.remove[[2]])
  # Remove the IDs from the dataset
  organism.df.interpro.clean<-organism.df.interpro%>%
    filter(transcript_is_canonical==1)%>%
    filter(interpro!="")#%>% #TODO:double check
    # filter(!interpro%in%interpro_ids.to.remove[[1]])
  # Store the IDs that weren't obtained
  # all.failed.interpro_ids[[organism.dataset.name.short]]<-interpro_ids.to.remove[[2]]
  # stopifnot(length(setdiff(unique(organism.df.feature.structure.all$ensembl_gene_id),
  #                          unique(organism.df.interpro.clean$ensembl_gene_id)))==0)
  # stopifnot(length(setdiff(unique(organism.df.interpro.clean$ensembl_gene_id),
  #                          unique(organism.df.feature.structure.all$ensembl_gene_id)))==0)
  all.domains.dfs[[organism.dataset.name.short]]<-organism.df.interpro.clean
  # Add the ensembl dataset to the list
  all.organism.dfs[[organism.dataset.name.short]]<-organism.df.feature.structure.all
  
  # for (rtable in c("organism.df.feature.structure.all","organism.df.interpro.clean")){
  for (rtable in c("organism.df.interpro.clean")){
    rtable.fname<-gsub("\\.","-",rtable)
    rtable.fname<-gsub("organism-df","bm",rtable.fname)
    rtable.path<-file.path(rtables.dir.path,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 organism.dataset.name.short,rtable.fname,
                                 paste0(coding.or.noncoding,".tsv"),sep = "-"))
    print(paste(organism.dataset.name.short, "table was saved as",rtable.path))
    write.table(get(rtable),
                file=rtable.path,
                sep = "\t",
                quote = F,
                row.names = F,
                col.names = T)
    
  }
  rm(organism.df.feature.structure.all)
  rm(organism.df)
  rm(organism.df.interpro)
  rm(organism.df.interpro.clean)
  rm(interpro_ids.to.remove)
  rm(ortholog.mapping)
}
all.failed.interpro_ids
# Add the mouse dataset to the list. Column names are those that are present in 
# all datasets (15 columns).
all.organism.dfs[["mmusculus"]]<-ensembl.df[,c("gene_name",
                                               intersect(colnames(ensembl.df),
                                                         colnames(all.organism.dfs[[1]])))]
saveRDS(all.organism.dfs,file = file.path(rdafiles.dir.path,
                                          paste(paste(format(Sys.time(),format="%Y%m%d"),
                                                            format(Sys.time(),format = "%H_%M_%S"),
                                                      sep = "_"),
                                                      "all-ensembl-organism-dfs.rds",sep="-")))
ortholog.matches<-all.ortholog.mappings %>% purrr::reduce(full_join, by = "mmusculus_gene_id")
write.table(ortholog.matches,file=file.path(rtables.dir.path,
                                            paste(paste(format(Sys.time(),format="%Y%m%d"),
                                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                                  "ensembl-ortholog-mappings.tsv",sep="-")),
            row.names = F,quote = T,sep = "\t")
all.domains.dfs[["mmusculus"]]<-interpro.df[,c(intersect(colnames(interpro.df),
                                                         colnames(all.domains.dfs[[1]])))]
saveRDS(all.domains.dfs,file = file.path(rdafiles.dir.path,
                                          paste(paste(format(Sys.time(),format="%Y%m%d"),
                                                      format(Sys.time(),format = "%H_%M_%S"),
                                                      sep = "_"),
                                                "all-ensembl-interpro_domains-dfs.rds",sep="-")))
save.image(file.path("./output/rdafiles",paste(paste(format(Sys.time(),format="%Y%m%d"),
                                            format(Sys.time(),format = "%H_%M_%S"),
                                            sep = "_"),
                                      "002-get-orthologs-from-mouse-ensembl-workspace.Rdata",sep = "-")))
workspace.date_time<-"20250402_16_28_30"
load(file.path("./output/rdafiles",workspace.date_time,
               "002-get-orthologs-from-mouse-ensembl-workspace.Rdata",sep = "-"))
