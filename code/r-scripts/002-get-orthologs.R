# BiocManager::install(c("vitkl/orthologsBioMART"), dependencies=T)
library(tidyverse)
library(biomaRt)
library(orthologsBioMART)
source("code/r-scripts/get_ensembl_data_feature_and_structure.R")
ensembl.df.date_time<-"20250319_15_19_51"
coding.or.noncoding<-"coding"
ensembl.df.file_name<-paste0(paste(ensembl.df.date_time,"mouse-bm-all",coding.or.noncoding,sep = "-"),".tsv")
rtables.dir.path<-"output/rtables"
rdafiles.dir.path<-"output/rdafiles"
ensembl.df <- read.table(file.path(rtables.dir.path,ensembl.df.file_name), 
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
all.orthologs.dfs<-list()
for(organism.name in names(downloaded.species.dataset_names)){
  organism.dataset.name<-unname(downloaded.species.dataset_names[organism.name])
  organism.dataset.name.short<-gsub("_gene_ensembl","",organism.dataset.name)
  
  # ensembl.organism<-useEnsembl(biomart="ensembl",
  #                              version = 112,
  #                              dataset = organism.dataset.name)
  # marts.from.to<-new.env()
  # marts.from.to$ensembl_from<-ensembl.mmusculus
  # marts.from.to$ensembl_to<-ensembl.organism
  organism.df<-findOrthologs(datasets_FROM_TO =#marts.from.to,
                             loadBIOMARTdatasets(from = "mmusculus_gene_ensembl",
                                                 to = paste(organism.dataset.name.short,"gene_ensembl",sep="_")),
                             from_filters = "ensembl_gene_id", from_values = ensembl.gene_ids,
                             to_attributes = my.attr.feature[!my.attr.feature%in%c('entrezgene_id')],
                             to_homolog_attribute = paste(organism.dataset.name.short,"homolog_ensembl_gene",sep = "_"),
                             from_gene_id_name = "mmusculus_gene_id",
                             to_gene_id_name = paste(organism.dataset.name.short,"gene_id",sep="_"))
  
  join_df_vector<-c("ensembl_gene_id", 
                    "ensembl_gene_id_version", 
                    "ensembl_transcript_id", 
                    "ensembl_transcript_id_version")
  orthologs.df<-organism.df[c("mmusculus_gene_id",paste(organism.dataset.name.short,"gene_id",sep="_"))]
  orthologs.df<-orthologs.df%>%distinct(.keep_all = T)
  all.orthologs.dfs[[organism.dataset.name.short]]<-orthologs.df
  organism.df<-organism.df%>%
    rename(ensembl_gene_id=paste(organism.dataset.name.short,"gene_id",sep="_"))
  organism.df.all<-organism.df
  # organism.df.all<-get_ensembl_data_feature_and_structure(feature_attributes_vec = my.attr.feature[!my.attr.feature%in%c('entrezgene_id')],
  #                                                         structure_attributes_vec = my.attr.structure,
  #                                                         feature_values_vec =unique(organism.df$ensembl_gene_id),
  #                                                         filter_feature_by = "ensembl_gene_id",
  #                                                         mart_obj = ensembl.organism,filter_structure_by = "ensembl_gene_id",
  #                                                         join_dfs_by = join_df_vector)
  all.organism.dfs[[organism.dataset.name.short]]<-organism.df.all
  rtable.path<-file.path(rtables.dir.path,
                         paste(paste(format(Sys.time(),format="%Y%m%d"),
                                     format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                               organism.dataset.name.short,"bm-all",coding.or.noncoding,paste0(".tsv"),sep = "-"))
  rtable.path<-gsub("-\\.tsv",".tsv",rtable.path)
  print(paste(organism.dataset.name.short, "table was saved as",rtable.path))
  write.table(organism.df.all,
              file=rtable.path,
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = T)
  rm(organism.df.all)
}
all.organism.dfs[["mmusculus"]]<-ensembl.df[,c("gene_name",intersect(colnames(ensembl.df),colnames(all.orthologs[[1]])))]
saveRDS(all.organism.dfs,file = file.path(rdafiles.dir.path,
                                          paste(paste(format(Sys.time(),format="%Y%m%d"),
                                                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                                      "all-ensembl-orthologs.rds",sep="-")))
ortholog.matches<-all.orthologs.dfs %>% purrr::reduce(full_join, by = "mmusculus_gene_id")
write.table(ortholog.matches,file=file.path(rtables.dir.path,
                                            paste(paste(format(Sys.time(),format="%Y%m%d"),
                                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                                  "ensembl-ortholog-matches.tsv",sep="-")),
            row.names = F,quote = T,sep = "\t")
