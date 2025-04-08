library(tidyverse)
library(biomaRt)
library(httr)
library(jsonlite)
source("code/r-scripts/check_interpro_hierarchy.R")

#' ---
#' title: "Obtaining Ensembl IDs and coordinates of genes associated with a 
#' pathway based on its Gene ontology (GO) ID"
#' author: "Amir Rakhimov"
#' date: "`r Sys.Date()`"
#' ---

#' ## Intro
#' This script identifies Ensembl IDs and gene coordinates of genes associated 
#' with the canonical NF-kappaB signaling pathway (GO ID: 0007249). The gene names and coordinates
#' are specific to mouse (*Mus musculus*).
#' 
#' Input: table of MGI IDs (Mouse Genome Informatics IDs) of genes associated
#' with the NF-kappaB signaling pathway (data/metadata/20250313_go_0007249-select-mus_musculus.txt)
#' NF-kappaB signaling pathway genes are collected from Gene Ontology:
#' 
#' * regulates_closure: GO:0007249	
#' * taxon_subset_closure_label: Mus musculus
#' 
#' I obtained it from the web version 
#' (https://amigo.geneontology.org/amigo/search/bioentity?q=*%3A*&fq=regulates_closure:%22GO%3A0007249%22&sfq=document_category:%22bioentity%22)
#' and selected the following columns:
#' 
#' 1. Gene/product (bioentity). For example, MGI:MGI:2142330	
#' 2. Gene/product (bioentity_label). Ppm1n
#' 3. Gene/product name (bioentity_name). Probable protein phosphatase 1N
#' 4. Type (type). protein
#' 5. PANTHER family (panther_family). PANTHER:PTHR47992
#' 
#' Output: tab-separated file with each gene's MGI ID, MGI symbol,
#' Entrez ID, Ensembl ID, chromosome name, start position, end position,
#' strand, +/-10kb coordinates.
#' 
#' Result: 287 coding genes (output/rtables/20250319_15_03_16-mouse-gene-locations-noncoding.tsv)
#' and 7 non-coding genes (output/rtables/20250319_15_03_16-mouse-gene-locations-noncoding.tsv)
#' 
#' Workflow:
#' 1. Prepare the Ensembl dataset
#' 2. Find genes and transcripts with Biomart based on MGI IDs
#' 3. Find genes and transcripts with Biomart based on other IDs like UniProt gene
#' 4. Join these tables
#' 5. Find detailed structure information (exon coordinates) for the joined table
#' 6. Separate into coding and non-coding genes
#' 7. Make a small table with gene coordinates
#' 
#' 
#' 
#' 
#' 
#' ```{r, setup, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/nfkappab_evo')
#' ```
#' ## Import the table from GO database and extract gene IDs.
metadata.dir.path<-"./data/metadata" # directory with metadata
rtables.dir.path<-"./output/rtables"
# Path to the table with MGI IDs
go_genes.table<-read.table(file.path(metadata.dir.path,"20250313_go_0007249-select-mus_musculus.txt"),
                         header = F,sep = "\t")
head(go_genes.table)

# exon_chrom_start Exon region start (bp) (can be 5'UTR start)
# Exons include UTRs.
# cdna_coding_start: don't include UTR. Coordinates are wrt the transcript w/ UTR.
# cds_start: don't include UTRs. Coordinates are wrt the processed transcript w/o UTR
# genomic_coding_start: exon_chrom_start+cDNA coding start

# exon_chrom_end Exon region end (bp)
# rank Exon rank in transcript (exon order)
# transcript_length Transcript length (including UTRs and CDS)

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
my.attr.feature.mgi<-c("mgi_id","mgi_symbol",my.attr.feature)
my.attr.feature.uniprot_gn<-c("uniprot_gn_id","mgi_symbol","ensembl_gene_id")
my.attr.feature.uniprot_iso<-c("mgi_symbol","uniprot_isoform",
                                    "ensembl_gene_id")

#' Rename the columns according to their contents: MGI ID, gene product,
#' gene symbol, source (database; almost all are from MGI), type (protein,
#' gene product), and PANTHER family ID.
colnames(go_genes.table)<-c("mgi_id","gene_product","gene_symbol","source","type",
                      "panther_id")
# Substitute "MGI:MGI:" with "MGI:" for subsequent steps.
go_genes.table$mgi_id<-gsub("MGI:MGI:","MGI:",go_genes.table$mgi_id)
go_genes.table$mgi_id<-gsub("PR:","",go_genes.table$mgi_id)

#' Separate the data into IDs that come from MGI database and IDs from
#' other databases.
#' MGI IDs:
mgi.genes.table<-subset(go_genes.table,source=="MGI")
# Extract the IDs as a vector
mgi.genes<-mgi.genes.table$mgi_id
#' Non-MGI IDs:
non_mgi.genes.table<-subset(go_genes.table,source!="MGI")
# Extract the IDs as a vector
non_mgi.genes<-non_mgi.genes.table$mgi_id

#' ## Prepare the Ensembl dataset to search with biomaRt
#' I used the mouse dataset.
ensembl.ds<- useEnsembl(biomart="ensembl", 
                     version = 112,
                     dataset = "mmusculus_gene_ensembl")
#' If you don't know the name of the dataset, use the `searchDatasets()` function
# searchDatasets(mart = ensembl, pattern = "musculus")

#' The `getBM()` function has three arguments that need to be introduced: 
#' `filters`, `values`, and `attributes.`
# filters = listFilters(ensembl.ds)
# head(filters)

# attributes = listAttributes(ensembl.ds)
# attributes[1:5,]

listAttributes(mart=ensembl.ds)%>%head
#' ## Search MGI ID isoforms.
#' Find Entrez IDs, Ensembl IDs, and coordinates of genes from the `mgi.genes` vector.
#' Filter based on `mgi_id`. The `mart` object is the `ensembl.ds` dataset that 
#' we got earlier.
bm.feature.mgi<-getBM(
  attributes=my.attr.feature.mgi,
  filters = "mgi_id",
  values=mgi.genes,
  mart=ensembl.ds)
head(bm.feature.mgi)

#' ## Search non-MGI ID isoforms.
#' Find information for genes from the `non_mgi.genes` vector.
#' Here, filter based on `uniprot_isoform` because there's one isoform from the
#' UniProt database. 
#' The `mart` object is the `ensembl.ds` dataset that we got earlier.
#' First, get the Ensembl gene id because the isoform might point to a 
#' non-canonical transcript (and we want canonical ones). 
#' Then, find the other attributes based on ensembl gene ids.
bm.feature.uniprot_iso<-getBM(
  attributes=my.attr.feature.uniprot_iso,
  filters = "uniprot_isoform",
  values=non_mgi.genes,
  mart=ensembl.ds
)
head(bm.feature.uniprot_iso)
# Here, we get the other attributes based on ensembl gene id.
bm.feature.uniprot_iso<-getBM(
  attributes=my.attr.feature,
  filters = "ensembl_gene_id",
  values=bm.feature.uniprot_iso$ensembl_gene_id,
  mart=ensembl.ds
)
head(bm.feature.uniprot_iso)

# Convert the chromosome name to character for joining in subsequent steps.
bm.feature.uniprot_iso$chromosome_name<-
  as.character(bm.feature.uniprot_iso$chromosome_name)

#' Find information for genes from the `non_mgi.genes` vector.
#' Here, filter based on `uniprot_gn_id` because we have a gene from the UniProt 
#' database.
bm.feature.uniprot_gn<-getBM(
  attributes=my.attr.feature.uniprot_gn,
  filters = "uniprot_gn_id",
  values=non_mgi.genes,
  mart=ensembl.ds
)
# Here, we get the other attributes based on ensembl gene id.
bm.feature.uniprot_gn<-getBM(
  attributes=my.attr.feature,
  filters = "ensembl_gene_id",
  values=bm.feature.uniprot_gn$ensembl_gene_id,
  mart=ensembl.ds
)
head(bm.feature.uniprot_gn)

# Convert the chromosome name to character for joining in subsequent steps.
bm.feature.uniprot_gn$chromosome_name<-
  as.character(bm.feature.uniprot_gn$chromosome_name)
head(bm.feature.uniprot_gn)

#' ## Join feature tables by `ensembl_gene_id.`
bm.feature.all<-bind_rows(bm.feature.mgi,
                            bm.feature.uniprot_gn,
                            bm.feature.uniprot_iso)%>%
  dplyr::distinct(ensembl_transcript_id_version,.keep_all = T)%>%
  dplyr::select(all_of(my.attr.feature.mgi))
head(bm.feature.all)

#' Get unique stable ensembl gene_ids that correspond to the combined feature dataset.
#' We will update it later because we will add genes from the structure data.
bm.feature.all.ensembl_gene_id<-bm.feature.all%>%
  distinct(ensembl_gene_id)%>%
  arrange%>%
  pull

#' Find attributes from the structure page by ensembl gene_id because
#' MGI ID doesn't exist in the structure page
bm.structure.all<-getBM(
  attributes=my.attr.structure,
  filters = "ensembl_gene_id",
  values=bm.feature.all.ensembl_gene_id,
  mart=ensembl.ds)
head(bm.structure.all)


#' Bind two tables: keep only canonical ensembl transcripts
bm.feature.structure.all<-bm.feature.all%>%
  full_join(bm.structure.all,
            by = join_by(ensembl_gene_id, 
                         ensembl_gene_id_version, 
                         ensembl_transcript_id, 
                         ensembl_transcript_id_version),
            relationship = "many-to-many")%>%
  drop_na(transcript_is_canonical)%>%
  rename("gene_name"="external_gene_name")

# The mgi_id are the same, but mgi_symbol are different
setdiff(bm.feature.structure.all$mgi_id,mgi.genes.table$mgi_id)
setdiff(mgi.genes.table$mgi_id,bm.feature.structure.all$mgi_id)
setdiff(bm.feature.structure.all$mgi_symbol,mgi.genes.table$gene_symbol)
setdiff(mgi.genes.table$gene_symbol,bm.feature.structure.all$mgi_symbol)
#' Remove duplicated IDs (Trim30b) 
bm.feature.structure.all<-bm.feature.structure.all%>%
  group_by(ensembl_transcript_id_version,mgi_id,exon_chrom_start)%>%
  distinct(ensembl_transcript_id_version,mgi_id,exon_chrom_start,.keep_all = TRUE)%>%
  ungroup()
#' Verify that we didn't lose any MGI IDs or Ensembl genes
bm.feature.structure.all.mgi_id.uniq<-unique(bm.feature.structure.all$mgi_id)
# Remove NA mgi_id for verification
bm.feature.structure.all.mgi_id.uniq<-bm.feature.structure.all.mgi_id.uniq[!is.na(bm.feature.structure.all.mgi_id.uniq)]
# Same for mgi.genes.table
mgi.genes.table.mgi_id.uniq<-unique(mgi.genes.table$mgi_id)
# We want the joined biomart table to have same amount of mgi id
stopifnot(length(bm.feature.structure.all.mgi_id.uniq)==length(mgi.genes.table.mgi_id.uniq))
# But it's possible that the joined biomart has more ensembl gene ids because
# we are getting them from structure data (uniprot ids). So the length of two
# vectors doesn't have to be the same: iomart table length should be more 
# than or equal to the length of mgi id table
stopifnot(length(unique(bm.feature.structure.all$ensembl_gene_id_version))>=
            length(unique(bm.feature.mgi$ensembl_gene_id_version)))


# Separate coding and non-coding sequences (pseudogenes) into two tables
bm.feature.structure.all.coding<-bm.feature.structure.all%>%
  group_by(ensembl_gene_id_version)%>%
  filter(cds_length==max(cds_length))%>%
  ungroup()

bm.feature.structure.all.noncoding<-bm.feature.structure.all%>%
  group_by(ensembl_gene_id_version)%>%
  filter(is.na(cds_length))%>%
  ungroup()

# Get InterPro domains. We can get them only from coding sequences (no pseudogenes!)
# Extract the ensembl gene ids of coding sequences
bm.feature.structure.all.coding.ensembl_gene_id<-bm.feature.structure.all.coding%>%
  distinct(ensembl_gene_id)%>%
  arrange%>%
  pull
# Get the protein domain data from biomart 
bm.interpro.all<-getBM(
  attributes=my.attr.feature.interpro,
  filters = "ensembl_gene_id",
  values=bm.feature.structure.all.coding.ensembl_gene_id,
  mart=ensembl.ds)
# Extract unique InterPro IDs
bm.interpro.all.interpro_id<-bm.interpro.all%>%
  filter(transcript_is_canonical==1)%>%
  filter(interpro!="")%>%
  distinct(interpro)%>%
  pull
# Identify families and homologous superfamilies
interpro_ids.to.remove<-check_interpro_hierarchy(bm.interpro.all.interpro_id)
# interpro_ids.to.remove.tryagain<-check_interpro_hierarchy(interpro_ids.to.remove[[2]])
# Select canonical transcripts and remove families and superfamilies
bm.interpro.all.clean<-bm.interpro.all%>%
  filter(transcript_is_canonical==1)%>%
  filter(interpro!="")
  filter(!interpro%in%interpro_ids.to.remove[[1]])
# We must not lose ensembl gene ids: all elements of feature+structure biomart dataset
# should be found in the cleaned interpro dataset, and vice versa 
# TODO: some proteins don't have any domains, only families annotation.
stopifnot(length(setdiff(unique(bm.feature.structure.all.coding$ensembl_gene_id),
                         unique(bm.interpro.all.clean$ensembl_gene_id)))==0) # This one particularly
stopifnot(length(setdiff(unique(bm.interpro.all.clean$ensembl_gene_id),
                         unique(bm.feature.structure.all.coding$ensembl_gene_id)))==0)
head(bm.interpro.all.clean)
#' Make a small version with coordinates and the most important information.
#' Add +/-10kb coordinates.
get_gene_table_locations<-function(gene_table){
  gene_table_locations<-gene_table%>%
    dplyr::select(chromosome_name,gene_name,start_position,
                  end_position,strand)%>%
    distinct(gene_name,.keep_all = TRUE)%>%
    rename("gene_start_bp"="start_position",
           "gene_end_bp"="end_position")%>%
    mutate(strand_new=ifelse(strand=="1","+","-"),
           gene_start_bp_minus_10k=gene_start_bp-10000,
           gene_end_bp_plus_10k=gene_end_bp+10000)
  return(gene_table_locations)
}
gene.locations.coding<-get_gene_table_locations(bm.feature.structure.all.coding)
gene.locations.noncoding<-get_gene_table_locations(bm.feature.structure.all.noncoding)

#' Save the tables
for (rtable in c("bm.feature.structure.all.coding","bm.feature.structure.all.noncoding",
                 "gene.locations.coding","gene.locations.noncoding","bm.interpro.all.clean")){
  rtable.fname<-gsub("\\.","-",rtable)
  rtable.path<-file.path(rtables.dir.path,
                         paste(paste(format(Sys.time(),format="%Y%m%d"),
                                     format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                               "mouse",paste0(rtable.fname,".tsv"),sep = "-"))
  print(paste(rtable, "was saved as",rtable.path))
  write.table(get(rtable),
              file=rtable.path,
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = T)

}
save.image(file="./output/rdafiles/001-get-genes-go-to-ensembl-workspace.Rdata")

# Final table format in the paper by Xia et al. 
# (https://www.nature.com/articles/s41586-024-07095-8):
# Chromosome.scaffold.name	Gene.name	Gene.start..bp.	Gene.end..bp.	Strand	Strand_new	Gene.start..bp.-10k	Gene.end..bp.+10k
# 1	CSF1	109910242	109930992	1	+	109900242	109940992

# Final table format in this script:
# chromosome_name gene_name       gene_start_bp   gene_end_bp     strand  strand_new      gene_start_bp_minus_10k gene_end_bp_plus_10k
# 13      F2r     95738311        95754995        -1      -       95728311        95764995

load(file="./output/rdafiles/001-get-genes-go-to-ensembl-workspace.Rdata")
