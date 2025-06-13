library(tidyverse)
library(ggpubr)
library(readxl)
library(tidyverse)
xl.file<-read_xlsx(path="./data/luminometer_data/20250612_try3_Amir_1.xlsx")
convert_luc_xl_to_tibble<-function(xl_file){
  # For single luciferase data (only one Read column)
  # Input: Excel file (raw)
  # Output: Tibble with "Sample" and "Read" columns
  row.all.na<-apply(xl_file, 1, function(x){all(is.na(x))})# Find rows with all NA
  xl_file_filtered<-xl_file[!row.all.na,] # Remove rows with all NA
  xl_file_filtered[1,1]<-"Measurement" # Remove "Measurement" column (first col)
  names(xl_file_filtered)<-xl_file_filtered[1,] # Column names are "Measurement",
                              # "Sample", "Read", "Remarks"
  xl_file_filtered<-xl_file_filtered[-1,] # Remove the row with names
  xl_file_filtered<-xl_file_filtered%>%
    select(-Remarks,-Measurement)%>%
    filter(grepl("Sample",Sample))%>% # Select only "Sample" and "Read" columns
    mutate(Sample=gsub(" \\/ Repl 1", "", Sample))%>%
    mutate(Read=as.numeric(Read))
  
  return(xl_file_filtered)
}


bind_two_luc_tibbles<-function(luc_file_1,luc_file_2){
  # Input: two tibbles with "Sample" and "Read" columns
  # Output: combined tibble
  # Add real sample number
  luc_file_1<-luc_file_1%>%
    mutate(sample_num=grep("[0-9]+",Sample))
  luc_file_2<-luc_file_2%>%
    mutate(sample_num=grep("[0-9]+",Sample)+nrow(luc_file_1))
  # Combine datasets
  luc_files_combined<-rbind(luc_file_1,luc_file_2)%>%
    select(-Sample)%>%
    rename("Sample"=sample_num)%>%
    relocate(Sample)
  return(luc_files_combined)
}

#' ---
#' title: "Analysing luminometer data in R"
#' author: "Amir Rakhimov"
#' date: "3 Mar, 2025"
#' ---

#' ## Intro
#' We are using two cell lines: **6W3, TIG111**
#' ```{r, setup, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/nfkappab_evo')
#' ```
#'            LAR II      Stop&Glo
#'              1             2         3 
#' o o o Ctrl   x             y         x/y   <- NF-kappaB-dependent firefly activity
#' o o o LPS
#' o o o TNF
#'
#' Ctrl (x/y) o o o -> get Ctrl Average = Ctrl activation -> bar in the barplot (value =1)
#' LPS (x/y) o o o -> get LPS Average -> bar in the barplot
#' TNF (x/y) o o o -> get TNF Average -> bar in the barplot
#' 
#'  LPS (x/y)
#' --------------  = Fold activation of LPS
#' Ctrl activation
#'
#'  TNF (x/y)
#' --------------  = Fold activation of TNF
#' Ctrl activation
#'
#'
getwd()
treatment.labels<-c("control"="Control",
                    "LPS"="LPS",
                    "TNF"="TNF-alpha")
# assay.name<-"20250122_H17_6W4"
# luc.df$cell_line<-gsub("H17","NMR_H17",luc.df$cell_line)
# luc.df$cell_line<-gsub("6W4","Ms_6W4",luc.df$cell_line)
# luc.df$cell_line<-gsub("TIG111","Human_TIG111",luc.df$cell_line)
# luc.df$cell_line<-gsub("DMR","DMR",luc.df$cell_line)
# luc.df$cell_line<-c(rep("NMR_H17",9),
#                     rep("Ms_6W4",9),
#                     rep("NMR_H17",9),
#                     rep("Ms_6W4",9))
# luc.df$plate<-c(rep("6h",18),
#                 rep("24h",18))
# luc.df$treatment<-rep(c(rep("control",3),
#                     rep("TNF",3),
#                     rep("LPS",3)),4)

# luc.df<-luc.df%>%
#   mutate(Read.2=ifelse(is.na(Read.2),round(Read.1/Ratio),Read.2))
# luc.df.ratio<-luc.df%>%
#   group_by(cell_line,plate,treatment)%>%
#   mutate(avg.ratio=mean(Ratio))
#############


# assay.name<-"20250204_luc1"

###############
# Dual luciferase data ####
assay.name<-"20250214_luc"
barplot.dir<-"./images/luciferase_assay_barplots"
# luc.df<-read.table(file.path("data",paste(assay.name,"luc.tsv",sep="_")))
luc.df<-read.table(file.path("data",paste(assay.name,"tsv",sep=".")),
                   header = T,sep = "\t")
luc.df$Sample<-gsub(" \\/ Repl 1", "", luc.df$Sample)
luc.df$cell_line<-c(rep("Ms_6W3",9),
                    rep("Human_TIG111",9),
                    rep("Ms_6W3",9),
                    rep("Human_TIG111",9))
luc.df$plate<-c(rep("6h",18),
                rep("24h",18))
luc.df$treatment<-rep(c(rep("control",3),
                    rep("TNF",3),
                    rep("LPS",3)),4)
head(luc.df)
# Get values for control wells
luc.df.ctrl.ratio<-luc.df%>%
  filter(treatment=="control")%>%
  group_by(cell_line,plate,treatment)%>%
  mutate(avg.ctrl.ratio=mean(Ratio))%>%
  ungroup%>%
  select(-Sample,-Read.1,-Read.2,-Ratio,-treatment)%>%
  distinct()

#' Calculate the fold change, average fold change, average Read 1, and 
#' average Read 2
luc.df.ratio<-luc.df%>%
  left_join(luc.df.ctrl.ratio)%>%
  group_by(cell_line,plate,treatment)%>%
  mutate(plate=factor(plate,levels=c("6h","24h")),
         cell_line=factor(cell_line,levels=c("Ms_6W3","Human_TIG111")))%>%
  mutate(fold=Ratio/avg.ctrl.ratio,
         avg.fold=mean(fold),
         avg.read1=mean(Read.1),
         avg.read2=mean(Read.2),
         avg.ratio=mean(Ratio))%>%
  ungroup()%>%
  mutate(treatment_numeric=as.numeric(factor(treatment))*2)

head(luc.df.ratio)
luc.assay.plot<-ggplot(luc.df.ratio,aes(x=treatment,
           y=avg.fold,
          fill=cell_line # for some reason, it prevents merging of error bars
))+
  # geom_bar(position="dodge",stat="identity")+
  geom_point(inherit.aes = FALSE, data=luc.df.ratio,
             aes(x=treatment_numeric,
                 y=fold,
                 color=cell_line,
                 shape=cell_line),
             stat="identity",
             size=2,
             stroke=1,
             # shape=21,
             position=position_dodge(width = 0.9))+
  scale_x_continuous(breaks = unique(luc.df.ratio$treatment_numeric),
                     labels = treatment.labels) + # continuous because I want
  # to increase the space between x axis ticks by converting into numeric and 
  # multiplying.
  # scale_x_discrete(label=treatment.labels)+
  scale_y_continuous(breaks = round(seq(0,
                                        max(luc.df.ratio$fold), by = 1),1))+
  facet_grid(~plate,
             labeller = as_labeller(c("6h" = "6 hours",
                                      "24h" = "24 hours")))+
  scale_shape_manual(values = c(16, 0))+
  geom_errorbar(data = luc.df.ratio,
                aes(x = treatment_numeric, 
                    ymin = avg.fold, 
                    ymax = avg.fold),
                width = 0.5, 
                color="black",
                position = position_dodge(width = 0.9),
                linewidth=1)+
  scale_color_manual(values = rep("black",length(unique(luc.df$cell_line))),
                     guide="none")+
  geom_hline(yintercept = 1,
             linetype="dashed",
             color="blue")+
  labs(x="",
       y="Fold change",
       shape="Cell line")+
  coord_cartesian(expand = FALSE,
                  ylim = c(0,max(luc.df.ratio$fold)+0.5))+
  theme_bw()

luc.assay.plot<-luc.assay.plot+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text = element_text(size=15),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        )
ggsave(filename = paste0(barplot.dir,"/",
                        paste(paste(format(Sys.time(),format="%Y%m%d"),
                                    format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "luc_assay",assay.name,sep="-"),".png"),
       plot = luc.assay.plot,
       width = 3000,height = 2000,
       units = "px",dpi=300,device = "png")

## Read.1 plot ####
read1.plot<-ggplot(luc.df.ratio,
       aes(x=treatment,
       y=Read.1,
       fill=cell_line # for some reason, it prevents merging of error bars
))+
  # geom_bar(position="dodge",stat="identity")+
  geom_point(inherit.aes = FALSE, data=luc.df.ratio,
             aes(x=treatment_numeric,
                 y=Read.1,
                 color=cell_line,
                 shape=cell_line),
             stat="identity",
             size=2,
             stroke=1,
             # shape=21,
             position=position_dodge(width = 0.9))+
  scale_x_continuous(breaks = unique(luc.df.ratio$treatment_numeric),
                     labels = treatment.labels) + # continuous because I want
  # to increase the space between x axis ticks by converting into numeric and 
  # multiplying.
  # scale_x_discrete(label=treatment.labels)+
  coord_cartesian(expand = FALSE,
                  ylim = c(0,max(luc.df.ratio$Read.1)+2e5))+
  scale_y_continuous(breaks = round(seq(0,
                                        max(luc.df.ratio$Read.1)+2e5, by = 1e6),1))+
  facet_grid(~plate,
             labeller = as_labeller(c("6h" = "6 hours",
                                      "24h" = "24 hours")))+
  scale_shape_manual(values = c(16, 0))+
  geom_errorbar(data = luc.df.ratio,
                aes(x = treatment_numeric, 
                    ymin = avg.read1, 
                    ymax = avg.read1),
                width = 0.5, 
                color="black",
                position = position_dodge(width = 0.9),
                linewidth=1)+
  scale_color_manual(values = rep("black",length(unique(luc.df$cell_line))),
                     guide="none")+
  geom_hline(yintercept = 1,
             linetype="dashed",
             color="blue")+
  labs(x="",
       y="Read 1",
       shape="Cell line")+
  theme_bw()
read1.plot<-read1.plot+theme(axis.text.x = element_text(size=15),
                 axis.text.y = element_text(size=15),
                 axis.title.y = element_text(size=15),
                 strip.text = element_text(size=15),
                 legend.text = element_text(size=20),
                 legend.title = element_text(size=20),
)


read2.plot<-ggplot(luc.df.ratio,
       aes(x=treatment,
           y=Read.2,
           fill=cell_line # for some reason, it prevents merging of error bars
       ))+
  # geom_bar(position="dodge",stat="identity")+
  geom_point(inherit.aes = FALSE, data=luc.df.ratio,
             aes(x=treatment_numeric,
                 y=Read.2,
                 color=cell_line,
                 shape=cell_line),
             stat="identity",
             size=2,
             stroke=1,
             # shape=21,
             position=position_dodge(width = 0.9))+
  scale_x_continuous(breaks = unique(luc.df.ratio$treatment_numeric),
                     labels = treatment.labels) + # continuous because I want
  # to increase the space between x axis ticks by converting into numeric and 
  # multiplying.
  # scale_x_discrete(label=treatment.labels)+
  coord_cartesian(expand = FALSE,
                  ylim = c(0,max(luc.df.ratio$Read.2)+0.1e8))+
  scale_y_continuous(breaks = round(seq(0,
                                        max(luc.df.ratio$Read.2)+0.1e8, by = 0.2e8),1))+
  facet_grid(~plate,
             labeller = as_labeller(c("6h" = "6 hours",
                                      "24h" = "24 hours")))+
  scale_shape_manual(values = c(16, 0))+
  geom_errorbar(data = luc.df.ratio,
                aes(x = treatment_numeric, 
                    ymin = avg.read2, 
                    ymax = avg.read2),
                width = 0.5, 
                color="black",
                position = position_dodge(width = 0.9),
                linewidth=1)+
  scale_color_manual(values = rep("black",length(unique(luc.df$cell_line))),
                     guide="none")+
  geom_hline(yintercept = 1,
             linetype="dashed",
             color="blue")+
  labs(x="",
       y="Read 2",
       shape="Cell line")+
  theme_bw()
read2.plot<-read2.plot+theme(axis.text.x = element_text(size=15),
                             axis.text.y = element_text(size=15),
                             axis.title.y = element_text(size=15),
                             strip.text = element_text(size=15),
                             legend.text = element_text(size=20),
                             legend.title = element_text(size=20),
)

## Ratio plot ####
ratio.plot<-ggplot(luc.df.ratio,
                   aes(x=treatment,
                       y=Ratio,
                       fill=cell_line # for some reason, it prevents merging of error bars
                   ))+
  # geom_bar(position="dodge",stat="identity")+
  geom_point(inherit.aes = FALSE, data=luc.df.ratio,
             aes(x=treatment_numeric,
                 y=Ratio,
                 color=cell_line,
                 shape=cell_line),
             stat="identity",
             size=2,
             stroke=1,
             # shape=21,
             position=position_dodge(width = 0.9))+
  scale_x_continuous(breaks = unique(luc.df.ratio$treatment_numeric),
                     labels = treatment.labels) + # continuous because I want
  # to increase the space between x axis ticks by converting into numeric and 
  # multiplying.
  # scale_x_discrete(label=treatment.labels)+
  coord_cartesian(expand = FALSE,
                  ylim = c(0,max(luc.df.ratio$Ratio)+0.05))+
  scale_y_continuous(breaks = round(seq(0,
                                        max(luc.df.ratio$Ratio), by = 0.1),1))+
  facet_grid(~plate,
             labeller = as_labeller(c("6h" = "6 hours",
                                      "24h" = "24 hours")))+
  scale_shape_manual(values = c(16, 0))+
  geom_errorbar(data = luc.df.ratio,
                aes(x = treatment_numeric, 
                    ymin = avg.ratio, 
                    ymax = avg.ratio),
                width = 0.5, 
                color="black",
                position = position_dodge(width = 0.9),
                linewidth=1)+
  scale_color_manual(values = rep("black",length(unique(luc.df$cell_line))),
                     guide="none")+
  geom_hline(yintercept = 1,
             linetype="dashed",
             color="blue")+
  labs(x="",
       y="Ratio",
       shape="Cell line")+
  theme_bw()
ratio.plot<-ratio.plot+theme(axis.text.x = element_text(size=15),
                             axis.text.y = element_text(size=15),
                             axis.title.y = element_text(size=15),
                             strip.text = element_text(size=15),
                             legend.text = element_text(size=20),
                             legend.title = element_text(size=20),
)

## Combined plot ####
combined.plot<-ggarrange(read1.plot,read2.plot,
          ratio.plot,luc.assay.plot,
          common.legend = TRUE, legend = "bottom"
)
# ggsave(filename = paste0(barplot.dir,"/",
#                          paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                      format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                "luc_assay-combined",assay.name,sep="-"),".png"),
#        plot = combined.plot,
#        width = 3500,height = 2500,
#        units = "px",dpi=300,device = "png")
combined.plot
# Use this format because we are knitting to R markdown.
# Also, use full path because pandoc doesn't understand the local/global path
# https://forum.posit.co/t/r-markdown-html-document-doesnt-show-image/41629
#' ```{r combined-plot, echo=FALSE, out.width = '100%'}
#' knitr::include_graphics("/home/rakhimov/projects/nfkappab_evo/images/luciferase_assay_barplots/20250304_16_46_49-luc_assay-combined-20250214_luc.png", error = FALSE)
#' ```


# Single luciferase data ####
luminometer.data.dir<-"./data/luminometer_data"
barplot.dir<-"./images/luciferase_assay_barplots"
# assay.date<-"20250423"
assay.date<-"20250507"
assay.date<-"20250612"
assay.name.1<-paste(assay.date,"p1",sep="_")
assay.name.2<-paste(assay.date,"p2",sep="_")
treatment.labels<-c("PBS"="PBS",
                    "murine_TNF"="Murine\nTNF-alpha",
                    "human_TNF"="Human\nTNF-alpha")
# luc.df<-read.table(file.path("data",paste(assay.name,"luc.tsv",sep="_")))
luc.df.1<-read.table(file.path(luminometer.data.dir,paste(assay.name.1,"tsv",sep=".")),
                   header = T,sep = "\t")
luc.df.2<-read.table(file.path(luminometer.data.dir,paste(assay.name.2,"tsv",sep=".")),
                   header = T,sep = "\t")
luc.df.1$Sample<-gsub(" \\/ Repl 1", "", luc.df.1$Sample)
luc.df.2$Sample<-gsub(" \\/ Repl 1", "", luc.df.2$Sample)

xl.file.1<-read_xlsx(path=file.path(luminometer.data.dir,
                                    paste(assay.date,"Amir_1.xlsx",sep="_")))
xl.file.2<-read_xlsx(path=file.path(luminometer.data.dir,
                                    paste(assay.date,"Amir_2.xlsx",sep="_")))
xl.1.filtered<-convert_luc_xl_to_tibble(xl.file.1)
xl.1.filtered$Read<-as.integer(xl.1.filtered$Read)
xl.2.filtered<-convert_luc_xl_to_tibble(xl.file.2)
xl.2.filtered$Read<-as.integer(xl.2.filtered$Read)

identical(xl.1.filtered,as_tibble(luc.df.1))
identical(xl.2.filtered,as_tibble(luc.df.2))

# Exclude samples: 20250507
luc.df.1<-luc.df.1[-c(13,22),]
xl.1.filtered<-xl.1.filtered[-c(13,22),]

# Add sample number
luc.df.1<-luc.df.1%>%
  mutate(sample_num=grep("[0-9]+",Sample))
luc.df.2<-luc.df.2%>%
  mutate(sample_num=grep("[0-9]+",Sample)+nrow(luc.df.1))

luc.df<-rbind(luc.df.1,luc.df.2)%>%
  select(-Sample)%>%
  rename("Sample"=sample_num)%>%
  relocate(Sample)

luc.df.xl<-bind_two_luc_tibbles(xl.1.filtered,xl.2.filtered)
identical(as_tibble(luc.df),luc.df.xl)


luc.df$cell_line<-c(rep("Ms_6W4",9),
                    rep("Human_TIG111",9),
                    rep("NMR",9),
                    rep("DMR",9),
                    rep("Elephant_LACF",9))
# Add treatment labels
# 202050423
# luc.df$treatment<-c("PBS","murine_TNF","human_TNF",
#                     rep("PBS",2),rep("murine_TNF",2),rep("human_TNF",2),
#                     rep(c(rep("PBS",3),rep("murine_TNF",3),rep("human_TNF",3)),4))
luc.df$treatment<-rep(c(rep("PBS",3),rep("murine_TNF",3),rep("human_TNF",3)),5)
head(luc.df)
# Get values for control wells
luc.df.ctrl.ratio<-luc.df%>%
  filter(treatment=="PBS")%>%
  group_by(cell_line,treatment)%>%
  mutate(avg.ctrl.ratio=mean(Read))%>%
  ungroup%>%
  select(-Sample,-Read,-treatment)%>%
  distinct()

#' Calculate the fold change, average fold change, and average Read
luc.df.ratio<-luc.df%>%
  left_join(luc.df.ctrl.ratio)%>%
  group_by(cell_line,treatment)%>%
  mutate(cell_line=factor(cell_line,levels=c("Ms_6W4","Human_TIG111",
                                             "NMR","DMR","Elephant_LACF")))%>%
  mutate(fold=Read/avg.ctrl.ratio,
         avg.fold=mean(fold),
         avg.read=mean(Read))%>%
  ungroup()%>%
  mutate(treatment_numeric=as.numeric(factor(treatment,levels=names(treatment.labels))))

head(luc.df.ratio)
luc.assay.plot<-ggplot(luc.df.ratio,aes(x=treatment,
                                        y=avg.fold,
                                        fill=cell_line # for some reason, it prevents merging of error bars
))+
  # geom_bar(position="dodge",stat="identity")+
  geom_point(inherit.aes = FALSE, data=luc.df.ratio,
             aes(x=treatment_numeric,
                 y=fold,
                 color=cell_line,
                 shape=cell_line),
             stat="identity",
             size=2,
             stroke=1,
             # shape=21,
             position=position_dodge(width = 0.9))+
  scale_x_continuous(breaks = unique(luc.df.ratio$treatment_numeric),
                     labels = factor(treatment.labels,levels =unname(treatment.labels))) + # continuous because I want
  # to increase the space between x axis ticks by converting into numeric and 
  # multiplying.
  # scale_x_discrete(label=treatment.labels)+
  # scale_y_continuous(breaks = round(seq(0,
  #                                       max(luc.df.ratio$fold), by = 10),1))+
  facet_wrap(~cell_line,
             labeller = as_labeller(c("Ms_6W4" = "6W4 (Mouse)",
                                      "Human_TIG111" = "TIG111 (Human)",
                                      "NMR"= "Naked mole-rat",
                                      "DMR" = "Damaraland mole-rat",
                                      "Elephant_LACF" = "LACF (Elephant)")),
             scales = "free",ncol = 2)+
  scale_shape_manual(values = c(0,1,2,5,16))+
  geom_errorbar(data = luc.df.ratio,
                aes(x = treatment_numeric, 
                    ymin = avg.fold, 
                    ymax = avg.fold),
                width = 0.5, 
                color="black",
                position = position_dodge(width = 0.9),
                linewidth=1)+
  scale_color_manual(values = rep("black",length(unique(luc.df$cell_line))),
                     guide="none")+
  geom_hline(yintercept = 1,
             linetype="dashed",
             color="blue")+
  labs(x="",
       y="Fold change",
       shape="Cell line")+
  theme_bw()

luc.assay.plot<-luc.assay.plot+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text = element_text(size=15),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
  )
ggsave(filename = paste0(barplot.dir,"/",
                         paste(paste(format(Sys.time(),format="%Y%m%d"),
                                     format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                               "luc_assay_pHAGE",assay.date,sep="-"),".png"),
       plot = luc.assay.plot,
       width = 3500,height = 3500,
       units = "px",dpi=300,device = "png")


luc.assay.plot.read<-ggplot(luc.df.ratio,aes(x=treatment,
                                        y=Read,
                                        fill=cell_line # for some reason, it prevents merging of error bars
))+
  # geom_bar(position="dodge",stat="identity")+
  geom_point(inherit.aes = FALSE, data=luc.df.ratio,
             aes(x=treatment_numeric,
                 y=Read,
                 color=cell_line,
                 shape=cell_line),
             stat="identity",
             size=2,
             stroke=1,
             # shape=21,
             position=position_dodge(width = 0.9))+
  scale_x_continuous(breaks = unique(luc.df.ratio$treatment_numeric),
                     labels = factor(treatment.labels,levels =unname(treatment.labels))) + # continuous because I want
  # to increase the space between x axis ticks by converting into numeric and 
  # multiplying.
  # scale_x_discrete(label=treatment.labels)+
  # scale_y_continuous(breaks = round(seq(0,
  #                                       max(luc.df.ratio$fold), by = 10),1))+
  facet_wrap(~cell_line,
             labeller = as_labeller(c("Ms_6W4" = "6W4 (Mouse)",
                                      "Human_TIG111" = "TIG111 (Human)",
                                      "NMR"= "Naked mole-rat",
                                      "DMR" = "Damaraland mole-rat",
                                      "Elephant_LACF" = "LACF (Elephant)")),
             scales = "free",ncol = 2)+
  scale_shape_manual(values = c(0,1,2,5,16))+
  geom_errorbar(data = luc.df.ratio,
                aes(x = treatment_numeric, 
                    ymin = avg.read, 
                    ymax = avg.read),
                width = 0.5, 
                color="black",
                position = position_dodge(width = 0.9),
                linewidth=1)+
  scale_color_manual(values = rep("black",length(unique(luc.df$cell_line))),
                     guide="none")+
  geom_hline(yintercept = 1,
             linetype="dashed",
             color="blue")+
  labs(x="",
       y="Read",
       shape="Cell line")+
  theme_bw()

luc.assay.plot.read<-luc.assay.plot.read+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text = element_text(size=15),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
  )
ggsave(filename = paste0(barplot.dir,"/",
                         paste(paste(format(Sys.time(),format="%Y%m%d"),
                                     format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                               "luc_assay_pHAGE_read",assay.date,sep="-"),".png"),
       plot = luc.assay.plot.read,
       width = 3500,height = 3500,
       units = "px",dpi=300,device = "png")

write.table(luc.df.ratio,file = paste0("./output/rtables/",
                                       paste(paste(format(Sys.time(),format="%Y%m%d"),
                                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                       "luc_assay_pHAGE_read",assay.date,sep = "-"),".tsv"),
            row.names = F)
