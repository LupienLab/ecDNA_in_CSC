rm(list = ls())

library(dplyr)
library(stringr)
library(data.table)
library(JunJunZai)

id = "G523_L-1"
motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alter_motifs <- toupper(c("", "", "SOX3,Sox11,Sox4,RUNX1,Sox2,Sox6", "Smad3,Smad4", "",                                      #5
                         "Twist2", "", "Oct6,POU2F2,POU2F1", "Klf12", "",                                                    #10
                         "SOX10,Sox21,Sox14,Sox12", "", "POU6F2", "POU3F2,Klf7", "Tcfap2e,SOX3,SOX4,SOX11"                   #15
                        ))

motif_df <- motif_table %>% mutate(motif = alter_motifs)
motif_df <- motif_df %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])

save(motif_df, file = paste0("./outputs/", id, ".motif.RData"))

#~~~
id = "G549_L-1"
motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alter_motifs <- toupper(c("", "SOX10", "SOX10", "TFAP2C,TFAP2A,Smad3", "",                                        #5
                         "FOSL2,FOSL1", "KLF10,Smad3,Smad4", "", "Pou4f3", "",                                #10
                         "", "", "Neurog1,Twist2", "Klf7,RUNX2,Smad3,Klf4,RUNX1", "",                             #15
                         "TFAP2A,Tcfap2b,TFAP2C", "Sox5", "TFAP2B,Tcfap2e,TFAP2A,TFAP2C,Tcfap2b", "", "",     #20
                         "", "", "", "", "NeuroD1",                                                               #25
                         "RUNX2", "", "Sox10,Sox14,Sox6,Sox3,SOX9", "", "",                                       #30
                         "POU6F2,Pou4f3,Pou6f1", "", "", "KLF16,KLF14,KLF9", "",                                  #35
                         "", "", "", "Tbx20", "",                                                                      #40
                         "NeuroD1", "Tcf3", "SOX4,Sox11"))

motif_df <- motif_table %>% mutate(motif = alter_motifs)
motif_df <- motif_df %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])
sort(unique(motif_df$motif))

save(motif_df, file = paste0("./outputs/", id, ".motif.RData"))

#~~~
id = "G583_L-1"
motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alter_motifs <- toupper(c("", "RUNX1,RUNX2", "", "AP-1,Fosl2", "SOX3,Sox6,Sox2,Sox12,Sox4,Sox10",                   #5
                         "", "RUNX1,RUNX2,Sox17", "Tbox:Smad", "Sox17", "",                                         #10
                         "Sox8,SOX10,Pou5f1::Sox2,Sox6"))

motif_df <- motif_table %>% mutate(motif = alter_motifs)
motif_df <- motif_df %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])
sort(unique(motif_df$motif))

save(motif_df, file = paste0("./outputs/", id, ".motif.RData"))

#~~~
id = "G583_L-2"
motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alter_motifs <- toupper(c("TCF3,TCF4", "TBX15", "", "", "",                             #5
                         "", "Tcf3", "SMAD3", "Sox4,SOX10", "",                             #10
                         "", "Sox17,Sox1", "Tcfap2e,Sox4,Sox11", "", "Tcfap2e,Tcf3",        #15
                         "", "POU6F2", "", "TCF4,TCF3,TBX1,TBX5,TBX15"                                 #19
                         ))

motif_df <- motif_table %>% mutate(motif = alter_motifs)
motif_df <- motif_df %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])
sort(unique(motif_df$motif))

save(motif_df, file = paste0("./outputs/", id, ".motif.RData"))

#~~~
id = "G583_L-3"
motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alter_motifs <- toupper(c("", "TCF4,TCF3", "", "", "",                                                                             #5
                         "Klf7,Klf12", "Tcf4,Tcf3", "Sox6,Sox11,Sox4,Sox3,Sox10,Sox15", "", "Tcfap2e,SOX11,Sox4",           #10
                         "Sox1,Sox17", "Sox5", "", "", "Smad2,Smad4,Smad3",                                                             #15
                         "Tcf21,Tcf12", "Klf4,KLF5,KLF1", "", "Tcfap2e,Tcf3"))         

motif_df <- motif_table %>% mutate(motif = alter_motifs)
motif_df <- motif_df %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])
sort(unique(motif_df$motif))

save(motif_df, file = paste0("./outputs/", id, ".motif.RData"))

#~~~ 
id = "G620_L-1"
motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/ATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alter_motifs <- toupper(c("", "", "", "POU6F2,Pou6f1", "KLF10,KLF9",                                                                               #5
                         "", "NeuroD1,Ascl2,Ascl1,Tcf3", "", "Sox6,Sox15,Sox14,Sox2,Sox17,Sox3,Sox10,Sox9", "",                                    #10
                         "Brn1", "", "", "", "",                                                                                                   #15
                         "Eomes,TCF3,TCF4", "Tbx5,TBX4", "", "", "",                                                                                        #20
                         "RUNX2,RUNX1,Klf4", "RUNX2,RUNX1", "POU5F1,Oct4,Sox17,POU2F2,Oct6", "Sox7,Sox15,Sox21,Sox8,Sox30,Sox9,Sox2", "POU6F1",    #25
                         "RUNX2,RUNX1", "", "Sox14", "Sox3,Sox12,Sox6,Tbox:Smad,Sox14,Sox4,SOX10"                                                            #30
                         ))

motif_df <- motif_table %>% mutate(motif = alter_motifs)
motif_df <- motif_df %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])
sort(unique(motif_df$motif))

save(motif_df, file = paste0("./outputs/", id, ".motif.RData"))


