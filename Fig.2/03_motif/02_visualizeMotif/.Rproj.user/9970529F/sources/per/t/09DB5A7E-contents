rm(list = ls())

library(dplyr)
library(stringr)
library(data.table)
library(JunJunZai)

id = "A7TK_T-1"

motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/homerResults/"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'
alterMotifs <- toupper(c("", "", "Tcf3,Tcf4", "Klf4", "",                                    #5
                         "", "Sox8", "", "Sox4,Sox11", "",                                   #10
                         "", "", "Tcf4", "", "FOSL1",                                        #15
                         "Eomes,Klf7", "Smad2,Smad4", "", "NeuroD1,Tcf12", "",               #20
                         "", "Sox5", "Pou6f1", "", "Klf7",                                   #25
                         "Smad3", "Sox13,Klf4,Eomes", "RUNX", "Smad4,Smad2", "",             #30
                         "", "Klf7,Klf4", "", "", "Smad2",                                   #35
                         "Sox3,Sox17,Sox14,Sox9,Sox2"))

motif_table <- motif_table %>% mutate(motif = alterMotifs) 
motif_table <- motif_table %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])
sort(unique(motif_table$motif))

save(motif_table, file = paste0("./outputs/", id, ".motif.RData"))

#~~~
id = "A7TK_T-2"

motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alterMotifs <- toupper(c("", "", "", "", "Smad3",                                                            #5
                         "Klf4", "", "Smad3", "", "RUNX2",                                                        #10
                         "Tcfap2a", "", "", "Smad2,Smad4,Smad3", "Tcfap2e,Tcfap2c,Tcfap2a,TFAP2A",           #15
                         "", "Sox11,Sox13,Sox4,Sox5,Sox8", "", "SOX10,Sox3", "AP-1,Tcf4",
                         "Sox12,Sox4,Sox11,Sox5"))               

motif_table <- motif_table %>% mutate(motif = alterMotifs)
motif_table <- motif_table %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])

save(motif_table, file = paste0("./outputs/", id, ".motif.RData"))

#~~~
id = "G933_T-1"

motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alterMotifs <- toupper(c("", "Tcfap2a,Tcfap2c,Tcfap2b,TFAP2A", "TFAP2A,Ascl2,Klf4", "Sox6,Sox2", "Klf7",    #5
                         "Klf7,Ascl2,Klf4", "Sox6,Sox12,Sox11,Sox7,Sox18"))

motif_table <- motif_table %>% mutate(motif = alterMotifs)
motif_table <- motif_table %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])

save(motif_table, file = paste0("./outputs/", id, ".motif.RData"))

#~~~
id = "G933_T-2"

motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alterMotifs <- toupper(c("Sox8", "", "TFAP2A", "Sox14,Sox9", "Oct2,Oct4,Pou3f3,Pou5f1,Pou2f3,Pou2f2",             #5
                         "", "", "Smad4", "Sox17,Sox10,Sox6,AP-1,Sox3,Sox30", "RUNX2,RUNX1",
                         "Smad3", "", "Eomes,RUNX2"))

motif_table <- motif_table %>% mutate(motif = alterMotifs)
motif_table <- motif_table %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])

save(motif_table, file = paste0("./outputs/", id, ".motif.RData"))

#~~~
id = "G933_T-3"

motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alterMotifs <- toupper(c("", "", "Smad3,TFAP2A,Smad2", "", "Tcfap2b,Tcfap2e",                                             #5
                         "", "AP-1", "", "RUNX2,Smad4,Smad2", "TFAP2A,Smad3",     #10
                         "", "", "Tcfap2e", "", "",                                     #15
                         "", "Smad2,Smad4", "", "Smad4,Smad2,Klf7,Smad3", "",                                                                                  #20
                         "", "", "", "", "Smad3",                                                                          #25
                         "", "", "Pou2f3,Oct4,Oct2", "", "RUNX1,RUNX2",                      #30
                         "", "Tcf3,Tcf4", "Sox14", "", ""))

motif_table <- motif_table %>% mutate(motif = alterMotifs)
motif_table <- motif_table %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])

save(motif_table, file = paste0("./outputs/", id, ".motif.RData"))


#~~~
id = "G946_T-1"

motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alterMotifs <- toupper(c("Pou5f1,Pou3f3,Oct4", "RUNX2,RUNX1", "", "Eomes", "NeuroD1",      #5
                         "Klf4,Eomes,Klf7", "", "Eomes", "", "Sox14,Sox7",                                              #10 
                         "", "", "Smad3,Ascl2", "", "Tcfap2e,Sox11,Sox4,Tcf3",                     #15
                         "RUNX2,RUNX1", "AP-1"))

motif_table <- motif_table %>% mutate(motif = alterMotifs)
motif_table <- motif_table %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])

save(motif_table, file = paste0("./outputs/", id, ".motif.RData"))


#~~~
id = "G958_T-1"

motif_files <- list.files(path = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/homerResults"), pattern = "^motif[0-9]+\\.motif$", full.names = TRUE)
motif_table <- loadHomerRes(homerDir = paste0("/Users/chupan/Documents/data_su2c/scATAC/motif/", id, "/results/"), novo = T, known = F, motifIndex = 1:length(motif_files))@'novo res'

alterMotifs <- toupper(c("", "Sox30", "", "AP-1", "",                                   #5
                         "Tcf3,Tcf4", "", "", "", "",
                         "Smad3", "RUNX1", "Tcf3,Tcf4", "", "",                                    #15
                         "Sox15", "Klf7", "", "Sox5,Sox8,Sox7,Sox18,Sox12,,Sox17", "Smad4,Smad2,Smad3,Tcf3",                               #20
                         "", "Tcfap2e,Tcf3", "Sox17,Klf4", "Smad4,Ascl2,TFAP2A", "Pou2f2",
                         "Smad3,NeuroD1,Eomes", "Smad3,SOX10", "Tcfap2e,Sox4,Sox11,Sox12,Sox10"))

motif_table <- motif_table %>% mutate(motif = alterMotifs)
motif_table <- motif_table %>% arrange(-neglogPvalue) %>% tidyr::separate_rows(motif, sep = ",") %>% mutate(rank = 1:dim(.)[1])

save(motif_table, file = paste0("./outputs/", id, ".motif.RData"))

