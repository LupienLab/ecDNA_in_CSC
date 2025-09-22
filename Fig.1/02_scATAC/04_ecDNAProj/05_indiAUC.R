#******** step 3: calculate the AUC of ROC and PRC curves for identifying normal cells as ecDNA-negative ********#

rm(list = ls())
    
library(pROC)
library(PRROC)
library(dplyr)
library(data.table)
    
id = "G958_T"
n = 1
    
content_score <- fread(paste0("/Users/chupan/Documents/data_su2c/scATAC/celltype/", id, "-", n, ".ecDNA_content.txt"), header = F, data.table = F, stringsAsFactors = F) %>% 
  {colnames(.) = c("barcode", "uscore");.} %>% mutate(uscore = 1 - uscore)
    
filter_cell <- fread(paste0("/Users/chupan/Documents/data_su2c/scATAC/celltype/", id, ".filtercell.txt"), header = F, data.table = F, stringsAsFactors = F) %>% {colnames(.) = "barcode";.} 
malign_cell <- fread(paste0("/Users/chupan/Documents/data_su2c/scATAC/celltype/", id, ".maligncell.txt"), header = F, data.table = F, stringsAsFactors = F) %>% {colnames(.) = "barcode";.} 
normal_cell <- anti_join(filter_cell, malign_cell, by = "barcode") %>% mutate(type = 1)
    
result <- content_score %>% left_join(normal_cell, by = "barcode") %>% mutate(type = if_else(!is.na(type), 1, 0)) 

roc_curve <- roc(result$type, result$uscore)
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))
    
#~~~ plot the (Receiver Operating Characteristic) ROC curve
plot(roc_curve, col = "blue", lwd = 2, cex.max = 1, cex.lab = 1, cex.axis = 1, main = paste0(id, "-", n, ", AUC =", format(auc_value, nsmall = 2)))
    
#~~~ calculate PRC (Precision-Recall Curve)
true_labels <- result$type
predicted_scores <- result$uscore
pr <- pr.curve(scores.class0 = predicted_scores[true_labels == 1], scores.class1 = predicted_scores[true_labels == 0], curve = TRUE)
    
#~~~ plot the PRC curve
plot(pr, col = "blue", lwd = 2, auc.main = TRUE)
    
