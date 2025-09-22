rm(list = ls())

library(dplyr)
library(ggplot2)
library(data.table)

id = "A7TK_T"
nth = 2

#~~~ input cell type information (ecDNA-positive vs. ecDNA-negative)
ecDNA_positive_cell <- read.table(paste0(id, ".outputs/", id, "-", nth, "_ecDNA_positive_cell.txt"), header = F, stringsAsFactors = F) %>% {colnames(.) = "barcode";.} %>% pull("barcode")
ecDNA_negative_cell <- read.table(paste0(id, ".outputs/", id, "-", nth, "_ecDNA_negative_cell.txt"), header = F, stringsAsFactors = F) %>% {colnames(.) = "barcode";.} %>% pull("barcode")

data <- data.frame(
  sample = rep(paste0(id, "-", nth), each = 1),
  category = rep(c("1", "2"), times = 1),
  number = c(length(ecDNA_positive_cell), length(ecDNA_negative_cell))
)

data_percentage <- data %>% group_by(sample) %>% mutate(percentage = number / sum(number) * 100)
data_percentage <- data_percentage %>% filter(sample == paste0(id, "-", nth))

if(id == "A7TK_T") colorPalette = c("1" = "salmon", "2" = "#FFA07A")
if(id == "G933_T") colorPalette = c("1" = "orange", "2" = "#FFDAB9")
if(id == "G946_T") colorPalette = c("1" = "orchid", "2" = "#F2D4E6")
if(id == "G958_T") colorPalette = c("1" = "purple", "2" = "#E6E6FA")
celltype = c("ecDNA(+)", "ecDNA(-)", "normal-like") 

pdf(file = paste0(id, ".outputs/", id, "-", nth, "_ecDNA_percentage.pdf"), width = 5, height = 1)

p <- ggplot(data_percentage, aes(x = sample, y = percentage, fill = category)) + theme_bw() + 
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = colorPalette, labels = celltype) +
  geom_hline(yintercept = data_percentage$percentage[2], linetype = "dashed", color = "black", linewidth = 0.25) + 
  geom_text(aes(label = paste0(round(percentage, 2), "%")), position = position_stack(vjust = 0.5), size = 3) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(title = NULL, x = NULL, y = "Percentage (%)", fill = "cell type") 

print(p)

dev.off()
