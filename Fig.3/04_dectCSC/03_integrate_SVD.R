#~~~~ step3: perform SVD on the integrated matrix ~~~#

rm(list = ls())

library(irlba)
library(Matrix)
library(dplyr)
library(stringr)
library(ggplot2)
set.seed(1)

peak_by_cell_GBM <- readRDS("./outputs/peak_by_cell_GBM.rds")
peak_by_cell_GSC <- readRDS("./outputs/peak_by_cell_GSC.rds")

peak_by_cell_matrix <- cbind(peak_by_cell_GBM, peak_by_cell_GSC)
peak_by_cell_matrix@x[peak_by_cell_matrix@x > 0] <- 1

#~~~ transpose matrix to have cells as rows
mat_t <- t(peak_by_cell_matrix)

#~~~ perform SVD using irlba (computing the first 5 components)
res <- irlba(mat_t, nv = 5)

scores <- res$u %*% diag(res$d)
colnames(scores) <- paste0("SV", 1:ncol(scores))

#~~~ compute variance explained
var_explained <- res$d^2 / sum(res$d^2)
ve_df <- data.frame(
  component = paste0("SV", seq_along(var_explained)),
  variance  = var_explained
)

res_df_23 <- scores %>% as.data.frame() %>% dplyr::select(SV2, SV3) %>% 
  mutate(barcode = colnames(peak_by_cell_matrix)) %>% mutate(type = str_extract(barcode, "(?<=_)[A-Z]"), sample = str_extract(barcode, "^[^#]+"))

#~~~ visualization
dir.create("outputs.SVD")
pdf(paste0("./outputs.SVD/SV23_integrate.pdf"), width = 3.5, height = 3.5)

ggplot(res_df_23, aes(x = SV2, y = SV3, colour = type)) +
  geom_point(size = 0.25) +
  theme_bw() +
  scale_color_manual(values = c("T" = "red", "L" = "purple"), name = "cell type", labels = c("GBM", "CSC")) +
  labs(title = "",
       x = sprintf("SV2: %.1f%% variance", 100 * var_explained[2]),
       y = sprintf("SV3: %.1f%% variance", 100 * var_explained[3])) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

saveRDS(res_df_23, file = paste0("outputs/integrate_SVD.rds"))

# Because you binarized the matrix (@x > 0 → 1), 
# most cells likely share many accessible peaks, so there’s a dominant global pattern. 
# SV1 is capturing that overwhelming signal (e.g., tumor vs GSC difference, library size bias, or global accessibility trends). 
# The other components capture residual heterogeneity.


res_df_12 <- scores %>% as.data.frame() %>% dplyr::select(SV1, SV2) %>% 
  mutate(barcode = colnames(peak_by_cell_matrix)) %>% mutate(type = str_extract(barcode, "(?<=_)[A-Z]"), sample = str_extract(barcode, "^[^#]+"))

pdf(paste0("./outputs.SVD/SV12_integrate.pdf"), width = 3.5, height = 3.5)

ggplot(res_df_12, aes(x = SV1, y = SV2, colour = type)) +
  geom_point(size = 0.25) +
  theme_bw() +
  scale_color_manual(values = c("T" = "red", "L" = "purple"), name = "cell type", labels = c("GBM", "CSC")) +
  labs(title = "",
       x = sprintf("SV1: %.1f%% variance", 100 * var_explained[1]),
       y = sprintf("SV2: %.1f%% variance", 100 * var_explained[2])) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"))

dev.off()

res_df_13 <- scores %>% as.data.frame() %>% dplyr::select(SV1, SV3) %>% 
  mutate(barcode = colnames(peak_by_cell_matrix)) %>% mutate(type = str_extract(barcode, "(?<=_)[A-Z]"), sample = str_extract(barcode, "^[^#]+"))

pdf(paste0("./outputs.SVD/SV13_integrate.pdf"), width = 3.5, height = 3.5)

ggplot(res_df_13, aes(x = SV1, y = SV3, colour = type)) +
  geom_point(size = 0.25) +
  theme_bw() +
  scale_color_manual(values = c("T" = "red", "L" = "purple"), name = "cell type", labels = c("GBM", "CSC")) +
  labs(title = "",
       x = sprintf("SV1: %.1f%% variance", 100 * var_explained[1]),
       y = sprintf("SV3: %.1f%% variance", 100 * var_explained[3])) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"))

dev.off()