rm(list = ls())

data <- c(G620_1 = 2, G583_3 = 277,  G583_2 = 185, G583_1 = 71, G549_1 = 25, G523_1 = 325)

# Create a color palette
library(RColorBrewer)
colors <- brewer.pal(length(data), "Blues")
#colors <- brewer.pal(length(data), "Reds")

# Enhanced Horizontal Barplot
bp <- barplot(data,
              horiz = TRUE,          # Horizontal bars
              col = "#6BAED6",          # Gradient color palette
              border = "darkblue",   # Dark blue border for bars
              xlab = "Count",        # X-axis label
              main = NULL, # Title
              las = 1,               # Horizontal axis labels
              xlim = c(0, max(data) + 50)) # Extend axis range for labels

# Add values on the bars
text(x = data,
     y = bp,
     labels = data,
     pos = 4,                # Position to the right of the bars
     cex = 0.8,              # Text size
     col = "black")          # Text color

# Add gridlines
abline(v = seq(0, max(data), by = 50), col = "lightgray", lty = "dotted")
