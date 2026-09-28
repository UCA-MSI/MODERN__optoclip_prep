#!/usr/bin/env Rscript
# =====================================================================
# peak_length_hist.R — P6 peak-length histogram (per biological sample)
# ---------------------------------------------------------------------
# Reproduces the draft ggplot block: peak length = V3 - V2 of the .pk
# file, histogram, saved as PDF. Driven by the Snakemake `script:`
# directive (uses the snakemake@ S4 object; no manual arg parsing).
# =====================================================================
suppressMessages(library(ggplot2))

pk_file  <- snakemake@input[["pk"]]
pdf_file <- snakemake@output[["pdf"]]
p        <- snakemake@params[["plot"]]

peaks <- read.delim(pk_file, header = FALSE)
peaks$V10 <- peaks$V3 - peaks$V2
peaks <- as.data.frame(peaks)

g <- ggplot(peaks, aes(V10)) +
  geom_histogram(binwidth = p$binwidth, fill = I("blue"), col = I("black")) +
  coord_cartesian(xlim = c(p$xlim_min, p$xlim_max)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

ggsave(pdf_file, plot = g,
       width = p$width_in, height = p$height_in, units = "in", dpi = p$dpi)
