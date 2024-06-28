# load packages
library(ggplot2)

# load csv file
args <- commandArgs(trailingOnly = TRUE)
seq_path = args[1]
seq <- read.csv(seq_path)

filename <- args[2]

# make plot
P = ggplot(seq, aes(x = aa_number, y = hydropathy)) +
  ggtitle(args[2]) +
  geom_line(aes(group = 1), size = 0.2) +
  geom_point(aes(color = group), size = 0.5) +
  xlim(0, 600) +  
  ylim(-4, 3) +
  scale_color_manual(values = c("signal_peptide" = "#AA4499", "cysteine" = "#964B00", "KEX2" = "#0072B2", "CPP" = "#50C878", "DREK" = "#FF0000", "DREK+CPP" = "#FDDA0D", "Black"))

outfile <- paste(filename, "png", sep=".")
ggsave(outfile, P, dpi=500, height = 7, width = 12)