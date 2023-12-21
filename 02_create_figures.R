# Recreation of Figure for Annelisa
library(tidyverse)

project_dir = "/hpc/pmc_vanheesch/projects/Jip/counts"

# Load quantification data
plotdata_wide <- read.delim(
  file = paste(project_dir,
               "mr1_hla_plotdata.txt", sep = "/"),
  sep = "\t",
  header = T
)

# Colour schemes
embryo_col <- c("#A5CAEB", "#F8B0AF")

# Make the various plots for the figures
p3 <-
  ggplot(data = plotdata_wide, aes(
    x = MR1,
    y = HLA_A,
    col = V2,
    shape = embryonal
  )) +
  geom_point() +
  scale_x_continuous(
    trans = 'log2',
    breaks = scales::trans_breaks("log2", function(x)
      2 ^ x),
    labels = scales::trans_format("log2",
                                  scales::math_format(2 ^ .x))
  ) +
  scale_y_continuous(
    trans = 'log2',
    breaks = scales::trans_breaks("log2", function(x)
      2 ^ x),
    labels = scales::trans_format("log2",
                                  scales::math_format(2 ^ .x))
  ) +
  labs(y = "HLA-A expression",
       x = "MR1 expression",
       col = "Tumor type",
       shape = "Tumor class") +
  scale_shape_manual(values = c(3, 18)) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 1)
  )

p4 <-
  ggplot(data = plotdata_wide, aes(
    x = MR1,
    y = HLA_B,
    col = V2,
    shape = embryonal
  )) +
  geom_point() +
  scale_x_continuous(
    trans = 'log2',
    breaks = scales::trans_breaks("log2", function(x)
      2 ^ x),
    labels = scales::trans_format("log2", scales::math_format(2 ^ .x))
  ) +
  scale_y_continuous(
    trans = 'log2',
    breaks = scales::trans_breaks("log2", function(x)
      2 ^ x),
    labels = scales::trans_format("log2", scales::math_format(2 ^ .x))
  ) +
  labs(y = "HLA-B expression",
       x = "MR1 expression",
       col = "Tumor type",
       shape = "Tumor class") +
  scale_shape_manual(values = c(3, 18)) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 1)
  )

p5 <-
  ggplot(data = plotdata_wide, aes(
    x = MR1,
    y = HLA_C,
    col = V2,
    shape = embryonal
  )) +
  geom_point() +
  scale_x_continuous(
    trans = 'log2',
    breaks = scales::trans_breaks("log2", function(x)
      2 ^ x),
    labels = scales::trans_format("log2", scales::math_format(2 ^ .x))
  ) +
  scale_y_continuous(
    trans = 'log2',
    breaks = scales::trans_breaks("log2", function(x)
      2 ^ x),
    labels = scales::trans_format("log2", scales::math_format(2 ^ .x))
  ) +
  labs(y = "HLA-C expression",
       x = "MR1 expression",
       col = "Tumor type",
       shape = "Tumor class") +
  scale_shape_manual(values = c(3, 18)) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 1)
  )

p1_violin <-
  ggplot(data = plotdata, aes(
    x = reorder(V2, CPM),
    y = CPM,
    fill = V2
  )) +
  geom_hline(yintercept = 20) +
  geom_violin(scale = "width") +
  ggbeeswarm::geom_quasirandom(
    fill = "white",
    col = "black",
    shape = 21,
    size = 1.2
  ) +
  ggpubr::stat_compare_means() +
  labs(y = "MR1 expression (counts per million)") +
  
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.background = element_rect(colour = "black", size = 1)
  )

p2 <-
  ggplot(data = plotdata, aes(x = embryonal, y = CPM, fill = embryonal)) +
  geom_hline(yintercept = 20) +
  geom_violin(scale = "width") +
  ggbeeswarm::geom_quasirandom(
    fill = "white",
    col = "black",
    shape = 21,
    size = 0.9
  ) +
  labs(y = "MR1 expression (counts per million)") +
  scale_fill_manual(values = embryo_col) +
  ggpubr::stat_compare_means() +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.background = element_rect(colour = "black", size = 1)
  )

# Correlations presented on the figures
cor.test(log2(plotdata_wide$MR1), log2(plotdata_wide$HLA_A), method = "pearson")
cor.test(log2(plotdata_wide$MR1), log2(plotdata_wide$HLA_B), method = "pearson")
cor.test(log2(plotdata_wide$MR1), log2(plotdata_wide$HLA_C), method = "pearson")

# Save plots
ggsave(
  filename = "MR1_CPM_per_tumortype_withpoints.pdf",
  device = "pdf",
  plot = p1_violin,
  height = 6,
  width = 8,
  path = project_dir
)

ggsave(
  filename = "MR1_CPM_embryonal_withpoints.pdf",
  device = "pdf",
  plot = p2,
  height = 3,
  width = 4,
  path = project_dir
)

ggsave(
  filename = "HLA-A_CPM_embryonal.pdf",
  device = "pdf",
  plot = p3,
  height = 6,
  width = 8,
  path = project_dir
)

ggsave(
  filename = "HLA-B_CPM_embryonal.pdf",
  device = "pdf",
  plot = p4,
  height = 6,
  width = 8,
  path = project_dir
)
ggsave(
  filename = "HLA-C_CPM_embryonal.pdf",
  device = "pdf",
  plot = p5,
  height = 6,
  width = 8,
  path = project_dir
)

ggsave(
  filename = "HLA_CPM_per_tumortype.pdf",
  device = "pdf",
  plot = ggpubr::ggarrange(
    p3,
    p4,
    p5,
    ncol = 3,
    nrow = 1,
    common.legend = T,
    legend = "right"
  ),
  height = 6,
  width = 18,
  path = project_dir
)
