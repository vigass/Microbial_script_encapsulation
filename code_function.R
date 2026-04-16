rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggsci)
# 绝对丰度转相对丰度
otu2oturel <- function(data_input) {
  
  if (any(colSums(data_input) == 0)) {
    stop("Some samples have total abundance = 0.")
  }
  
  sweep(data_input, 2, colSums(data_input), "/") %>%
    as.data.frame()
}

# 获取top taxa
get_top_taxa <- function(input, n = 20) {
  if (n > nrow(input)) {
    warning("n is larger than number of taxa; returning all taxa.")
    n <- nrow(input)
  }
  
  input %>%
    rowMeans() %>%
    sort(decreasing = TRUE) %>%
    head(n) %>%
    names()
}

# 画图准备
make_plot_df <- function(input, sample_info, top_taxa) {
  
  input %>%
    rownames_to_column("Taxon") %>%
    pivot_longer(
      cols = -Taxon,
      names_to = "Sample_ID",
      values_to = "Relative_Abundance"
    ) %>%
    left_join(sample_info, by = "Sample_ID") %>%
    mutate(
      Taxon = if_else(Taxon %in% top_taxa, Taxon, "Others")
    ) %>%
    group_by(Sample_ID, Groups, Taxon) %>%
    summarise(
      Relative_Abundance = sum(Relative_Abundance),
      .groups = "drop"
    )
  
}
make_group_df <- function(plot_df_sample) {
  plot_df_sample %>%
    group_by(Groups, Taxon) %>%
    summarise(
      Relative_Abundance = mean(Relative_Abundance),
      .groups = "drop"
    )
}

# 画图
plot_taxa_sample_bar <- function(plot_df_sample, fill_colors) {
  my_theme <- theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 12, face = "italic"),
      strip.background = element_rect(fill = "gray90"), 
      text = element_text(family = "Arial")
    )
  p <- ggplot(
    plot_df_sample,
    aes(
      x = Sample_ID,
      y = Relative_Abundance,
      fill = Taxon
    )
  ) +
    geom_bar(stat = "identity", width = 0.6) +
    facet_wrap(~Groups, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = fill_colors) + 
    my_theme +
    theme(
      axis.text.x = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 1))
  
  return(p)
}

plot_taxa_group_bar <- function(plot_df_group, fill_colors) {
  
  my_theme <- theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 12, face = "italic"),
      strip.background = element_rect(fill = "gray90")
    )
  
  p <- ggplot(
    plot_df_group,
    aes(
      x = Groups,
      y = Relative_Abundance,
      fill = Taxon
    )
  ) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = fill_colors) +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(fill = guide_legend(ncol = 1))
  
  return(p)
}


# 数据导入
sample_info <- readxl::read_xlsx("~/Downloads/sample_map.xlsx") %>%
  setNames(c("Sample_ID", "Groups"))

otu <- readxl::read_xlsx("~/Downloads/Species_count_format.xlsx")

otu_matrix <- otu %>%
  column_to_rownames("Species")

otu_rel <- otu2oturel(otu_matrix)
top_taxa <- get_top_taxa(otu_rel, n = 20)
plot_df_sample <- make_plot_df(otu_rel, sample_info, top_taxa)
plot_df_group <- make_group_df(plot_df_sample)
# plot_df_group <- plot_df_sample %>%
#   group_by(Groups, Taxon) %>%
#   summarise(
#     Relative_Abundance = mean(Relative_Abundance),
#     .groups = "drop"
#   )

final_cols <- c(
  setNames(
    pal_d3("category20")(length(top_taxa)),
    top_taxa
  ),
  Others = "grey60"
)

group_levels <- c(
  "Int8_25d", "Int8_35d",
  "Int9_25d", "Int9_35d",
  "NC_25d",   "NC_35d",
  "PC_25d",   "PC_35d"
)

plot_df_sample <- plot_df_sample %>% 
  mutate(
    Taxon = factor(Taxon, levels = rev(c(top_taxa, "Others"))),
    Groups = factor(Groups, levels = group_levels)
  )

p1 <- plot_taxa_sample_bar(plot_df_sample, final_cols)
plot(p1)


plot_df_group <- plot_df_group %>% 
  mutate(
    Taxon = factor(Taxon, levels = rev(c(top_taxa, "Others"))),
    Groups = factor(Groups, levels = group_levels)
  )
p2 <- plot_taxa_group_bar(plot_df_group, final_cols)
plot(p2)
