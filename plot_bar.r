rm(list = ls())

library(tidyverse)
library(data.table)
library(readxl)
library(ggsci)

#========================
# 1. Read Data
#========================
sample_info <- read_xlsx("~/Downloads/4.11sample_map.xlsx") %>%
  setNames(c("Sample_ID", "Groups"))
sample_info$Sample_ID <- sample_info$Groups
sample_info <- sample_info %>% 
  distinct()

sample_info <- sample_info %>% 
  mutate(Groups = case_when(
    str_detect(Groups, "^NC") ~ "NC", 
    str_detect(Groups, "^PC") ~ "PC", 
    str_detect(Groups, "^Int") ~ "INT",
    TRUE ~ Groups
  ))
sample_info$Groups <- factor(sample_info$Groups, 
                             levels = c("NC", "PC", "INT"))

otu <- read_xlsx("~/Downloads/4.11全taxonomy_Species_abund_Group.xlsx")

otu_matrix <- otu %>%
  column_to_rownames("Species")

#========================
# 2. Relative Abundance
#========================
otu_rel <- sweep(otu_matrix, 2, colSums(otu_matrix), "/") %>%
  as.data.frame()

#========================
# 3. Top20 Taxa
#========================
top_num <- 20
top_taxa <- otu_rel %>%
  rowMeans() %>%
  sort(decreasing = TRUE) %>%
  head(top_num) %>%
  names()

#========================
# 4. Long Format
#========================
plot_df <- otu_rel %>%
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

#========================
# 5. Factor Levels
#========================
sample_levels <- c(
  "NC_25d",   "NC_35d",
  "PC_25d",   "PC_35d", 
  "Int8_25d", "Int8_35d",
  "Int9_25d", "Int9_35d"
)

group_levels <- c("NC", "PC", "INT")

plot_df <- plot_df %>%
  mutate(
    Taxon = factor(Taxon, levels = rev(c(top_taxa, "Others"))),
    Groups = factor(Groups, levels = group_levels), 
    Sample_ID = factor(Sample_ID, levels = sample_levels) 
  )

#========================
# 6. Colors
#========================
final_cols <- c(
  setNames(
    pal_d3("category20")(length(top_taxa)),
    top_taxa
  ),
  Others = "grey60"
)

#========================
# 7. Sample-level Barplot
#========================
p_sample <- ggplot(
  plot_df,
  aes(
    x = Sample_ID,
    y = Relative_Abundance,
    fill = Taxon
  )
) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = final_cols, name = "Species") +
  # facet_wrap(~Groups, scales = "free_x", nrow = 1) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    legend.text = element_text(face = "italic"), 
    text = element_text(family = "Arial")
  ) + 
  guides(fill = guide_legend(ncol = 1))

plot(p_sample)

# ggsave(
#   "~/Desktop/Part_time_job/2026_04_08/plot_bar_sample.tiff",
#   p_sample,
#   dpi = 600,
#   width = 10,
#   height = 7
# )
# otu_rel1 <- otu_rel %>% 
#   rownames_to_column("Species") %>% 
#   select("Species", everything())
# write_csv(otu_rel1,file = "~/Desktop/Part_time_job/2026_04_08/otu_rel.csv")

# #========================
# # 8. Group Mean Barplot
#========================
group_df <- plot_df %>%
  group_by(Groups, Taxon) %>%
  summarise(
    Relative_Abundance = mean(Relative_Abundance),
    .groups = "drop"
  )

p_group <- ggplot(
  group_df,
  aes(
    x = Groups,
    y = Relative_Abundance,
    fill = Taxon
  )
) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = final_cols, name = "Species") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.text = element_text(face = "italic"),
    text = element_text(family = "Arial")
  ) +
  guides(fill = guide_legend(ncol = 1))

plot(p_group)

# ggsave(
#   "~/Desktop/Part_time_job/2026_04_08/plot_bar_group.tiff",
#   p_group,
#   dpi = 600,
#   width = 10,
#   height = 7
# )
# 
# #========================
# # 9. Save Output
# #========================
save(
  otu,
  otu_rel,
  sample_info,
  file = "~/Desktop/Part_time_job/2026_04_08/output.rdata"
)
