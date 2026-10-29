


library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(readxl)
library(paletteer)
library(rstatix)
library(ggpubr)


setwd("/mnt/4tb/home/mgimenez/Matias/Cryosmetics")

ps_all <- readRDS("Phyloseq_nometa_EMU_all.RDS")

#ps_all <- readRDS("Phyloseq_nometa_EMU_all_March26.RDS")

# Keep patients with both T0 and T1
meta <- as(sample_data(ps_all), "data.frame")

patients_keep <- meta %>%
  group_by(Patient) %>%
  summarise(n_time = n_distinct(Time)) %>%
  filter(n_time == 2) %>%
  pull(Patient)


ps_paired <- prune_samples(meta$Patient %in% patients_keep, ps_all)

ps_all <- ps_paired



# Add metadata

dataset <- read_excel("Tabla_metadatos_Cryos.xlsx")

mtd <- dataset[2:29,c(1:3)]
colnames(mtd) <- c("Patients", "Sex", "Age")

mtd$Patients <- gsub("[.]", "_", mtd$Patients)


sd <- data.frame(sample_data(ps_all))

mtd2 <- mtd %>%
  rename(Patient = Patients)

sd_new <- sd %>%
  left_join(mtd2, by = "Patient")

rownames(sd_new) <- rownames(sd)

sd_new$Age[which(sd_new$Patient=="CRYO_046")] <- 21
sd_new$Sex[which(sd_new$Patient=="CRYO_046")] <- "F"

sd_new$Age[which(sd_new$Patient=="CRYO_027")] <- 19
sd_new$Age[which(sd_new$Patient=="CRYO_030")] <- 67
sd_new$Age[which(sd_new$Patient=="CRYO_031")] <- 30
sd_new$Age[which(sd_new$Patient=="CRYO_032")] <- 35
sd_new$Age[which(sd_new$Patient=="CRYO_034")] <- 32


sample_data(ps_all) <- sample_data(sd_new)


# ==============================
# 0. Batch correction (MMUPHin, Time as covariate)
# ==============================

library(MMUPHin)

correct_batch <- TRUE

if (correct_batch) {

  meta_batch <- data.frame(sample_data(ps_all))
  meta_batch$Batch <- factor(meta_batch$Batch)
  meta_batch$Time <- factor(meta_batch$Time)

  counts <- as(otu_table(ps_all), "matrix")
  if (!taxa_are_rows(ps_all)) counts <- t(counts)

  fit_batch <- adjust_batch(
    feature_abd = counts,
    batch = "Batch",
    covariates = "Time",
    data = meta_batch,
    control = list(verbose = FALSE)
  )

  otu_table(ps_all) <- otu_table(
    round(fit_batch$feature_abd_adj),
    taxa_are_rows = TRUE
  )
  ps_all <- prune_taxa(taxa_sums(ps_all) > 0, ps_all)
}


library(phyloseq)
library(tidyverse)
library(paletteer)

# ==============================
# 1. Age groups
# ==============================

meta <- as(sample_data(ps_all), "data.frame")

meta <- meta %>%
  mutate(
    Age = as.numeric(Age),
    Age_group = case_when(
      Age <= 30 ~ "<= 30 years",
      Age >= 31 & Age <= 40 ~ "31-40 years",
      Age >= 41 ~ ">= 41 years",
      TRUE ~ NA_character_
    ),
    Age_group = factor(Age_group, levels = c("<= 30 years", "31-40 years", ">= 41 years"))
  )

rownames(meta) <- sample_names(ps_all)
sample_data(ps_all) <- sample_data(meta)

# ==============================
# 2. T0 samples, relative abundance
# ==============================

ps_t0 <- prune_samples(
  sample_data(ps_all)$Time == "T0" &
    !is.na(sample_data(ps_all)$Age),
  ps_all
)

ps_rel <- transform_sample_counts(ps_t0, function(x) x / sum(x))

# ==============================
# 3. Long format
# ==============================

df_species <- psmelt(ps_rel) %>%
  mutate(
    Species_plot_name = as.character(OTU)
  )

# ==============================
# 4. Top species by mean abundance
# ==============================

top15_species <- df_species %>%
  group_by(Species_plot_name) %>%
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 12) %>%
  pull(Species_plot_name)

# ==============================
# 5. Remaining taxa as Others
# ==============================

df_barplot <- df_species %>%
  mutate(
    Species_plot = ifelse(
      Species_plot_name %in% top15_species,
      Species_plot_name,
      "Others"
    ),
    Species_plot = factor(
      Species_plot,
      levels = c(top15_species, "Others")
    )
  ) %>%
  group_by(Sample, Age_group, Species_plot) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

# ==============================
# 6. Core microbiome (T0)
# ==============================

ps_rel_t0 <- transform_sample_counts(
  ps_t0,
  function(x) x / sum(x)
)

core_ps_t0 <- core(
  ps_rel_t0,
  detection = 0,
  prevalence = 0.99
)

core_taxa <- taxa_names(core_ps_t0)

core_abund_df <- psmelt(core_ps_t0) %>%
  mutate(
    Species = as.character(OTU)
  )

taxa_order <- core_abund_df %>%
  group_by(Species) %>%
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) %>%
  pull(Species)

core_abund_df <- core_abund_df %>%
  mutate(
    Species = factor(Species, levels = rev(taxa_order))
  )



# ==============================
# 7. Species colors
# ==============================

all_species <- unique(c(top15_species, taxa_order))

species_cols <- setNames(
  c(as.character(paletteer::paletteer_d(
    "rcartocolor::Antique",
    n = 12
  )),"#bf9fad"),
  all_species
)

species_cols_barplot <- c(
  species_cols,
  "Others" = "grey75"
)


# ==============================
# 8. Fig 1A: barplot
# ==============================

Fig_1a <- ggplot(
  df_barplot,
  aes(x = Sample, y = Abundance, fill = Species_plot)
) +
  geom_col(width = 0.85) +
  facet_grid(
    ~ Age_group,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(values = species_cols_barplot) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

Fig_1a


# ==============================
# 9. Fig 1B: Shannon by age group
# ==============================

alpha_shannon <- estimate_richness(
  ps_t0,
  measures = "Shannon"
)

alpha_shannon$Sample <- rownames(alpha_shannon)

meta_t0 <- as(sample_data(ps_t0), "data.frame")
meta_t0$Sample <- rownames(meta_t0)

alpha_shannon <- alpha_shannon %>%
  left_join(
    meta_t0 %>% select(Sample, Age_group),
    by = "Sample"
  ) %>%
  mutate(
    Age_group = factor(
      Age_group,
      levels = c("<= 30 years", "31-40 years", ">= 41 years")
    )
  )

p_shannon <- alpha_shannon %>%
  pairwise_wilcox_test(
    Shannon ~ Age_group,
    p.adjust.method = "BH"
  ) %>%
  slice(c(1, 3, 2)) %>%  # longest bracket on top
  add_xy_position(x = "Age_group", step.increase = 0.25)

Fig_1B <- ggplot(
  alpha_shannon,
  aes(x = Age_group, y = Shannon, fill = Age_group)
) +
  
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    color = "gray25",
    alpha = 0.65
  ) +
  
  ggbeeswarm::geom_quasirandom(
    width = 0.15,
    size = 2.2,
    alpha = 0.7,
    color = "gray30"
  ) +
  
  stat_pvalue_manual(
    p_shannon,
    label = "p.adj",
    hide.ns = FALSE,
    tip.length = 0.01,
    size = 3
  ) +
  
  scale_fill_manual(values = c(
    "<= 30 years" = "#FED789FF",
    "31-40 years" = "#023743FF",
    ">= 41 years" = "#453947FF"
  )) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black")
  ) +
  
  labs(
    y = "Shannon diversity"
  )

Fig_1B


library(phyloseq)
library(vegan)
library(tidyverse)
library(paletteer)

# ==============================
# 10. CLR + Aitchison distance
# ==============================

prev_t0 <- apply(otu_table(ps_t0), 1, function(x) mean(x > 0))
ps_t0_prev <- prune_taxa(prev_t0 >= 0.10, ps_t0)

ps_clr <- transform_sample_counts(
  ps_t0_prev,
  function(x) log(x + 1) - mean(log(x + 1))
)

aitchison_dist <- phyloseq::distance(
  ps_clr,
  method = "euclidean"
)

# ==============================
# 11. NMDS
# ==============================

set.seed(123)
ordu <- ordinate(
  ps_clr,
  method = "NMDS",
  distance = aitchison_dist
)

ord_df <- as.data.frame(ordu$points)
ord_df$Sample <- rownames(ord_df)

meta_t0 <- as(sample_data(ps_clr), "data.frame")
meta_t0$Sample <- rownames(meta_t0)

ord_df <- left_join(
  ord_df,
  meta_t0,
  by = "Sample"
)

# ==============================
# 12. PERMANOVA
# ==============================

set.seed(123)
perm_age <- adonis2(
  aitchison_dist ~ Age_group,
  data = meta_t0,
  permutations = 999
)

perm_age

p_perm <- signif(perm_age$`Pr(>F)`[1], 3)
r2_perm <- signif(perm_age$R2[1], 3)

# ==============================
# 13. Fig 1C: NMDS plot
# ==============================

Fig_1c <- ggplot(
  ord_df,
  aes(x = MDS1, y = MDS2)
) +
  
  stat_ellipse(
    aes(color = Age_group),
    linewidth = 1,
    level = 0.80,
    alpha = 0.8
  ) +
  
  geom_point(
    aes(fill = Age_group),
    shape = 21,
    size = 4,
    color = "gray20",
    alpha = 0.65,
    stroke = 0.6
  ) +
  
  scale_fill_manual(values = c(
    "<= 30 years" = "#FED789FF",
    "31-40 years" = "#023743FF",
    ">= 41 years" = "#453947FF"
  )) +
  
  scale_color_manual(values = c(
    "<= 30 years" = "#FED789FF",
    "31-40 years" = "#023743FF",
    ">= 41 years" = "#453947FF"
  )) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.5,
    label = paste0(
      "PERMANOVA\nR² = ",
      r2_perm,
      "\np = ",
      p_perm
    ),
    size = 4
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  
  labs(
    subtitle = paste("NMDS stress =", round(ordu$stress, 3)),
    x = "NMDS1",
    y = "NMDS2"
  )

Fig_1c


library(microbiome)
library(phyloseq)
library(tidyverse)
library(ggbeeswarm)

# ==============================
# 14. Fig 1D: core taxa (prevalence 99%)
# ==============================

ps_rel_t0 <- transform_sample_counts(
  ps_t0,
  function(x) x / sum(x)
)

core_ps_t0 <- core(
  ps_rel_t0,
  detection = 0,
  prevalence = 0.99
)

core_taxa <- taxa_names(core_ps_t0)

core_taxa
length(core_taxa)

core_abund_df <- psmelt(core_ps_t0) %>%
  mutate(
    Species = as.character(OTU),
    Species = factor(Species, levels = core_taxa)
  )

taxa_order <- core_abund_df %>%
  group_by(Species) %>%
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) %>%
  pull(Species)

core_abund_df <- core_abund_df %>%
  mutate(
    Species = factor(
      Species,
      levels = rev(taxa_order)
    )
  )

Fig_1d<- ggplot(
  core_abund_df,
  aes(x = Species, y = Abundance, fill = Species)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    color = "gray25",
    alpha = 0.85
  ) +
  
  ggbeeswarm::geom_quasirandom(
    width = 0.18,
    size = 1.8,
    alpha = 0.65,
    color = "gray30"
  ) +
  
  scale_fill_manual(values = species_cols) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1)
  ) +
  
  coord_flip() +
  
  theme_minimal(base_size = 12) +
  
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black")
  ) +
  
  labs(
    x = NULL,
    y = "Relative abundance"
  )

Fig_1d

library(patchwork)

# ==============================
# 15. Final layout
# ==============================


top_row <- Fig_1a +
  Fig_1B +
  plot_layout(widths = c(2, 1))

bottom_row <- Fig_1c +
  Fig_1d +
  plot_layout(widths = c(1, 1.3))

Final_Figure <- top_row /
  bottom_row +
  plot_annotation(tag_levels = "A")

Final_Figure


ggsave(
  file.path(
    "/mnt/4tb/home/mgimenez/Matias/Cryosmetics/Cryosmetics",
    "Figure_1.png"
  ),
  Final_Figure,
  width = 15,
  height = 10,
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)
