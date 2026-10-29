
library(tidyverse)
library(phyloseq)
library(vegan)
library(rstatix)
library(forcats)
library(rstatix)
library(ggpubr)

setwd("/mnt/4tb/home/mgimenez/Matias/Cryosmetics")

ps_all <- readRDS("Phyloseq_nometa_EMU_all.RDS")

# Keep patients with both T0 and T1
meta <- as(sample_data(ps_all), "data.frame")

patients_keep <- meta %>%
  group_by(Patient) %>%
  summarise(n_time = n_distinct(Time)) %>%
  filter(n_time == 2) %>%
  pull(Patient)

ps_paired <- prune_samples(meta$Patient %in% patients_keep, ps_all)

ps_all <- ps_paired

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

df_reads <- data.frame(
  Sample = sample_names(ps_all),
  Reads = sample_sums(ps_all)
)

meta <- as(sample_data(ps_all), "data.frame")
meta$Sample <- rownames(meta)

# ==============================
# 1. Shannon diversity
# ==============================

ps_rel <- transform_sample_counts(ps_all, function(x) x / sum(x))

alpha_df <- estimate_richness(ps_rel, measures = "Shannon")

alpha_df$Sample <- rownames(alpha_df)

meta <- as(sample_data(ps_rel), "data.frame")

alpha_df <- merge(alpha_df, meta, by.x = "Sample", by.y = "row.names")

colnames(alpha_df)[4] <- "Time"
colnames(alpha_df)[9] <- "Patient"

alpha_wide <- alpha_df %>%
  select(Patient, Time, Shannon) %>%
  pivot_wider(names_from = Time, values_from = Shannon)

alpha_wide <- alpha_wide %>%
  filter(!is.na(T0) & !is.na(T1))

alpha_wide <- alpha_wide %>%
  mutate(
    delta = T1 - T0
  )

# ==============================
# 2. Shannon change groups (|delta| > 0.2)
# ==============================
threshold <- 0.2

alpha_wide <- alpha_wide %>%
  mutate(
    Change = case_when(
      delta > threshold ~ "Increase",
      delta < -threshold ~ "Decrease",
      TRUE ~ "Stable"
    )
  )

alpha_long <- alpha_wide %>%
  pivot_longer(cols = c(T0, T1), names_to = "Time", values_to = "Shannon")

alpha_df <- alpha_long %>%
  mutate(
    Time = factor(Time, levels = c("T0", "T1"))
  )

# ==============================
# 3. Fig 2A: Shannon T0 vs T1 by group
# ==============================
pvals <- alpha_df %>%
  group_by(Change) %>%
  wilcox_test(Shannon ~ Time, paired = TRUE) %>%
  mutate(
    p_label = paste0("p = ", signif(p, 3))
  )

y_pos <- alpha_df %>%
  group_by(Change) %>%
  summarise(y = max(Shannon, na.rm = TRUE) * 1.08)

pvals <- left_join(pvals, y_pos, by = "Change")

fill_cols <- c(
  "T0" = "#b087c4",
  "T1" = "#87afc4"
)

ggplot(alpha_df, aes(x = Time, y = Shannon, fill = Time)) +
  
  geom_violin(
    trim = FALSE,
    alpha = 0.45,
    color = NA,
    width = 0.9
  ) +
  
  geom_boxplot(
    outlier.shape = NA,
    width = 0.4,
    alpha = 0.75,
    color = "grey30"
  ) +
  
  geom_line(
    aes(group = Patient),
    color = "grey55",
    alpha = 0.5,
    linewidth = 0.6
  ) +
  
  geom_jitter(
    width = 0.08,
    size = 2.8,
    alpha = 0.7,
    color = "grey35"
  ) +
  
  facet_wrap(~ Change, scales = "free_y") +
  
  geom_text(
    data = pvals,
    aes(x = 1.5, y = y, label = p_label),
    inherit.aes = FALSE,
    color="grey25",
    size = 3.85
  ) +
  
  scale_fill_manual(values = fill_cols) +
  
  labs(
    x = NULL,
    y = "Shannon"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "grey25"),
    axis.title = element_text(color = "grey25")
  ) -> Fig2a
  
  Fig2a

# ==============================
# 4. Fig 2B: delta Shannon distribution
# ==============================

delta_df <- alpha_df %>%
  distinct(Patient, Change, delta)

library(ggridges)
library(paletteer)

ggplot(delta_df, aes(x = delta, y = Change, fill = Change)) +
  geom_density_ridges(
    alpha = 0.70,
    color = "grey30",
    linewidth = 0.5,
    bandwidth=0.15
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey35"
  ) +
  labs(
    x = expression(Delta~Shannon~"(T1 - T0)"),
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  scale_fill_paletteer_d("lisa::EdwardHopper") -> Fig2b

Fig2b

# ==============================
# 5. Add change groups to sample metadata
# ==============================

library(dplyr)

change_map <- alpha_wide %>%
  select(Patient, Change)

meta <- as(sample_data(ps_all), "data.frame")

meta$Sample <- rownames(meta)

meta2 <- meta %>%
  left_join(change_map, by = "Patient")

rownames(meta2) <- meta2$Sample

meta2$Sample <- NULL

meta2 <- meta2[sample_names(ps_all), , drop = FALSE]

all(rownames(meta2) == sample_names(ps_all))

sample_data(ps_all) <- sample_data(meta2)

ps_all <- prune_samples(!is.na(sample_data(ps_all)$Change), ps_all)

ps <- prune_samples(!is.na(sample_data(ps_all)$Change), ps_all)

# ==============================
# 6. Prevalence filter (>= 10%)
# ==============================
prev_threshold <- 0.10

prev_taxa <- apply(otu_table(ps), 1, function(x) sum(x > 0) / length(x))

ps <- prune_taxa(prev_taxa >= prev_threshold, ps)

# ==============================
# 7. CLR + Aitchison distance
# ==============================
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

ps_clr <- transform_sample_counts(ps_rel, function(x) {
  log(x + 1) - mean(log(x + 1))
})

aitchison_dist <- distance(ps_clr, method = "euclidean")

# ==============================
# 8. NMDS
# ==============================
set.seed(123)
ordu <- ordinate(ps_clr, method = "NMDS", distance = aitchison_dist)
ordu$stress

ord_df <- as.data.frame(ordu$points)
ord_df$Sample <- rownames(ord_df)

meta <- as(sample_data(ps_clr), "data.frame")
meta$Sample <- rownames(meta)

ord_df <- merge(ord_df, meta, by = "Sample")

# ==============================
# 9. Fig 2C: NMDS plot
# ==============================
ggplot(ord_df, aes(x = MDS1, y = MDS2)) +
  
  geom_path(
    aes(group = Patient),
    color = "gray55",
    alpha = 0.45,
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(fill = Change, shape = Time),
    size = 3.5,
    color = "gray25",
    stroke = 0.6,
    alpha = 0.7
  ) +
  
  stat_ellipse(
    aes(group = Time, linetype = Time),
    level = 0.85,
    linewidth = 0.8,
    color = "gray50"
  ) +
  
  scale_shape_manual(
    values = c(
      "T0" = 21,
      "T1" = 24
    )
  ) +
  
  scale_fill_paletteer_d("lisa::EdwardHopper")+
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        color = "gray25"
      )
    )
  ) +
  
  theme_minimal(base_size = 13) +
  
  labs(
    subtitle = paste("Stress =", round(ordu$stress, 3)),
    x = "NMDS1",
    y = "NMDS2"
  ) -> Fig2c

Fig2c

# ==============================
# 10. Pairwise distances
# ==============================

aitch_mat <- as.matrix(aitchison_dist)

meta <- as(sample_data(ps_clr), "data.frame")
meta$Sample <- rownames(meta)

dist_df <- as.data.frame(as.table(aitch_mat)) %>%
  rename(Sample1 = Var1,
         Sample2 = Var2,
         Distance = Freq) %>%
  filter(Sample1 != Sample2) %>%
  left_join(meta %>% select(Sample, Patient, Time, Change),
            by = c("Sample1" = "Sample")) %>%
  rename(Patient1 = Patient,
         Time1 = Time,
         Change1 = Change) %>%
  left_join(meta %>% select(Sample, Patient, Time, Change),
            by = c("Sample2" = "Sample")) %>%
  rename(Patient2 = Patient,
         Time2 = Time,
         Change2 = Change)

dist_df <- dist_df %>%
  mutate(pair_id = paste(
    pmin(as.character(Sample1), as.character(Sample2)),
    pmax(as.character(Sample1), as.character(Sample2)),
    sep = "_"
  )) %>%
  distinct(pair_id, .keep_all = TRUE)

dist_3groups <- dist_df %>%
  mutate(
    comparison = case_when(
      Patient1 != Patient2 & Time1 == "T0" & Time2 == "T0" ~ "Between individuals T0",
      Patient1 != Patient2 & Time1 == "T1" & Time2 == "T1" ~ "Between individuals T1",
      Patient1 == Patient2 & Time1 != Time2 ~ "Within individual T0-T1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(comparison)) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c(
        "Between individuals T0",
        "Within individual T0-T1",
        "Between individuals T1"
      )
    )
  )

dist_3groups2 <- dist_3groups %>%
  mutate(
    comparison = case_when(
      Time1 == "T0" & Time2 == "T0" ~ "Between individuals T0",
      Time1 == "T1" & Time2 == "T1" ~ "Between individuals T1",
      Time1 != Time2 & Patient1 == Patient2 ~ "Within individuals T0-T1",
      TRUE ~ NA_character_
    ),
    comparison = factor(
      comparison,
      levels = c(
        "Between individuals T0",
        "Between individuals T1",
        "Within individuals T0-T1"
      )
    )
  ) %>%
  filter(!is.na(comparison), !is.na(Distance))

# ==============================
# 11. Fig 2D: mean between-individual vs within-individual distance
# ==============================

mean_between_ind <- dist_df %>%
  filter(Patient1 != Patient2) %>%
  filter(
    (Time1 == "T0" & Time2 == "T0") |
      (Time1 == "T1" & Time2 == "T1")
  ) %>%
  bind_rows(
    dist_df %>%
      filter(Patient1 != Patient2) %>%
      filter(
        (Time1 == "T0" & Time2 == "T0") |
          (Time1 == "T1" & Time2 == "T1")
      ) %>%
      transmute(
        Sample1 = Sample2,
        Patient1 = Patient2,
        Time1 = Time2,
        Change1 = Change2,
        Sample2 = Sample1,
        Patient2 = Patient1,
        Time2 = Time1,
        Change2 = Change1,
        Distance = Distance
      )
  ) %>%
  group_by(Sample1, Patient1, Time1, Change1) %>%
  summarise(
    mean_distance = mean(Distance, na.rm = TRUE),
    n_comparisons = n(),
    .groups = "drop"
  ) %>%
  transmute(
    Sample = Sample1,
    Patient = Patient1,
    Change = Change1,
    comparison = paste0("Mean between individuals ", Time1),
    Distance = mean_distance
  )

within_patient <- dist_3groups2 %>%
  filter(comparison == "Within individuals T0-T1") %>%
  transmute(
    Sample = Sample1,
    Patient = Patient1,
    Change = Change1,
    comparison = "Within individuals T0-T1",
    Distance = Distance
  )

mean_plus_within <- bind_rows(
  mean_between_ind,
  within_patient
) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c(
        "Mean between individuals T0",
        "Mean between individuals T1",
        "Within individuals T0-T1"
      )
    )
  )

# Paired Wilcoxon by patient
mean_plus_within <- mean_plus_within %>%
  arrange(comparison, Patient)

pairwise_res <- mean_plus_within %>%
  pairwise_wilcox_test(
    Distance ~ comparison,
    paired = TRUE,
    p.adjust.method = "BH"
  ) %>%
  slice(c(1, 3, 2)) %>%  # longest bracket on top
  add_xy_position(x = "comparison", step.increase = 0.4)

ggplot(mean_plus_within, aes(x = comparison, y = Distance, fill = comparison)) +
  
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    color = "gray25",
    linewidth = 0.5,
    alpha = 0.85
  ) +
  
  ggbeeswarm::geom_quasirandom(
    width = 0.15,
    size = 2.2,
    alpha = 0.7,
    color = "gray30"
  ) +
  
  stat_pvalue_manual(
    pairwise_res,
    label = "p.adj",
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 3.35
  ) +
  labs(y = "Community divergence")+
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"))+
  scale_fill_manual(values=c("#b087c4", "#87afc4","#A4804CFF")) -> Fig2d
 
 Fig2d

 
 
# ==============================
# 12. Final layout
# ==============================
 library(patchwork)
 
 top_row <- Fig2a + Fig2b +
   plot_layout(widths = c(2, 1))
 
 bottom_row <- Fig2c + Fig2d +
   plot_layout(widths = c(1, 1))
 
 Fig2_final <- top_row / bottom_row +
   plot_annotation(tag_levels = "A")

 Fig2_final

 ggsave(
   file.path(
     "/mnt/4tb/home/mgimenez/Matias/Cryosmetics/Cryosmetics",
      "Figure_2.png"
   ),
   Fig2_final,
   width = 15,
   height = 10,
   dpi = 300,
   bg = "white",
   device = ragg::agg_png
 )

