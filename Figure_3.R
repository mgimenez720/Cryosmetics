

library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(readxl)
library(paletteer)
library(dplyr)
library(Hmisc)
library(igraph)
library(tidygraph)
library(ggraph)
library(patchwork)
library(rstatix)
library(ggpubr)


setwd("/mnt/4tb/home/mgimenez/Matias/Cryosmetics")



ps_all <- readRDS("EMU_meta_physeq.RDS")


# ==============================
# 2. Filtrar T0 y transformar a relativa
# ==============================

ps_t0 <- prune_samples(
  sample_data(ps_all)$Time == "T0" &
    !is.na(sample_data(ps_all)$Age),
  ps_all
)


# Abundancia relativa
ps_t0_rel <- transform_sample_counts(ps_t0, function(x) x / sum(x))

# Matriz taxa x samples
otu <- as(otu_table(ps_t0_rel), "matrix")

if (!taxa_are_rows(ps_t0_rel)) {
  otu <- t(otu)
}

# occupancy y abundancia media
taxa_ecology_t0 <- data.frame(
  taxon = rownames(otu),
  occupancy = rowSums(otu > 0) / ncol(otu),
  mean_abundance = rowMeans(otu),
  median_abundance = apply(otu, 1, median),
  max_abundance = apply(otu, 1, max)
) %>%
  mutate(
    category = case_when(
      occupancy >= 0.80 ~ "Core",
      occupancy >= 0.30 & occupancy < 0.80 ~ "Intermediate",
      occupancy < 0.30 ~ "Accessory"
    )
  )


taxa_for_network_t0 <- taxa_ecology_t0 %>%
  filter(
    occupancy >= 0.10,
    mean_abundance >= 0.0001
  )


# 1. Extraer nombres de taxa que pasan el filtro
taxa_keep_t0 <- taxa_for_network_t0$taxon

# 2. Subsetear phyloseq T0 con esos taxa
ps_net_t0 <- prune_taxa(taxa_keep_t0, ps_t0)

# 3. Remover taxa y muestras vacías
ps_net_t0 <- prune_taxa(taxa_sums(ps_net_t0) > 0, ps_net_t0)
ps_net_t0 <- prune_samples(sample_sums(ps_net_t0) > 0, ps_net_t0)

# 4. Extraer matriz taxa x samples
otu_t0 <- as(otu_table(ps_net_t0), "matrix")

if (!taxa_are_rows(ps_net_t0)) {
  otu_t0 <- t(otu_t0)
}

# 5. CLR con pseudocount
otu_t0_clr <- t(apply(otu_t0 + 1, 2, function(x) {
  log(x / exp(mean(log(x))))
}))

# ahora queda samples x taxa

# 6. Correlaciones Spearman
cor_t0 <- Hmisc::rcorr(
  as.matrix(otu_t0_clr),
  type = "spearman"
)

rho_t0 <- cor_t0$r
p_t0   <- cor_t0$P

# 7. Corrección BH de p-values
p_adj_t0 <- matrix(
  p.adjust(p_t0, method = "BH"),
  nrow = nrow(p_t0),
  ncol = ncol(p_t0),
  dimnames = dimnames(p_t0)
)

# 8. Pasar a tabla de edges
edge_table_t0 <- as.data.frame(as.table(rho_t0)) %>%
  rename(
    edge_from = Var1,
    edge_to = Var2,
    rho = Freq
  ) %>%
  mutate(
    pvalue = as.vector(p_t0),
    padj = as.vector(p_adj_t0)
  ) %>%
  filter(
    edge_from != edge_to,
    rho > 0.6,
    padj < 0.05
  )

# 9. Evitar duplicados A-B / B-A
edge_table_t0 <- edge_table_t0 %>%
  rowwise() %>%
  mutate(
    pair_id = paste(sort(c(edge_from, edge_to)), collapse = "___")
  ) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(-pair_id)

# 10. Crear red
g_t0_pos <- graph_from_data_frame(
  d = edge_table_t0,
  directed = FALSE
)

# 11. Agregar atributos de nodos
node_attr_t0 <- taxa_for_network_t0 %>%
  select(taxon, occupancy, mean_abundance, median_abundance, max_abundance, category)

V(g_t0_pos)$occupancy <- node_attr_t0$occupancy[
  match(V(g_t0_pos)$name, node_attr_t0$taxon)
]

V(g_t0_pos)$mean_abundance <- node_attr_t0$mean_abundance[
  match(V(g_t0_pos)$name, node_attr_t0$taxon)
]

V(g_t0_pos)$category <- node_attr_t0$category[
  match(V(g_t0_pos)$name, node_attr_t0$taxon)
]

# 12. Agregar métricas de centralidad
V(g_t0_pos)$degree <- degree(g_t0_pos)

V(g_t0_pos)$strength <- strength(
  g_t0_pos,
  weights = E(g_t0_pos)$rho
)

V(g_t0_pos)$betweenness <- betweenness(
  g_t0_pos,
  directed = FALSE,
  weights = 1 / E(g_t0_pos)$rho
)

V(g_t0_pos)$eigenvector <- eigen_centrality(
  g_t0_pos,
  directed = FALSE,
  weights = E(g_t0_pos)$rho
)$vector


node_table_t0 <- data.frame(
  taxon = V(g_t0_pos)$name,
  category = V(g_t0_pos)$category,
  occupancy = V(g_t0_pos)$occupancy,
  mean_abundance = V(g_t0_pos)$mean_abundance,
  degree = V(g_t0_pos)$degree,
  strength = V(g_t0_pos)$strength,
  betweenness = V(g_t0_pos)$betweenness,
  eigenvector = V(g_t0_pos)$eigenvector
)


# 1. Subset T1
ps_t1 <- subset_samples(ps_all, Time == "T1")
ps_t1 <- prune_taxa(taxa_sums(ps_t1) > 0, ps_t1)
ps_t1 <- prune_samples(sample_sums(ps_t1) > 0, ps_t1)

# 2. Abundancia relativa para calcular occupancy y mean abundance
ps_t1_rel <- transform_sample_counts(ps_t1, function(x) x / sum(x))

otu_t1_rel <- as(otu_table(ps_t1_rel), "matrix")

if (!taxa_are_rows(ps_t1_rel)) {
  otu_t1_rel <- t(otu_t1_rel)
}

# 3. Tabla ecológica T1
taxa_ecology_t1 <- data.frame(
  taxon = rownames(otu_t1_rel),
  occupancy = rowSums(otu_t1_rel > 0) / ncol(otu_t1_rel),
  mean_abundance = rowMeans(otu_t1_rel),
  median_abundance = apply(otu_t1_rel, 1, median),
  max_abundance = apply(otu_t1_rel, 1, max)
) %>%
  mutate(
    category_t1 = case_when(
      occupancy >= 0.80 ~ "Core",
      occupancy >= 0.30 ~ "Intermediate",
      TRUE ~ "Accessory"
    )
  )

taxa_for_network_t1 <- taxa_ecology_t1 %>%
  filter(
    occupancy >= 0.10,
    mean_abundance >= 0.0001
  )


# 4. Extraer taxa que pasan filtro
taxa_keep_t1 <- taxa_for_network_t1$taxon

# 5. Subsetear phyloseq T1
ps_net_t1 <- prune_taxa(taxa_keep_t1, ps_t1)

# 6. Remover taxa y muestras vacías
ps_net_t1 <- prune_taxa(taxa_sums(ps_net_t1) > 0, ps_net_t1)
ps_net_t1 <- prune_samples(sample_sums(ps_net_t1) > 0, ps_net_t1)

# 7. Extraer matriz taxa x samples
otu_t1 <- as(otu_table(ps_net_t1), "matrix")

if (!taxa_are_rows(ps_net_t1)) {
  otu_t1 <- t(otu_t1)
}

# 8. CLR con pseudocount
otu_t1_clr <- t(apply(otu_t1 + 1, 2, function(x) {
  log(x / exp(mean(log(x))))
}))

# ahora queda samples x taxa

# 9. Correlaciones Spearman
cor_t1 <- Hmisc::rcorr(
  as.matrix(otu_t1_clr),
  type = "spearman"
)

rho_t1 <- cor_t1$r
p_t1   <- cor_t1$P

# 10. Corrección BH de p-values
p_adj_t1 <- matrix(
  p.adjust(p_t1, method = "BH"),
  nrow = nrow(p_t1),
  ncol = ncol(p_t1),
  dimnames = dimnames(p_t1)
)

# 11. Pasar a tabla de edges
edge_table_t1 <- as.data.frame(as.table(rho_t1)) %>%
  rename(
    edge_from = Var1,
    edge_to = Var2,
    rho = Freq
  ) %>%
  mutate(
    pvalue = as.vector(p_t1),
    padj = as.vector(p_adj_t1)
  ) %>%
  filter(
    edge_from != edge_to,
    rho > 0.6,
    padj < 0.05
  )

# 12. Evitar duplicados A-B / B-A
edge_table_t1 <- edge_table_t1 %>%
  rowwise() %>%
  mutate(
    pair_id = paste(sort(c(edge_from, edge_to)), collapse = "___")
  ) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(-pair_id)

# 13. Crear red
g_t1_pos <- graph_from_data_frame(
  d = edge_table_t1,
  directed = FALSE
)

# 14. Agregar atributos de nodos
node_attr_t1 <- taxa_for_network_t1 %>%
  select(taxon, occupancy, mean_abundance, median_abundance, max_abundance, category_t1)

V(g_t1_pos)$occupancy_t1 <- node_attr_t1$occupancy[
  match(V(g_t1_pos)$name, node_attr_t1$taxon)
]

V(g_t1_pos)$mean_abundance_t1 <- node_attr_t1$mean_abundance[
  match(V(g_t1_pos)$name, node_attr_t1$taxon)
]

V(g_t1_pos)$category_t1 <- node_attr_t1$category_t1[
  match(V(g_t1_pos)$name, node_attr_t1$taxon)
]

# 15. Centralidades
V(g_t1_pos)$degree_t1 <- degree(g_t1_pos)

V(g_t1_pos)$strength_t1 <- strength(
  g_t1_pos,
  weights = E(g_t1_pos)$rho
)

V(g_t1_pos)$betweenness_t1 <- betweenness(
  g_t1_pos,
  directed = FALSE,
  weights = 1 / E(g_t1_pos)$rho
)

V(g_t1_pos)$eigenvector_t1 <- eigen_centrality(
  g_t1_pos,
  directed = FALSE,
  weights = E(g_t1_pos)$rho
)$vector


node_table_t1 <- data.frame(
  taxon = V(g_t1_pos)$name,
  category_t1 = V(g_t1_pos)$category_t1,
  occupancy_t1 = V(g_t1_pos)$occupancy_t1,
  mean_abundance_t1 = V(g_t1_pos)$mean_abundance_t1,
  degree_t1 = V(g_t1_pos)$degree_t1,
  strength_t1 = V(g_t1_pos)$strength_t1,
  betweenness_t1 = V(g_t1_pos)$betweenness_t1,
  eigenvector_t1 = V(g_t1_pos)$eigenvector_t1
)


#Crear categorías basadas en la ocupancia en T0

taxa_category_baseline <- taxa_ecology_t0 %>%
  select(taxon, category) %>%
  rename(category_baseline = category)


V(g_t0_pos)$category_baseline <- taxa_category_baseline$category_baseline[
  match(V(g_t0_pos)$name, taxa_category_baseline$taxon)
]

V(g_t1_pos)$category_baseline <- taxa_category_baseline$category_baseline[
  match(V(g_t1_pos)$name, taxa_category_baseline$taxon)
]

V(g_t0_pos)$category_baseline[is.na(V(g_t0_pos)$category_baseline)] <- "Not detected in T0"
V(g_t1_pos)$category_baseline[is.na(V(g_t1_pos)$category_baseline)] <- "Not detected in T0"





node_table_t0 <- data.frame(
  taxon = V(g_t0_pos)$name,
  degree = degree(g_t0_pos),
  strength = strength(g_t0_pos, weights = E(g_t0_pos)$rho),
  betweenness = betweenness(g_t0_pos, directed = FALSE, weights = 1 / E(g_t0_pos)$rho),
  time = "T0"
)

node_table_t1 <- data.frame(
  taxon = V(g_t1_pos)$name,
  degree = degree(g_t1_pos),
  strength = strength(g_t1_pos, weights = E(g_t1_pos)$rho),
  betweenness = betweenness(g_t1_pos, directed = FALSE, weights = 1 / E(g_t1_pos)$rho),
  time = "T1"
)


taxa_category_baseline <- taxa_ecology_t0 %>%
  select(taxon, category) %>%
  rename(category_baseline = category)

node_centrality <- bind_rows(node_table_t0, node_table_t1) %>%
  left_join(taxa_category_baseline, by = "taxon") %>%
  mutate(
    category_baseline = ifelse(
      is.na(category_baseline),
      "Not detected in T0",
      category_baseline
    ),
    category_baseline = factor(
      category_baseline,
      levels = c("Core", "Intermediate", "Accessory", "Not detected in T0")
    ),
    time = factor(time, levels = c("T0", "T1"))
  )


wilcox_res_degree <- node_centrality %>%
  group_by(category_baseline) %>%
  wilcox_test(
    degree ~ time
  ) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

wilcox_res_degree <- wilcox_res_degree %>%
  add_xy_position(
    x = "category_baseline"
  )

wilcox_res_degree




#### 
#Negative correlation graphs
####
 
 
 edge_table_t0_neg <- as.data.frame(as.table(rho_t0)) %>%
   rename(
     edge_from = Var1,
     edge_to = Var2,
     rho = Freq
   ) %>%
   mutate(
     pvalue = as.vector(p_t0),
     padj = as.vector(p_adj_t0)
   ) %>%
   filter(
     edge_from != edge_to,
     rho < -0.6,
     padj < 0.05
   ) %>%
   rowwise() %>%
   mutate(
     pair_id = paste(sort(c(edge_from, edge_to)), collapse = "___")
   ) %>%
   ungroup() %>%
   distinct(pair_id, .keep_all = TRUE) %>%
   select(-pair_id)
 
 
 edge_table_t1_neg <- as.data.frame(as.table(rho_t1)) %>%
   rename(
     edge_from = Var1,
     edge_to = Var2,
     rho = Freq
   ) %>%
   mutate(
     pvalue = as.vector(p_t1),
     padj = as.vector(p_adj_t1)
   ) %>%
   filter(
     edge_from != edge_to,
     rho < -0.6,
     padj < 0.05
   ) %>%
   rowwise() %>%
   mutate(
     pair_id = paste(sort(c(edge_from, edge_to)), collapse = "___")
   ) %>%
   ungroup() %>%
   distinct(pair_id, .keep_all = TRUE) %>%
   select(-pair_id)
 
 #Crear grafo T0
 
 g_t0_neg <- graph_from_data_frame(
   d = edge_table_t0_neg,
   directed = FALSE
 )
 
 #Crear grafo T1
 g_t1_neg <- graph_from_data_frame(
   d = edge_table_t1_neg,
   directed = FALSE
 )
 
 #Agregar categoría a los taxa
 
 taxa_category_baseline <- taxa_ecology_t0 %>%
   select(taxon, category) %>%
   rename(category_baseline = category)
 
 V(g_t0_neg)$category_baseline <- taxa_category_baseline$category_baseline[
   match(V(g_t0_neg)$name, taxa_category_baseline$taxon)
 ]
 
 V(g_t1_neg)$category_baseline <- taxa_category_baseline$category_baseline[
   match(V(g_t1_neg)$name, taxa_category_baseline$taxon)
 ]
 
 V(g_t0_neg)$category_baseline[is.na(V(g_t0_neg)$category_baseline)] <- "Not detected in T0"
 V(g_t1_neg)$category_baseline[is.na(V(g_t1_neg)$category_baseline)] <- "Not detected in T0"
 
 #Calcular centralidad negativa
 
 
 V(g_t0_neg)$neg_degree <- degree(g_t0_neg)
 V(g_t1_neg)$neg_degree <- degree(g_t1_neg)
 
 V(g_t0_neg)$neg_strength <- strength(
   g_t0_neg,
   weights = abs(E(g_t0_neg)$rho)
 )
 
 V(g_t1_neg)$neg_strength <- strength(
   g_t1_neg,
   weights = abs(E(g_t1_neg)$rho)
 )
 
 
 #######
 # Netowrks visualization
 #######
 
 category_colors <- c(
   "Core" = "#1b9e77",
   "Intermediate" = "#7570b3",
   "Accessory" = "#d95f02",
   "Not detected in T0" = "grey70"
 )
 
 
 edge_limits <- range(
   c(
     E(g_t0_pos)$rho,
     E(g_t1_pos)$rho,
     abs(E(g_t0_neg)$rho),
     abs(E(g_t1_neg)$rho)
   ),
   na.rm = TRUE
 )
 
 edge_scale_common <- scale_edge_width(
   name = expression("|rho|"),
   range = c(0.2, 1.3),
   limits = edge_limits,
   breaks = c(0.6, 0.7, 0.8, 0.9, 1.0)
 )
 
 
 set.seed(123)
 
 p_net_t0 <- ggraph(g_t0_pos, layout = "fr") +
   geom_edge_link(
     aes(width = rho),
     alpha = 0.25,
     color = "grey20"
   ) +
   geom_node_point(
     aes(
       color = category_baseline
     ),
     alpha = 0.9,
     size=2
   ) +
   scale_color_manual(values = category_colors) +
   edge_scale_common +
   theme_void() +
   labs(
     #  title = "Positive co-occurrence network - T0",
     color = "Prevalence",
     edge_width = "|rho|"
   )
 
 p_net_t0
 
 
 
 set.seed(123)
 
 p_net_t1 <- ggraph(g_t1_pos, layout = "fr") +
   geom_edge_link(
     aes(width = rho),
     alpha = 0.25,
     color = "grey20"
   ) +
   geom_node_point(
     aes(
       color = category_baseline
     ),
     alpha = 0.9,
     size=2
   ) +
   scale_color_manual(values = category_colors) +
   edge_scale_common +
   theme_void() +
   theme(legend.position = "right")+
   labs(
     #  title = "Positive co-occurrence network - T1",
     color = "Prevalence",
     edge_width = "|rho|"
   )+
   guides(
     edge_width = "none", 
     color="none"
   )
 
 p_net_t1
 
 set.seed(123)
 layout_t0_neg <- create_layout(g_t0_neg, layout = "fr")
 
 set.seed(123)
 layout_t1_neg <- create_layout(g_t1_neg, layout = "fr")
 
 p_neg_t0 <- ggraph(layout_t0_neg) +
   geom_edge_link(
     aes(width = abs(rho)),
     alpha = 0.35,
     color = "grey40"
   ) +
   geom_node_point(
     aes(
       #  size = neg_degree,
       color = category_baseline
     ),
     alpha = 0.9
   ) +
   scale_color_manual(values = category_colors) +
   edge_scale_common +
   theme_void() +
   theme(legend.position = "none")+
   labs(
     # title = "Negative co-occurrence network - T0",
     color = "T0 category",
     edge_width = "|rho|"
   )+
   guides(
     edge_width = "none", 
     color="none"
   )
 
 p_neg_t0
 
 p_neg_t1 <- ggraph(layout_t1_neg) +
   geom_edge_link(
     aes(width = abs(rho)),
     alpha = 0.35,
     color = "grey40"
   ) +
   geom_node_point(
     aes(
       # size = neg_degree,
       color = category_baseline
     ),
     alpha = 0.9
   ) +
   scale_color_manual(values = category_colors) +
   edge_scale_common +
   scale_size_continuous(range = c(2, 8)) +
   theme_void() +
   theme(legend.position = "none")+
   labs(
     #title = "Negative co-occurrence network - T1",
     color = "T0 category",
     edge_width = "|rho|"
   )+
   guides(
     edge_width = "none", 
     color="none"
   )
 
  p_neg_t1

 #Boxplot
 
 neg_node_t0 <- data.frame(
   taxon = V(g_t0_neg)$name,
   neg_degree = degree(g_t0_neg),
   Time = "T0"
 )
 
 neg_node_t1 <- data.frame(
   taxon = V(g_t1_neg)$name,
   neg_degree = degree(g_t1_neg),
   Time = "T1"
 )
 
 neg_nodes <- bind_rows(neg_node_t0, neg_node_t1) %>%
   left_join(taxa_category_baseline, by = "taxon") %>%
   mutate(
     category_baseline = ifelse(is.na(category_baseline),
                                "Not detected in T0",
                                category_baseline)
   )
 
 
 # Define factors
 neg_nodes <- neg_nodes %>%
   mutate(
     category_baseline = factor(
       category_baseline,
       levels = c(
         "Core",
         "Intermediate",
         "Accessory",
         "Not detected in T0"
       )
     ),
     Time = factor(
       Time,
       levels = c("T0", "T1")
     )
   )
 
 #statistical differences between negative degrees distributions
 
 wilcox_res_neg <- neg_nodes %>%
   group_by(category_baseline) %>%
   wilcox_test(
     neg_degree ~ Time
   ) %>%
   adjust_pvalue(method = "BH") %>%
   add_significance("p.adj")
 
 wilcox_res_neg
 
 wilcox_res_neg <- wilcox_res_neg %>%
   add_xy_position(
     x = "category_baseline"
   )
 
 fill_cols <- c(
   "T0" = "#b087c4",
   "T1" = "#87afc4"
 )
 
 
 ggplot(
   node_centrality,
   aes(
     x = category_baseline,
     y = degree,
     fill = time
   )
 ) +
   geom_boxplot(
     outlier.shape = NA,
     alpha = 0.6,
     position = position_dodge(width = 0.8)
   ) +
   geom_jitter(
     aes(color = time),
     position = position_jitterdodge(
       jitter.width = 0.15,
       dodge.width = 0.8
     ),
     #alpha = 0.45,
     size = 2
   ) +
   stat_pvalue_manual(
     wilcox_res_degree,
     label = "p.adj",
     tip.length = 0.01
   ) +
   theme_classic() +
   labs(
     x = "Taxa category",
     y = "Degree",
     fill = "Time",
     color = "Time"
   ) +
   scale_fill_manual(values=fill_cols)+  
   scale_color_manual(values=fill_cols) -> pos_boxplot
 
 pos_boxplot
 
 ggplot(
   neg_nodes,
   aes(
     x = category_baseline,
     y = neg_degree,
     fill = Time
   )
 ) +
   geom_boxplot(
     outlier.shape = NA,
     alpha = 0.6
   ) +
   geom_jitter(
     aes(color = Time),
     position = position_jitterdodge(
       jitter.width = 0.15,
       dodge.width = 0.8
     ),
     size = 2
   ) +
   stat_pvalue_manual(
     wilcox_res_neg,
     label = "p.adj",
     tip.length = 0.01
   ) +
   theme_classic() +
   theme(legend.position = "none")+
   labs(
     x = "Taxa category",
     y = "Negative degree",
     fill = "Time",
     color = "Time"
   ) +
   scale_fill_manual(values = fill_cols, guide = "none") +
   scale_color_manual(values = fill_cols, guide = "none") -> neg_boxplot
 

  neg_boxplot 
 
  
  
  final_fig <- 
    (p_net_t0+ p_net_t1) /
    (p_neg_t0 + p_neg_t1) /
    (pos_boxplot + neg_boxplot) +
    plot_layout(
      ncol = 1,
      heights = c(1.3, 1.3, 1),
      guides = "collect"
    ) +
    plot_annotation(
      tag_levels = list(c("A", "B", "C", "D", "E", "F")),
      theme = theme(
        plot.tag = element_text(
          size = 18,
          face = "bold"
        )
      )
    ) &
    theme(
      legend.position = "right"
    )
  
  final_fig
  