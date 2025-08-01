# create X and Y matrix of data (all samples)
out_m_hl %>% dplyr::select(accession, Genus, abundance) %>%
  filter(!(accession %in% c("SRR32638406")) & !(Genus == "g__")) %>% # outlier based on original nmds plot
  group_by(accession, Genus) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(abundance = as.numeric(abundance)) %>%
  pivot_wider(names_from = Genus, values_from = abundance, id_cols = accession) %>%
  column_to_rownames("accession") -> Y_matrix
out_m_hl %>% dplyr::select(accession, Genus, species_s, abundance, flux.CH4, mound_ID) %>%
  filter(!(accession %in% c("SRR32638406")) & !(Genus == "g__")) %>% # outlier based on original nmds plot
  group_by(accession, Genus, species_s, flux.CH4, mound_ID) %>%
  summarise(abundance = sum(abundance)) %>%
  pivot_wider(names_from = Genus, values_from = abundance)-> X_matrix

# Transform data
spe.h <- decostand(Y_matrix, "hellinger")

# Calculate distance matrix (bray curtis)
spe_distmat <- vegdist(spe.h, method = "bray")

# run PCoA
pcoa_output <- vegan::wcmdscale(spe_distmat, k = 2, eig = TRUE, add = TRUE)

# Extract PCoA scores
pcoa_scores <- as.data.frame(pcoa_output$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
pcoa_scores$accession <- rownames(pcoa_scores)

# join with metadata
merged_df <- merge(pcoa_scores, X_matrix, by = "accession")

# test community differences by species and flux 
mod <- adonis2(spe_distmat ~ species_s + flux.CH4, data = merged_df, permutations = 999)

# dispersion checks 
# For species
bd_species <- betadisper((spe_distmat), merged_df$species_s)
anova(bd_species) #ns

# For flux groups, you might bin flux to check
merged_df$flux_bin <- cut(merged_df$flux.CH4, breaks = 5)  # or quantile-based
bd_flux <- betadisper((spe_distmat), merged_df$flux_bin)
anova(bd_flux) #ns

# genera vectors
fit <- envfit(pcoa_output, spe.h, permutations = 999)
vectors <- as.data.frame(scores(fit, display = "vectors"))
vectors$Genus <- rownames(vectors)
vectors$Dim1 <- vectors$Dim1 * 0.1
vectors$Dim2 <- vectors$Dim2 * 0.1

# visualise 
a <- ggplot(merged_df, aes(x = PCoA1, y = PCoA2, color = species_s)) +
  geom_point(size = 3, aes(color = species_s)) +
  stat_ellipse(type = "t", level = 0.95, linetype = "dashed", aes(color = species_s)) +
  
  # Add genus vectors
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
               arrow = arrow(length = unit(0.25, "cm")),
               inherit.aes = FALSE, color = "black") +
  
  # Add genus labels
  geom_text(data = vectors, aes(x = Dim1, y = Dim2, label = Genus),
            inherit.aes = FALSE, color = "black", size = 3, vjust = 1.2) +
  
  labs(
    x = paste0("PCoA1 (", round(pcoa_output$eig[1] / sum(pcoa_output$eig) * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_output$eig[2] / sum(pcoa_output$eig) * 100, 1), "%)")
  ) +
  scale_color_viridis_d(option = "D") +
  theme_classic()

# Adding metadata to pcoa scores
pcoa_scores_df <- left_join(merged_df, seq_metadata_x, by = "accession")
pcoa_scores_df <- left_join(pcoa_scores_df, sample_metadata_x, by = "sample_name")
pcoa_scores_df <- left_join(pcoa_scores_df, flux_relabund_pmoa)
pcoa_scores_df$variable_emitters <- as.factor(ifelse(is.na(pcoa_scores_df$variable_emitters), "n", pcoa_scores_df$variable_emitters))

# sum of all methanotrophs per sample
pcoa_scores_df <- pcoa_scores_df %>%
  mutate(meth_sum = rowSums(across(starts_with("g__")))) %>%
  mutate(shannon = diversity(dplyr::select(., starts_with("g__")), index = "shannon"))

# calculate relative abundance for all columns beginning with g__ and total methanotroph count
pcoa_scores_df <- left_join(pcoa_scores_df, total_reads)

pcoa_scores_df <- pcoa_scores_df %>%
  mutate(across(starts_with("g__"), ~ . / total_reads, .names = "relabund_{col}")) %>%
  mutate(relabund_meth_sum = (across(starts_with("meth_sum"))[[1]]) / total_reads)

# remove if there are NAs in pmoa
pcoa_scores_df <- pcoa_scores_df %>%
  filter(!is.na(resampled_relabund))

# TEST PCoA axes influence on relabund PMOA
mod_int <- lmer(resampled_relabund_scaled ~ PCoA1 * PCoA2 + (1|mound_ID), data = pcoa_scores_df)
summary(mod_int)

b_pred <- ggpredict(mod_int, terms = c("PCoA1", "PCoA2"))  %>% plot()

sm <- "#fce725"
med <- '#20908c'
lg <- '#440153'
v <- c(-0.02, 0, 0.02)

b <- ggplot(pcoa_scores_df, aes(x = PCoA1, y = resampled_relabund_scaled, color = PCoA2)) +
  geom_jitter() +
  
  # -0.02 group
  geom_ribbon(data = filter(b_pred, group == '-0.02'),
              aes(x = x, ymin = conf.low, ymax = conf.high),
              inherit.aes = FALSE, fill = sm, alpha = 0.2) +
  geom_line(data = filter(b_pred, group == '-0.02'),
            aes(x = x, y = predicted),
            inherit.aes = FALSE, colour = sm) +
  
  # 0 group
  geom_ribbon(data = filter(b_pred, group == '0'),
              aes(x = x, ymin = conf.low, ymax = conf.high),
              inherit.aes = FALSE, fill = med, alpha = 0.2) +
  geom_line(data = filter(b_pred, group == '0'),
            aes(x = x, y = predicted),
            inherit.aes = FALSE, colour = med) +
  
  # 0.02 group
  geom_ribbon(data = filter(b_pred, group == '0.02'),
              aes(x = x, ymin = conf.low, ymax = conf.high),
              inherit.aes = FALSE, fill = lg, alpha = 0.2) +
  geom_line(data = filter(b_pred, group == '0.02'),
            aes(x = x, y = predicted),
            inherit.aes = FALSE, colour = lg) +
  
  scale_color_viridis_c(option = "D") +
  theme_classic()

ggarrange(a, b)

# TEST by individual genus (relative abundance)
# list of unique things begining with "relabund"
relabund_cols <- grep("^relabund_g", names(pcoa_scores_df), value = TRUE)

# loop through each relative abundance column and run a model
for (x in relabund_cols) {
  cat("Model for", x, ":\n")
  mod <- lmer(flux.CH4 ~ get(x) + species_s + (1|mound_ID), data = pcoa_scores_df)
  print(summary(mod))
}
# TEST effect of pmoA abundance on CH4 emission (across mounds)
pcoa_scores_df$mound_ID <- as.factor(pcoa_scores_df$mound_ID)
mod <- lmer(flux.CH4 ~ resampled_relabund_scaled + species_s +  PCoA1 * PCoA2 + relabund_meth_sum +  (1|mound_ID), data = pcoa_scores_df)
summary(mod)

# TEST for variable mounds 
# test community composition, pmoa, total methanotroph abundance, species on flux (response) 
mod <- lmer(flux.CH4 ~ resampled_relabund_scaled + species_s + PCoA1 * PCoA2  + relabund_meth_sum + (1|mound_ID), data = pcoa_scores_df)
summary(mod)
