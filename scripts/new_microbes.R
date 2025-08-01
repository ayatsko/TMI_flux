bac <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/metagenomes/Bact_tax_abun.csv")
sample_metadata <- read.csv("/Users/abbeyyatsko/Downloads/sample_metadata.csv")
seq_metadata <- read.csv("/Users/abbeyyatsko/Downloads/metadata-15165487-processed-ok.csv")
fluxes <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/finalfluxes/all_fluxes_mound.csv")
methanotroph <- read.csv("/Users/abbeyyatsko/Downloads/methanotroph_genera.csv")
flux_relabund_pmoa <- read.csv("/Users/abbeyyatsko/Downloads/flux_relabund.csv")

# subset seq_metadata
seq_metadata_x <- seq_metadata[, c("accession", "sample_name")]
sample_metadata_x <- sample_metadata[, c("sequencing_id", "sample_id", "extract_ID")]
colnames(sample_metadata_x)[2] <- "sample_name"
ids <- left_join(
  seq_metadata_x,
  sample_metadata_x,
  by = "sample_name"
)

# join with flux data 
fluxes <- fluxes[fluxes$campaign == "may22", ]
fluxes$extract_ID <- paste(fluxes$sample,"-",fluxes$tag, "-", fluxes$position, sep = "")
x <- fluxes[, c("ID_measurement", "extract_ID", "flux.CH4", "avg_respT", "species_s")]

y <- left_join(ids, x, by = "extract_ID")

#pivot longer for bac so that all SRR are rows 
bac_long <- bac %>%
  pivot_longer(
    cols = starts_with("SRR"),
    names_to = "accession",
    values_to = "abundance"
  )

out <- left_join(bac_long, y, by = "accession") 

# filter out for Genus in methanotroph 
methanotroph$genus <- paste0("g__", methanotroph$genus)
out_m <- out %>%
  filter(Genus %in% methanotroph$genus) %>%
  filter(!is.na(flux.CH4))

out_m$mound_ID <- sub("-.*", "", out_m$extract_ID)

# nmds 
library(vegan)

# create matrix of data 
out_m %>% dplyr::select(accession, Species, abundance) %>%
  filter(!(accession %in% c("SRR32638406", "SRR32638398")) & !(Species == "s__")) %>% # this one is a clear outlier on nmds plot
  group_by(accession, Species) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(abundance = as.numeric(abundance)) %>%
  pivot_wider(names_from = Species, values_from = abundance, id_cols = accession) %>%
  column_to_rownames("accession") -> Y_matrix

out_m %>% dplyr::select(accession, Genus, abundance) %>%
  filter(!(Genus == "g__")) %>% 
  group_by(accession, Genus) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(abundance = as.numeric(abundance)) %>%
  pivot_wider(names_from = Genus, values_from = abundance, id_cols = accession) %>%
  column_to_rownames("accession") -> Y_matrix

out_m %>% dplyr::select(accession, Species, species_s, abundance, flux.CH4, mound_ID) %>%
  filter(!accession %in% c("SRR32638406", "SRR32638398") & !(Species == "s__")) %>% # this one is a clear outlier on nmds plot
  group_by(accession, Species, species_s, flux.CH4, mound_ID) %>%
  summarise(abundance = sum(abundance)) %>%
  pivot_wider(names_from = Species, values_from = abundance)-> X_matrix

out_m %>% dplyr::select(accession, Genus, species_s, abundance, flux.CH4, mound_ID) %>%
  filter(!(Genus == "g__")) %>% # this one is a clear outlier on nmds plot
  group_by(accession, Genus, species_s, flux.CH4, mound_ID) %>%
  summarise(abundance = sum(abundance)) %>%
  pivot_wider(names_from = Genus, values_from = abundance)-> X_matrix

# Transform data
spe.h <- decostand(Y_matrix, "hellinger")

# Calculate distance matrix
spe_distmat <- vegdist(spe.h, method = "bray") %>%  as.matrix()

# Running NMDS in vegan (metaMDS)
nmds_output <- metaMDS(spe_distmat, distance = "bray", k = 6, maxit = 999, trymax = 500, wascores = T)

# Rownames
accession <- rownames(Y_matrix) 

# Extracting NMDS scores
data.scores <- as.data.frame(scores(nmds_output, display="sites"))  %>% 
  as_tibble() %>% 
  add_column(accession)

# Adding metadata to NMDS scores
nmds_scores_df <- left_join(data.scores, X_matrix, by = "accession") 
nmds_scores_df <- left_join(nmds_scores_df, seq_metadata_x, by = "accession")
nmds_scores_df <- left_join(nmds_scores_df, sample_metadata_x, by = "sample_name")
nmds_scores_df <- left_join(nmds_scores_df, flux_relabund_pmoa)

nmds_scores_df$variable_emitters <- as.factor(ifelse(is.na(nmds_scores_df$variable_emitters), "n", nmds_scores_df$variable_emitters))

ggplot(data=nmds_scores_df) + 
  geom_point(aes(x=NMDS1,y=NMDS2, color = mound_ID, shape = species_s), size = 4)+
  stat_ellipse(aes(x=NMDS1,y=NMDS2, color = mound_ID), alpha = 0.6, level = 0.95) +
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)

mod_adon <- adonis2(Y_matrix ~ species_s + flux.CH4, permutations = 999, 
                    data = X_matrix, 
                    by = "terms")

mod <- betadisper(vegdist(decostand(Y_matrix, 'hellinger'), 'jaccard'),
                  interaction(nmds_scores_df$mound_ID, nmds_scores_df$variable_emitters), bias.adjust=T)
anova(mod)

ggplot(nmds_scores_df, aes(x = NMDS1, y = resampled_relabund)) + 
  geom_point() +
  geom_smooth(method = "loess", se = FALSE)

mod <- lm(resampled_relabund ~ NMDS1 * NMDS2, data = nmds_scores_df)
summary(mod)

ggpredict(mod, terms = c("NMDS1", "NMDS2")) %>%
  plot() +
  labs(x = "NMDS1", y = "Predicted Relabund PMOA") +
  theme_minimal()

# Extract species scores
spp.fit <- envfit(nmds_output, Y_matrix)

# ggplot plot - species scores
spp.scrs <- as.data.frame(scores(spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) 
spp.scrs <- cbind(spp.scrs, pval = spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

nmds_scores_df$NMDS1.2 <- sqrt(nmds_scores_df$NMDS1^2 + nmds_scores_df$NMDS2^2)

# SEM too 
vars_needed <- c("flux.CH4", "species_s", "NMDS1.2", "resampled_relabund")
clean_data <- na.omit(nmds_scores_df[ , vars_needed])

mod1 <- lm(NMDS1.2 ~ species_s, data = clean_data)
mod2 <- lm(resampled_relabund ~ NMDS1.2, data = clean_data)
mod3 <- lm(flux.CH4 ~ NMDS1.2 + species_s + resampled_relabund, data = clean_data)

sem_model <- psem(
  mod1,
  mod2,
  mod3
)

summary(mod3)

# flux ~ relabund 
mod <- lm(flux.CH4 ~ resampled_relabund + NMDS1*NMDS2 + avg_respT + species_s, data = nmds_scores_df)
summary(mod)

ggplot(nmds_scores_df, aes(x = resampled_relabund, y = flux.CH4)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "pmoA", y = "Flux CH4") +
  theme_minimal()


# WAPL 
library(rioja)

model <- WAPLS(Y_matrix, nmds_scores_df$flux.CH4, npls = 5)
summary(model)

model_cv <- crossval(model)
summary(model_cv)

# yj <- yeojohnson(nmds_scores_df$flux.CH4)
# nmds_scores_df$flux_yj <- yj$x.t

performance(model)
rand.t.test(model_cv)

# sum of all methanotrophs per sample
nmds_scores_df <- nmds_scores_df %>%
  mutate(meth_sum = rowSums(across(starts_with("g__")))) %>%
  mutate(shannon = diversity(dplyr::select(., starts_with("g__")), index = "shannon"))

mod <- lm(flux.CH4 ~ shannon, data = nmds_scores_df) # flux.CH4 as predictor
summary(mod)

# flux ~ meth_sum
mod <- lm(flux.CH4 ~ meth_sum, data = nmds_scores_df)
# pmoa ~ meth_sum 
mod <- lm(resampled_relabund_scaled ~ meth_sum, data = nmds_scores_df)
summary(mod)

