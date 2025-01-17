# ko number
integer_cols <- sapply(ko, is.integer)
total_reads <- colSums(ko[, integer_cols])

total_reads_ko <- data.frame(
  Sample = names(total_reads),
  Sum = total_reads,
  row.names = NULL
)

# plot distribution of total reads per sample
ggplot(total_reads_ko, aes(x = reorder(Sample, -Sum), y = Sum)) +
  geom_bar(stat = "identity") +
  labs(title = "KO table - reads per sample", x = "Total reads", y = "Frequency") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ko_ok <- ko[ , !names(ko) %in% c("A108","A109", "A107", "A38", "A40")]
ko_ok <- ko_ok[, -1]

ko_ok_long <- ko_ok %>% 
  pivot_longer(cols = -Annotation, names_to = "sequencing_id", values_to = "reads")

# resample for lower limit  
total_reads_ko <- total_reads_ko[!(total_reads_ko$Sample %in% c("A108","A109", "A107", "A38", "A40")), ]
lower_limit <- min(total_reads_ko$Sum)

# relative abundance 
ko_ok_long <- left_join(ko_ok_long, total_reads_ko, by = c("sequencing_id" = "Sample"))
ko_ok_long$relabund <- ko_ok_long$reads/ko_ok_long$Sum
ko_ok_long %>% group_by(sequencing_id) %>% summarise(Sum_prop = sum(relabund)) -> out

# pmoa gene only - target methanotrophy 
pmoa <- ko_ok_long[grepl("pmoA", ko_ok_long$Annotation), ]

# relative abundance of pmoa across samples
ggplot(pmoa, aes(x = sequencing_id, y = relabund, fill = sequencing_id)) +
  geom_bar(stat = "identity") +
  labs(title = "Relative abundance of methane pathway per sample", x = "Sample", y = "Relative abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# join pmoa relabund with flux data
flux_relabund <- left_join(pmoa, metadata, by = c("sequencing_id"))
flux_relabund$species <- as.factor(flux_relabund$species)

# remove samples with no methane flux data
flux_relabund <- flux_relabund[!is.na(flux_relabund$flux.CH4), ]

# standardise data by rescaling with mean = 0 and sd = 1 
flux_relabund <- flux_relabund %>%
  mutate(resampled_relabund_scaled = standardize(relabund))

flux_relabund$resampled_relabund_scaled <- as.numeric(flux_relabund$resampled_relabund_scaled)

# consider only positive values in the model 
flux_relabund_pos <- flux_relabund[flux_relabund$flux.CH4 > 0, ]

# does greater relative abundance of methane metabolism pathway correlate with lower methane flux?
ggplot(flux_relabund_pos, aes(x = resampled_relabund_scaled, y = log(flux.CH4))) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_classic() + 
  xlab("Relative abundance of pmoA (mean centred)") + 
  ylab("CH4 Flux")

flux_relabund_pos$species_s <- as.factor(flux_relabund_pos$species_s)
mod <- lmer(log(flux.CH4) ~ resampled_relabund_scaled + (1|sample), data = flux_relabund_pos)
summary(mod)
Anova(mod)

# species difference in relative abundance
ggplot(flux_relabund_pos, aes(x = species_s, y = relabund)) +
  geom_boxplot() +
  theme_classic() +
  xlab("") + 
  ylab("Relative abundance of pmoA")

mod <- aov(relabund ~ species_s, data = flux_relabund_pos)
summary(mod)
emmeans(mod, pairwise ~ species_s)



