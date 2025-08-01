## Libraries
library(dplyr)
library(glmmTMB)
library(ggplot2)
library(ggeffects)
library(DHARMa)
library(bestNormalize)
library(lme4)
library(car)
library(emmeans)
library(MuMIn)

## Data
p_fluxes_f <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/finalfluxes/p_fluxes_f.csv")
all_fluxes_mound_avg <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/finalfluxes/all_fluxes_mound_avg.csv")
all_fluxes_mound <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/finalfluxes/all_fluxes_mound.csv")
all_fluxes_resample_reps_avg <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/finalfluxes/all_fluxes_resample_reps_avg.csv")
flux_relabund_pmoa <- read.csv("/Users/abbeyyatsko/Downloads/flux_relabund.csv")

# Result 1: CH4 emissions at the individual termite and mound-level
## individual 
p_fluxes_f$species <- as.factor(p_fluxes_f$species)
mod <- aov(TEF_ug_g_h ~ species, data = p_fluxes_f)
summary(mod) # posthoc shows significant differences between all 
emmeans(mod, pairwise ~ species, adjust = "tukey") 

## mound 
all_fluxes_mound$ID_mound <- as.factor(all_fluxes_mound$ID_mound)
all_fluxes_mound$species_s <- as.factor(all_fluxes_mound$species_s)
mod <- lmer(flux.CH4 ~ species_s + (1 | ID_mound), data = all_fluxes_mound)
mod <- lmer(mean_mound_fluxCH4 ~ species_s + (1 | sample), data = all_fluxes_resample_reps_avg)
Anova(mod) # posthoc shows N. magnus is significantly higher than other two species

# Result 2: Relationship between mound microbial communities and CH4 emission
# for POSITIVE ONLY fluxes: flux_relabund_pos <- flux_relabund_pmoa[flux_relabund_pmoa$flux.CH4 > 0, ]
mod <- lmer(log(flux.CH4) ~ resampled_relabund_scaled + species_s + (1|sample), data = flux_relabund_pos)
Anova(mod)

mod <- glmmTMB(flux.CH4 ~ resampled_relabund_scaled + species_s  + (1 | sample),
        data = flux_relabund_pmoa,
        family = tweedie())
summary(mod)

# for ALL fluxes: 
# Apply Yeo-Johnson transformation - allows for zero and negative values
yj <- yeojohnson(flux_relabund_pmoa$flux.CH4)
flux_relabund_pmoa$flux_yj <- yj$x.t

# Fit the model using transformed response
mod_yj <- lmer(flux.CH4 ~ resampled_relabund_scaled + species_s + (1 | sample), data = flux_relabund_pmoa)
DHARMa::simulateResiduals(mod_yj, plot = TRUE)
Anova(mod_yj)

# nonlinear doesnt really make model checks a bit better
mod2 <- glmmTMB(flux_yj ~ poly(resampled_relabund_scaled, 2) + species_s + (1 | sample),
                data = flux_relabund_pmoa,
                family = gaussian())

# Result 3: Influence of mound wall thickness on CH4 emission
wt <- all_fluxes_mound[!is.na(all_fluxes_mound$wall_thickness_mm), ]
wt$sample <- as.factor(wt$sample)
wt$wall_thickness_centered <- ave(wt$wall_thickness_mm, wt$species_s, FUN = scale)

# for POSITIVE ONLY fluxes: wt_pos <- wt[wt$flux.CH4 > 0, ]
mod <- lmer(log(flux.CH4) ~ wall_thickness_centered*species_s + (1 | ID_mound), data = wt_pos)
Anova(mod)

# for ALL fluxes: 
yj <- yeojohnson(wt$flux.CH4)
wt$flux_yj <- yj$x.t
mod_yj <- glmmTMB(flux_yj ~ wall_thickness_centered + species_s + (1 | sample),
               data = wt,
               family = gaussian())
DHARMa::simulateResiduals(mod_yj, plot = TRUE)
summary(mod_yj)

# viz for wt
preds <- ggpredict(mod_yj, terms = c("wall_thickness_centered", "species_s"))
ggplot() +
  geom_point(data = wt, aes(x = wall_thickness_centered, y = flux_yj, color = species_s), alpha = 0.6) +
  geom_line(data = preds, aes(x = x, y = predicted, color = group), size = 1.2) +
  geom_ribbon(data = preds, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(x = "Wall Thickness (centered)", y = "CH4 Flux", color = "Species", fill = "Species",
       title = "Predicted CH4 Flux by Wall Thickness and Species") +
  theme_minimal()

# Result 4: Influence of mound volume on CH4 emission
all_fluxes_mound_avg$volume_m3_final_centered <- ave(all_fluxes_mound_avg$volume_m3_final, all_fluxes_mound_avg$species_s, FUN = scale)
all_fluxes_mound_avg <- all_fluxes_mound_avg[all_fluxes_mound_avg$ID_mound != "MD8-may22", ] # average flux is negative

# for ALL fluxes: 
mod <- lmer(mound_CH4flux ~ volume_m3_final_centered*species_s + (1 | sample), data = all_fluxes_mound_avg)
Anova(mod)

yj <- yeojohnson(all_fluxes_mound_avg$mound_CH4flux)
all_fluxes_mound_avg$flux_yj <- yj$x.t
mod <- glmmTMB(flux_yj ~ volume_m3_final_centered*species_s + (1 | sample),
        data = all_fluxes_mound_avg,
        family = gaussian())
summary(mod)
DHARMa::simulateResiduals(mod, plot = TRUE)

# nonlinear makes the model checks a bit better
mod2 <- glmmTMB(flux_yj ~ poly(volume_m3_final_centered, 2)*species_s + (1 | sample), 
                data = all_fluxes_mound_avg,
                family = gaussian())

# Result 5: Seasonal change in mound CH4 emission and relationship with temperature
mod <- lmer(mean_mound_fluxCH4 ~ temp * campaign + (1 | sample) + (1 | species_s), data = all_fluxes_resample_reps_avg) # model does better when species is taken out as random effect
# mod <- lmer(mean_mound_fluxCH4 ~ temp + campaign*species_s + (1 | sample) , data = all_fluxes_resample_reps_avg)
Anova(mod)
emm <- emmeans(mod, ~ campaign)
pairs(emm)

ggplot(all_fluxes_resample_reps_avg, 
       aes(x = campaign, y = mean_mound_fluxCH4, color = species_s, group = species_s)) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(width = 0.2)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.2)) +
  theme_minimal() +
  labs(y = "Mean CH4 flux", x = "Campaign", color = "Species")

# BUILD GLOBAL MODEL (and explore model selection)
# add average pmoA abundance to all_fluxes_resample_reps_avg
# summarise counts by flux_relabund_pmoa$sample
flux_relabund_pmoa$sample <- as.factor(flux_relabund_pmoa$sample)
flux_relabund_pmoa_avg <- flux_relabund_pmoa %>%
  group_by(sample) %>%
  summarise(#count = n(),
            avg_relabund = mean(resampled_relabund, , na.rm = TRUE))

# merge with all_fluxes_resample_reps_avg
all_fluxes_resample_reps_avg <- merge(all_fluxes_resample_reps_avg, flux_relabund_pmoa_avg, by = "sample", all.x = TRUE)

# extract wt and sample when wt is not na
wt_avg <- all_fluxes_resample_reps_avg[!is.na(all_fluxes_resample_reps_avg$wall_thickness_mm), ]
wt_avg <- wt_avg[, c("sample", "wall_thickness_mm")]
colnames(wt_avg) <- c("sample", "wall_thickness_mm_avg")
all_fluxes_resample_reps_avg <- left_join(all_fluxes_resample_reps_avg, wt_avg, by = "sample")

# variable structuring
all_fluxes_resample_reps_avg$species_s <- as.factor(all_fluxes_resample_reps_avg$species_s)
all_fluxes_resample_reps_avg$sample <- as.factor(all_fluxes_resample_reps_avg$sample)
all_fluxes_resample_reps_avg$campaign <- as.factor(all_fluxes_resample_reps_avg$campaign)

# center all species-level variables
all_fluxes_resample_reps_avg$volume_m3_final_scaled <- ave(all_fluxes_resample_reps_avg$volume_m3_final, all_fluxes_resample_reps_avg$species_s, FUN = scale)
all_fluxes_resample_reps_avg$wall_thickness_mm_avg_scaled <- ave(all_fluxes_resample_reps_avg$wall_thickness_mm_avg, all_fluxes_resample_reps_avg$species_s, FUN = scale)
all_fluxes_resample_reps_avg$avg_relabund_scaled <- ave(all_fluxes_resample_reps_avg$avg_relabund, all_fluxes_resample_reps_avg$species_s, FUN = scale)

# global model (without avg_relabund)
mod <- lmer(mean_mound_fluxCH4 ~ species_s + campaign + temp + # species_s*campaign
              volume_m3_final_scaled + wall_thickness_mm_avg_scaled + 
              (1 | sample), data = all_fluxes_resample_reps_avg)

# global model (with avg_relabund)
all_fluxes_resample_reps_avg_relabund <- all_fluxes_resample_reps_avg[!is.na(all_fluxes_resample_reps_avg$avg_relabund), ]
mod <- lmer(mean_mound_fluxCH4 ~ species_s * campaign + temp + 
              volume_m3_final_scaled + wall_thickness_mm_avg_scaled + avg_relabund_scaled +
              (1 | sample), data = all_fluxes_resample_reps_avg_relabund)

Anova(mod)
DHARMa::simulateResiduals(mod, plot = TRUE)

all_fluxes_resample_reps_avg_relabund %>%
  group_by(species_s) %>%
  summarise(count = n_distinct(sample))

# model selection using remove one variable at a time 
mod_minus_species <- update(mod, . ~ . - species_s)
mod_minus_campaign <- update(mod, . ~ . - campaign)
mod_minus_temp <- update(mod, . ~ . - temp)
mod_minus_vol <- update(mod, . ~ . - volume_m3_final_scaled)
mod_minus_wt <- update(mod, . ~ . - wall_thickness_mm_avg_scaled)
# mod_minus_pmoa <- update(mod, . ~ . - avg_relabund_scaled)

anova(mod, mod_minus_species)
anova(mod, mod_minus_campaign)
anova(mod, mod_minus_temp)
anova(mod, mod_minus_vol) # potentially remove
anova(mod, mod_minus_wt) # potentially remove
anova(mod, mod_minus_pmoa) # potentially remove

# model selection using AIC
# Allow lmer models to be used in dredge
options(na.action = "na.fail")

# Run all-subsets model selection
model_set <- dredge(mod)

# View top models
head(model_set)

# SEM attempt  
library(piecewiseSEM)
library(lme4) 

vars_needed <- c("mean_mound_fluxCH4", "species_s", "volume_m3_final_scaled",
                 "wall_thickness_mm_avg_scaled", "temp", "campaign", "sample")
clean_data <- na.omit(all_fluxes_resample_reps_avg[ , vars_needed])

mod1 <- lm(volume_m3_final_scaled ~ species_s, data = clean_data)
mod2 <- lm(wall_thickness_mm_avg_scaled ~ species_s, data = clean_data)
mod3 <- lmer(mean_mound_fluxCH4 ~ volume_m3_final_scaled + 
               wall_thickness_mm_avg_scaled + species_s +
               temp + campaign + (1 | sample), data = clean_data)

sem_model <- psem(
  mod1,
  mod2,
  mod3
)

summary(sem_model) 


## SEM with pmoa 
vars_needed <- c("mean_mound_fluxCH4", "species_s", "volume_m3_final_scaled",
                 "wall_thickness_mm_avg_scaled", "avg_relabund_scaled", "temp", "campaign", "sample")
clean_data <- na.omit(all_fluxes_resample_reps_avg_relabund[ , vars_needed])

mod1 <- lm(volume_m3_final_scaled ~ species_s, data = clean_data)
mod2 <- lm(wall_thickness_mm_avg_scaled ~ species_s, data = clean_data)
mod3 <- lmer(mean_mound_fluxCH4 ~ volume_m3_final_scaled + 
               wall_thickness_mm_avg_scaled + avg_relabund_scaled + species_s +
               temp + campaign + (1 | sample), data = clean_data)
mod4 <- lm(avg_relabund_scaled ~ wall_thickness_mm_avg_scaled + volume_m3_final_scaled, data = clean_data)

sem_model <- psem(
  mod1,
  mod2,
  mod3, 
  mod4
)
summary(sem_model)

# average mound flux by species 
avg_flux_by_species <- all_fluxes_mound_avg %>%
  group_by(species_s) %>%
  summarise(avg_mound_fluxCH4 = mean(mean_mound_fluxCH4, na.rm = TRUE),
            sd_mound_fluxCH4 = sd(mound_CH4flux, na.rm = TRUE),
            n_samples = n_distinct(sample))
