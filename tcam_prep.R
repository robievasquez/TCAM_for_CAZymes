# Impute microbiome data
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Load libraries
library(tidyverse)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.5, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

# Function to impute missing values using the average of the week before and after
impute_missing <- function(x) {
  n <- length(x)
  if (n < 3) return(x)  # If there are less than 3 elements, return as is
  
  for (i in 2:(n - 1)) {
    if (is.na(x[i]) & !is.na(x[i - 1]) & !is.na(x[i + 1])) {
      x[i] <- mean(c(x[i - 1], x[i + 1]), na.rm = TRUE)
    } else if (is.na(x[i]) & !is.na(x[i - 1])){
      x[i] <- x[i - 1] + (x[i - 1] - x[i - 2])
      if(x[i] < 0) x[i] <- x[i - 1]
    }
  }
  return(x)
}

# Open data
cayman_t <- rio::import("cayman/tables/families_cpm_table.tsv")
meta <- readRDS("meta_microbiome_run1.RDS")

dim(cayman_t) #307 features, 140 samples

#low count samples to exclude
lowcountids <- meta$ID[which(meta$LowCount == TRUE)]
lowcountids #"L98"  "L110" "L111"

#Tidy cayman data
dim(cayman_t) #307 features, 140 samples
names(cayman_t)
cayman_t <- cayman_t %>%
  dplyr::rename(name = family) %>%
  dplyr::select(-any_of(lowcountids))

names(cayman_t)
dim(cayman_t) # features x samples after removing low-count samples
rownames(cayman_t) <- cayman_t$name
cayman_t$name <- NULL

# Transpose to sample x family table and append metadata
df <- as.data.frame(t(as.matrix(cayman_t)))
df$ID <- rownames(df)
df <- df %>%
  mutate(across(-ID, as.numeric)) %>%
  left_join(meta, by = "ID") %>%
  filter(!is.na(MouseID), !is.na(Age_ints), !is.na(Age_weeks))

head(df)[1:5,1:5]

meta_cols <- c("ID", "MouseID", "Age_weeks", "Age_ints", "Genotype", "Sex", "GenotypePerSex", "LowCount")
# Keep only true family columns (present in the transposed family table), not extra metadata fields from `meta`
microbiome_cols <- intersect(names(df), rownames(cayman_t))
head(microbiome_cols); tail(microbiome_cols)

# Filter out low-prevalence CAZyme families before imputation so prevalence is based on observed data
prevalence_threshold <- 0.10
microbiome_prevalence <- colMeans(df[, microbiome_cols] > 0, na.rm = TRUE)
microbiome_cols <- microbiome_cols[microbiome_prevalence >= prevalence_threshold]
message(length(microbiome_cols), " CAZyme families retained after prevalence filtering (>=10% of samples)")

# Expand to a complete mouse x timepoint grid so missing timepoints become explicit rows
mouse_meta <- df %>% distinct(MouseID, Genotype, Sex, GenotypePerSex)
timepoints <- df %>% distinct(Age_ints, Age_weeks)
df_complete <- mouse_meta %>%
  tidyr::crossing(timepoints) %>%
  left_join(df, by = c("MouseID", "Age_ints", "Age_weeks", "Genotype", "Sex", "GenotypePerSex")) %>%
  arrange(MouseID, Age_ints)

# Pivot to long format for imputation
dftot_long <- df_complete %>%
  pivot_longer(cols = all_of(microbiome_cols), names_to = "microbe", values_to = "abundance") %>%
  arrange(MouseID, Age_ints, microbe) %>%
  mutate(Age_weeks = fct_inorder(as.factor(Age_weeks)), 
         MouseID = fct_inorder(as.factor(MouseID)))

# Apply the imputation function to each microbiome variable
dftot_long_imputed <- dftot_long %>%
  group_by(MouseID, microbe) %>%
  arrange(Age_weeks) %>%
  mutate(abundance = impute_missing(abundance)) %>%
  ungroup()

# And pivot back to wide format
df_imputed_wide <- dftot_long_imputed %>%
  pivot_wider(names_from = microbe, values_from = abundance) %>%
  select(ID, MouseID, Age_weeks, Age_ints, Genotype, Sex, GenotypePerSex, all_of(microbiome_cols)) %>%
  mutate(ID = case_when(is.na(ID) ~ str_c("I", row_number()), 
            .default = ID))

# Filter out terminal timepoints with too many missings (see script 2_1)
df_imputed_wide <- df_imputed_wide %>% filter(!Age_ints %in% c(18, 26))

## This filter is applied last so that the timepoints could still be used to impute neigbouring timepoints ##
any(is.na(df_imputed_wide)) # FALSE - There's no missing values left

# Print the resulting data frame and check for 2 mice
df %>% # BEFORE IMPUTATION, including week 8 and 18
       filter(MouseID == "80") %>% 
  select(all_of(microbiome_cols[15]), Age_weeks) %>% head(n = 30)
df_imputed_wide %>% # AFTER IMPUTATION
       filter(MouseID == "80") %>% 
  select(all_of(microbiome_cols[15]), Age_weeks, GenotypePerSex) %>% 
       head(n = 8)

df %>% # BEFORE IMPUTATION, including week 8 and 18
       filter(MouseID == "79") %>% 
  select(all_of(microbiome_cols[15]), Age_weeks) %>% head(n = 30)
df_imputed_wide %>% # AFTER IMPUTATION
       filter(MouseID == "79") %>% 
  select(all_of(microbiome_cols[15]), Age_weeks, GenotypePerSex) %>% 
       head(n = 8)

# Save the imputed df
write.csv(df_imputed_wide, "new_tcam/imputed_microbiome_data_2.csv", row.names = FALSE)
saveRDS(df_imputed_wide, "new_tcam/imputed_microbiome_data_2.RDS")

# Lineplot before and after imp as example
microbiome_var <- microbiome_cols[1]  # Select the first microbiome variable for plotting
data_before_imputation <- df %>%
  arrange(Age_ints) %>%
  mutate(Age_weeks = fct_inorder(as.factor(Age_weeks))) %>%
  filter(MouseID == "79") %>%
  select(MouseID, Age_weeks, all_of(microbiome_var), Age_ints) %>%
  filter(!Age_ints %in% c(18))
  
data_after_imputation <- df_imputed_wide %>% filter(MouseID == "79") %>% 
       select(MouseID, Age_weeks, any_of(microbiome_var), Age_ints)

(pl1 <- ggplot(data = data_before_imputation, 
       aes(x = Age_ints, y = .data[[microbiome_var]], na.rm = FALSE)) +
  geom_line(aes(group = MouseID)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  labs(title = "Sample 79 before imputation",
       x = "Age (weeks)",
       y = str_c("Relative abundance ", microbiome_var)) +
  theme_Publication())

(pl2 <- ggplot(data = data_after_imputation, 
        aes(x = Age_ints, y = .data[[microbiome_var]], na.rm = FALSE)) +
  geom_line(aes(group = MouseID)) +
  geom_point() +
  labs(title = "Sample 79 after imputation",
       x = "Age (weeks)",
       y = str_c("Relative abundance ", microbiome_var)) +
  theme_Publication())

ggpubr::ggarrange(pl1, pl2, labels = c("A", "B"), ncol = 2, nrow = 1)

# Save the plot
ggsave("new_tcam/sample_line_plot_2.pdf", width = 10, height = 6)
