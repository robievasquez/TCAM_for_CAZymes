# TCAM plots ALS project
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ggsci)

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


plot_loadings <- function(data, title) {
  data$comp <- data[[2]]
  ggplot(data, aes(x = reorder(X, comp), y = comp)) +
    geom_bar(stat = "identity", 
        fill = c(rep("firebrick", 10), rep("royalblue",10))) +
    coord_flip() +
    labs(title = title, x = "", y = "Loading") +
    theme_Publication()
}

## Load data
tcam <- read.csv("pythonoutput_2/df_plot.csv") %>% 
  select(1:17) %>% 
  mutate(Genotype = fct_relevel(Genotype, "WT", after = 0L),
          GenotypePerSex = fct_relevel(GenotypePerSex, "Female WT", after = 0L),
          GenotypePerSex = fct_relevel(GenotypePerSex, "Male WT", after = 1L))
head(tcam)
load <- read.csv("pythonoutput_2/df_loadings.csv")
head(load)
mb <- readRDS("imputed_microbiome_data_2.RDS")
meta <- readRDS("meta_microbiome.RDS")

## PERMANOVAs ##
tcam <- tcam %>% filter(Age_ints == 6) 
f1n <- names(tcam)[7]
f2n <- names(tcam)[8]

set.seed(1234)
# 1. PERMANOVA for overall genotype effect
tcam_dist <- tcam %>% select(all_of(c(f1n, f2n))) %>% dist(method = "euclidean")
(permanova_genotype <- adonis2(tcam_dist ~ Genotype, data = tcam, permutations = 999, method = "euclidean", by = "term"))

# 2. PERMANOVA for genotype in females
tcam_fem <- tcam %>% filter(Sex == "Female")
tcam_fem_dist <- tcam_fem %>% select(all_of(c(f1n, f2n))) %>% dist(method = "euclidean")
(permanova_fem <- adonis2(tcam_fem_dist ~ Genotype, data = tcam_fem, permutations = 999, method = "euclidean", by = "term"))

# 3. PERMANOVA for genotype in males
tcam_male <- tcam %>% filter(Sex == "Male")
tcam_male_dist <- tcam_male %>% select(all_of(c(f1n, f2n))) %>% dist(method = "euclidean")
(permanova_male <- adonis2(tcam_male_dist ~ Genotype, data = tcam_male, permutations = 999, method = "euclidean", by = "term"))

# 4. PERMANOVA for genotype in males
tcam_dist <- tcam %>% select(all_of(c(f1n, f2n))) %>% dist(method = "euclidean")
(permanova_interaction <- adonis2(tcam_dist ~ Genotype * Sex, data = tcam, permutations = 999, method = "euclidean", by = "term"))

annotation_text <- str_c("genotype*sex: p = ", permanova_interaction$`Pr(>F)`[3], "\n",
                    "genotype (F): p = ", permanova_fem$`Pr(>F)`[1], "\n",
                    "genotype (M): p = ", permanova_male$`Pr(>F)`[1])

### TCAM plot ###
(tcampl <- ggplot(data = tcam, aes(x = .data[[f1n]], y = .data[[f2n]], 
                        color = Genotype, fill = Genotype, shape = Sex)) +
        #stat_ellipse(geom = "polygon", alpha = 0.3) +
        geom_point(size = 5) +
        ggtitle('TCAM') +
        scale_color_manual(values = pal_nejm()(6)[c(6,3)]) +
        scale_fill_manual(values = pal_nejm()(6)[c(6,3)]) +
        theme_minimal() +
        labs(x=str_c(str_replace(f1n, "[.]", " "), "%"),
            y=str_c(str_replace(f2n, "[.]", " "), "%")) +
        annotate("text", x = Inf, y = Inf, 
                label = annotation_text,
                hjust = 1, vjust = 1, 
                size = 4) +
        theme_Publication() +
        theme(legend.title = element_blank()))
ggsave("routput_2/f1f2_scatter.pdf", width = 6, height = 6) 


### PERMANOVA ###
# 1. PERMANOVA for Sex effect within WT
tcam_wt <- tcam %>% filter(Genotype == "WT")
tcam_wt_dist <- tcam_wt %>% select(all_of(c(f1n, f2n))) %>%dist(method = "euclidean")
(permanova_sex_wt <- adonis2(tcam_wt_dist ~ Sex, data = tcam_wt, permutations = 999, method = "euclidean", by = "term"))

# 2. PERMANOVA for Sex effect within TDP-43
tcam_tdp <- tcam %>% filter(Genotype == "TDP43") 
tcam_tdp_dist <- tcam_tdp %>% select(all_of(c(f1n, f2n))) %>% dist(method = "euclidean")
(permanova_sex_tdp <- adonis2(tcam_tdp_dist ~ Sex, data = tcam_tdp, permutations = 999, method = "euclidean", by = "term"))

pvalues <- c(permanova_sex_wt$`Pr(>F)`[1], permanova_sex_tdp$`Pr(>F)`[1])
genotype <- c("WT", "TDP43")
res <- data.frame(Genotype = genotype, pvalue = pvalues)

## TCAM plot ##
(tcamgenotype <- ggplot(data = tcam, aes(x = .data[[f1n]], y = .data[[f2n]], 
                        color = Sex, fill = Sex, shape = Sex)) +
        stat_ellipse(geom = "polygon", alpha = 0.3) +
        geom_point(size = 3, alpha = 0.75) +
        ggtitle('TCAM: sex differences') +
        scale_color_manual(values = pal_nejm()(2)) +
        scale_fill_manual(values = pal_nejm()(2)) +
        theme_minimal() +
        labs(x=str_c(str_replace(f1n, "[.]", " "), "%"),
            y=str_c(str_replace(f2n, "[.]", " "), "%")) +
        geom_text(data = res, aes(x = Inf, y = Inf, label = paste0("sex: p = ", round(pvalue, 3))),
                              hjust = 1.1, vjust = 1.1, size = 4, inherit.aes = FALSE) +
        theme_Publication() +
        facet_wrap(~Genotype) +
        theme(legend.title = element_blank()))
ggsave(tcamgenotype, filename = "routput_2/f1f2_scatter_pergenotype.pdf", width = 6, height = 6)   

## Loading of component plots ##
head(load)
names(load)[1:10]
last <- nrow(load)
min10 <- last - 9
f1 <- load %>% select(X, all_of(f1n)) %>% arrange(-.data[[f1n]])
f1 <- f1[c(1:10, min10:last),]
f2 <- load %>% select(X, all_of(f2n)) %>% arrange(-.data[[f2n]])
f2 <- f2[c(1:10, min10:last),]

(f1plot <- plot_loadings(f1, "TCAM F1"))
ggsave("routput_2/loading_pc1.pdf", width = 6, height = 7)

plot_loadings(f2, "TCAM F2")
ggsave("routput_2/loading_pc2.pdf", width = 6, height = 7)

## Lineplots ##
f1 <- load %>% select(X, all_of(f1n)) %>% arrange(-.data[[f1n]])
f1 <- f1[c(1:10, (last-9):last),]
f1

pseudocounts <- mb %>%  select(any_of(f1$X)) %>%
  summarise(across(everything(), ~ min(.x[.x > 0], na.rm = TRUE) / 2))
mbsel <- mb %>% ungroup(.) %>% select(any_of(f1$X), ID, Genotype, Age_ints, Sex, GenotypePerSex) %>%
  mutate( across(f1$X, ~log10(.x + pseudocounts[[cur_column()]]))
          ) %>%
  filter(Age_ints <= 16) %>%
  rename_at(c(f1$X), ~str_remove(.x, " \\(.*\\)$"))

colnames(mbsel)

microbe_cols <- setdiff(names(mbsel), c("ID", "Genotype", "Age_ints", "Sex", "GenotypePerSex"))
plot_microbe_cols <- microbe_cols[vapply(mbsel[microbe_cols], function(x) any(!is.na(x)), logical(1))]

means_ses <- mbsel %>%
  group_by(Genotype, Age_ints) %>%
  summarise(across(all_of(plot_microbe_cols),
          list(mean = ~mean(.x, na.rm = TRUE),
            se = ~sd(.x, na.rm = TRUE) / sqrt(n())),
          .names = "{.col}__{.fn}"),
            .groups = "drop")

means_ses_long <- means_ses %>%
  filter(Age_ints %in% c(6, 8, 10, 12, 14, 16)) %>%
  pivot_longer(
    cols = -c(Genotype, Age_ints),
    names_to = c("microbe", ".value"),
    names_sep = "__"
  ) %>%
  mutate(microbe_label = str_remove(microbe, " \\(.*\\)$")) %>%
  filter(!is.na(mean))

mbsel_long <- mbsel %>%
  pivot_longer(
    cols = c(-ID, -Genotype, -Age_ints, -GenotypePerSex, -Sex),
    names_to = "microbe",
    values_to = "value"
  ) %>%
  mutate(microbe_label = str_remove(microbe, " \\(.*\\)$")) %>%
  filter(!is.na(value))

unique(mbsel_long$microbe)

genotype_lineplots <- ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                                    group = Genotype, color = Genotype)) +
                      geom_line() +
                      geom_point(size = 0.5) +
                      geom_jitter(data = mbsel_long, aes(x = Age_ints, y = value, color = Genotype),
                                    width = 0.1, alpha = 0.7) +
                      scale_color_manual(values = pal_nejm()(8)[c(3,6)]) +
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                      labs(title = "Abundance over time - TDP43 vs WT",
                          x = "Age (weeks)",
                          y = "log10(CPM + 1)") +
                      scale_x_continuous(breaks = c(6,8, 10,12,14,16)) +
                      facet_wrap(~ microbe_label, scales = "free", nrow = 5, drop = TRUE) +
                      theme_Publication() +
                      theme(strip.text = element_text(size = 8))
genotype_lineplots

ggsave("routput_2/lineplots_genotype_log10.pdf", width = 15, height = 15)

# The same but then for sex differences within TDP43
means_ses <- mbsel %>%
  filter(Genotype == "TDP43", Age_ints %in% c(6, 8, 10, 12, 14, 16)) %>%
  group_by(Sex, Age_ints) %>%
  summarise(across(all_of(plot_microbe_cols),
          list(mean = ~mean(.x, na.rm = TRUE),
            se = ~sd(.x, na.rm = TRUE) / sqrt(n())),
          .names = "{.col}__{.fn}"),
      .groups = "drop")

means_ses

means_ses_long <- means_ses %>% # to long
  pivot_longer(cols = c(-Sex, -Age_ints), 
                        names_to = c("microbe", ".value"), names_sep = "__") %>%
  mutate(microbe_label = str_remove(microbe, " \\(.*\\)$")) %>%
  filter(!is.na(mean))

mbsel_long <- mbsel %>%
  filter(Age_ints %in% c(6,8, 10,12,14,16)) %>% # only plot up to week 16, since week 18 has many missing values
  pivot_longer(cols = c(-ID, -Genotype, -Age_ints, -GenotypePerSex, -Sex), 
                        names_to = "microbe", values_to = "value") %>%
  mutate(microbe_label = str_remove(microbe, " \\(.*\\)$")) %>%
  filter(!is.na(value))

(sexdiff_lineplots <- ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                                    group = Sex, color = Sex)) +
                      geom_line() +
                      geom_jitter(data = mbsel_long %>% filter(Genotype == "TDP43"), aes(x = Age_ints, y = value, 
                                    color = Sex),
                                    width = 0.1, alpha = 0.7) +
                      scale_color_nejm() +
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                      labs(title = "Abundance over time for TDP43 - males vs females",
                          x = "Age (weeks)",
                          y = "log10(CPM + 1)") +
                      scale_x_continuous(breaks = c(6,8, 10,12,14,16)) +
                      facet_wrap(~ microbe_label, scales = "free", nrow = 5, drop = TRUE)) +
                      theme_Publication()

ggsave("routput_2/lineplots_tdp43_sex.pdf", width = 15, height = 15)

# The same but then for sex differences within WT
means_ses <- mbsel %>% filter(Genotype == "WT", Age_ints %in% c(6,8, 10,12,14,16)) %>% ## Calculate means and standard deviations
  group_by(Sex, Age_ints) %>%
  summarise(across(all_of(plot_microbe_cols), list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~(sd(.x, na.rm = TRUE) / sqrt(n()))), .names = "{.col}__{.fn}"),
            .groups = "drop")
means_ses

means_ses_long <- means_ses %>% # to long
  pivot_longer(cols = c(-Sex, -Age_ints), 
                        names_to = c("microbe", ".value"), names_sep = "__") %>%
  mutate(microbe_label = str_remove(microbe, " \\(.*\\)$")) %>%
  filter(!is.na(mean))
mbsel_long <- mbsel %>%
  pivot_longer(cols = c(-ID, -Genotype, -Age_ints, -GenotypePerSex, -Sex), 
                        names_to = "microbe", values_to = "value") %>%
  mutate(microbe_label = str_remove(microbe, " \\(.*\\)$")) %>%
  filter(!is.na(value))

sexdiff_lineplots_wt <- ggplot(means_ses_long, aes(x = Age_ints, y = mean, 
                                    group = Sex, color = Sex)) +
                      geom_line() +
                      geom_jitter(data = mbsel_long %>% filter(Genotype == "WT"), aes(x = Age_ints, y = value, 
                                    color = Sex),
                                    width = 0.1, alpha = 0.7) +
                      scale_color_nejm() +
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                      labs(title = "Abundance over time for WT - males vs females",
                          x = "Age (weeks)",
                          y = "log10(relative abundance)") +
                      scale_x_continuous(breaks = c(6,8, 10,12,14,16)) +
                      facet_wrap(~ microbe_label, scales = "free", nrow = 5, drop = TRUE) +
                      theme_Publication() +
                      theme(strip.text = element_text(size = 8))

sexdiff_lineplots_wt

ggsave("routput_2/lineplots_wt_sex.pdf", width = 15, height = 15)
