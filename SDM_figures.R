# This script is designed to provide the figures for the article and 
# includes an ANOVA for the local performance of SDMs.
# WD, packages:
setwd("Z:/jseuren/Stageonderzoek")
library(reshape2)
library(dplyr)
library(ggplot2)
library(stats)
library(patchwork)
library(ggtext)
# Loading in data in the performance of SDMs at local scale:
confusion_matrix_data <- read.csv("./output/modelprestatie/250m/data_confusion_matrix_250.csv")
# We only need TSS data for the anova and figures, so extract
# relevant columns.
cm_TSS <- confusion_matrix_data[,c(1,2,12)]
colnames(cm_TSS) <- c("Municipality", "Species", "TSS")
# convert species names to scientific names
cm_TSS <- cm_TSS %>%
  mutate(Species = recode(Species, "adder" = "Vipera berus", "gladde_slang" =         "Coronella austriaca", "hazelworm" = "Anguis fragilis",                      "levendbarende_hagedis" = "Zootoca vivipara", "muurhagedis" =                "Podarcis muralis", "ringslang" = "Natrix helvetica",
        "zandhagedis" = "Lacerta agilis"))
# extact sens and spec data:
cm_sens_spec <- confusion_matrix_data[,c(1,2,10,11)]
colnames(cm_sens_spec) <- c("Municipality", "Species", "Sensitivity", "Specificity")
cm_sens_spec <- cm_sens_spec %>%
  mutate(Species = recode(Species, "adder" = "Vipera berus", "gladde_slang" =         "Coronella austriaca", "hazelworm" = "Anguis fragilis",                      "levendbarende_hagedis" = "Zootoca vivipara", "muurhagedis" =                "Podarcis muralis", "ringslang" = "Natrix helvetica",
                          "zandhagedis" = "Lacerta agilis"))
####################################################################
# First, look at linear model and ANOVA to find out if the species or the municipality of the local applications are more important for determining the performance (=TSS) of the ensemble models at local scale.
# set as factor
cm_TSS$Species <- as.factor(cm_TSS$Species)
cm_TSS$Municipality <- as.factor(cm_TSS$Municipality)
TSS_model <- lm(TSS ~ Municipality + Species, data = cm_TSS)
# check contributions to regression
summary(TSS_model) # summary of the model, adjusted R squared of 0.6353
contrasts(cm_TSS$Municipality)
contrasts(cm_TSS$Species)
plot(TSS_model) # plots look decent, quite nice model fit.
# Now to perform a two-way ANOVA on this linear model to see how important species and municipality are for explaining the variance in TSS-values.
anova(TSS_model)
############################################################################
# Now to generate the figures. Figure 1: TSS per species (boxplot), plot national performance as 'dot' in the same figure. Figure 2 (or 1b): TSS per municipality (boxplot).
# First, add columns with overall performance and sample size based on species column.
cm_TSS$National_TSS <- ifelse(cm_TSS$Species == "Vipera berus", 0.942, 
                              ifelse(cm_TSS$Species == "Coronella austriaca", 0.919,
                              ifelse(cm_TSS$Species == "Anguis fragilis", 0.828,
                              ifelse(cm_TSS$Species == "Zootoca vivipara", 0.814,
                              ifelse(cm_TSS$Species == "Podarcis muralis", 0.825, 
                              ifelse(cm_TSS$Species == "Natrix helvetica", 0.757, 
                              ifelse(cm_TSS$Species == "Lacerta agilis", 0.898, NA)))))))
# sample size as vector
sample_size <- cm_TSS %>% group_by(Species) %>% summarize(n = n())
cm_TSS <- cm_TSS %>%
  left_join(sample_size) %>%
  mutate(n = paste0(Species,"\n", "n = ", n))
cm_TSS$n <- c('n = 5','n = 9','n = 7','n = 5','n = 5','n = 6','n = 9','n = 8','n = 7','n = 5','n = 5','n = 6','n = 9','n = 8','n = 7','n = 5','n = 9','n = 9','n = 8','n = 7','n = 5','n = 6','n = 9','n = 8','n = 7','n = 5','n = 6','n = 9','n = 8','n = 7','n = 5','n = 5','n = 6','n = 9','n = 1','n = 8','n = 7','n = 8','n = 6','n = 9','n = 8')
# dummy dataframe for national value:
dummy_data <- data.frame(Legend = "National TSS Value", x = NA, y = NA)
# and dataframe with national values:
Species <- c("Vipera berus","Coronella austriaca", "Anguis fragilis",                      "Zootoca vivipara", "Podarcis muralis","Natrix helvetica",
             "Lacerta agilis")
TSS <- c(0.942, 0.919, 0.828, 0.814, 0.825, 0.757, 0.898)
national_TSS <- data.frame(Species, TSS)
labs_species <- c('A. fragilis', 'C. austriaca', 'L. agilis', 'N. helvetica', 'P. muralis', 'V. berus', 'Z. vivipara')
sample_sizes_spec <- c('6', '5', '7', '8', '1', '5', '9')
x_labels <- paste0(labs_species, "\n(n = ", sample_sizes_spec, ")")
# figure:
cm_plot <- ggplot(data = cm_TSS, mapping = aes(x = Species, y = TSS)) +
  geom_boxplot(fill = 'lightgray', outlier.shape = NA) +
  theme_bw() +
  # Keep jittered points black but map fill for dodging
  #geom_jitter(aes(fill = Species),  # Add fill mapping for proper dodging
   #           color = 'black',  # Keep points black
    #          position = position_jitterdodge(0.1),
     #         alpha = 0.5) +
  labs(x = "Species", y = "TSS-value") +
  # Apply custom x-axis labels with formatting
  scale_x_discrete(labels = c(
    "<i>A. fragilis</i><br><b>(n = 6)</b>", 
    "<i>C. austriaca</i><br><b>(n = 5)</b>",
    "<i>L. agilis</i><br><b>(n = 7)</b>",
    "<i>N. helvetica</i><br><b>(n = 8)</b>",
    "<i>P. muralis</i><br><b>(n = 1)</b>",
    "<i>V. berus</i><br><b>(n = 5)</b>",
    "<i>Z. vivipara</i><br><b>(n = 9)</b>"
  )) +
  theme(axis.title = element_text(colour = 'grey2', size = 32, face = 'bold'), 
        # Format x-axis text using ggtext::element_markdown
        axis.text.x = element_markdown(size = 24, colour = 'grey2', face = 'bold'),
        axis.text.y = element_text(colour = 'grey2', size = 24, face = 'bold'),
        axis.ticks = element_blank()) +
  theme(legend.background = element_rect(fill = 'white', colour = 'black',
                                         linewidth = 1, linetype = 'solid'),
        legend.title = element_text(size = 28, face = "bold")) +
  theme(legend.text = element_text(size = 25, face = "italic"),
        legend.key.size = unit(5, 'line'),
        legend.key = element_rect(colour = "black")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_text(margin = margin(t = 20))) +
  # Add blue dots for national TSS values (blue in plot and legend)
  geom_point(data = national_TSS, aes(x = Species, y = TSS, color = "National TSS Value"),
             size = 7, show.legend = TRUE) +
  # Use scale_color_manual to set the blue dot for "National TSS Value"
  scale_color_manual(name = NULL, values = c("National TSS Value" = "blue")) +
  # Ensure species colors remain black in the legend
  guides(
    fill = guide_legend(order = 1), 
    color = guide_legend(order = 2, override.aes = list(size = 7)),
    shape = guide_legend(override.aes = list(color = "black"))
  )

# save:
ggsave("./output/figuren/verslag_250/cm_species2.png", cm_plot, dpi = 300,
       height = 18, width = 22)
############################################################################
# second figure, per municipality.
# columns of national data in vector
Species <- c("Vipera berus","Coronella austriaca", "Anguis fragilis",                      "Zootoca vivipara", "Podarcis muralis","Natrix helvetica",
             "Lacerta agilis")
TSS <- c(0.942, 0.919, 0.828, 0.814, 0.825, 0.757, 0.898)
Municipality <- c('NL','NL' ,'NL' ,'NL' ,'NL' ,'NL' ,'NL')
# turn into dataframe
national_data <- data.frame(Species, TSS, Municipality)
# take relevant columns from initial TSS data
mun_data <- cm_TSS[,1:3]
# merge both dataframes
mun_TSS <- rbind(mun_data, national_data)
# add sample size:
mun_TSS$n <- c('3','3','3','6','6','6','6','6','6','6','6', '6', '6','6','6','2','2','3', '3', '3','5','5', '5', '5','5','5','5', '5', '5','5','7','7', '7', '7','7','7', '7','1', '3', '3', '3', '7','7', '7', '7','7','7', '7')
############################################################################
# Now put municipality into figure:
# Define the sample sizes
sample_sizes <- c('3', '6', '6', '2', '3','5', '5', '7', '1', '3', '7')

# Get the unique municipality names in the correct order
municipality_levels <- c(setdiff(unique(mun_TSS$Municipality), "NL"), "NL")

# Combine municipality names with sample sizes (horizontally aligned below)
municipality_labels <- paste0(municipality_levels, "\n(n = ", sample_sizes, ")")

# Update Municipality factor to ensure proper ordering
mun_TSS$Municipality <- factor(mun_TSS$Municipality, levels = municipality_levels)

#plot:
cm_plot_gemeente <- ggplot(data = mun_TSS, mapping = aes(x = Municipality, y = TSS)) +
  # Boxplot with conditional fill
  geom_boxplot(aes(fill = Municipality == "NL"), outlier.shape = NA) +
  scale_fill_manual(values = c("FALSE" = "lightgray", "TRUE" = "blue")) +  # Set colors
  theme_bw() +
  # Jitter points, mapping fill for proper dodging
  #geom_jitter(aes(fill = Municipality == "NL"),
   #           color = 'black',  # Keep points black
    #          position = position_jitterdodge(0.1),
     #         alpha = 0.5) +
  labs(x = "Municipality", y = "TSS") +
  scale_x_discrete(labels = municipality_labels) +  # Apply custom labels
  theme(axis.title = element_text(colour = 'grey2', size = 32, face = "bold"),
        axis.text = element_text(colour = 'grey2', size = 24, face = "bold"),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5, face = "bold", size = 21),
        axis.ticks = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.background = element_rect(fill = 'white', colour = 'black',
                                         linewidth = 1, linetype = 'solid'),
        legend.title = element_text(size = 28, face = "bold")) +
  theme(legend.text = element_text(size = 25, face = "bold"),
        legend.key.size = unit(5, 'line'),
        legend.key = element_rect(colour = "black")) +
  theme(axis.title.x = element_text(margin = margin(t = 20))) +
  theme(legend.position = 'none')
# save:
ggsave("./output/figuren/verslag_250/cm_municipality.png", cm_plot_gemeente, dpi = 300, height = 18, width = 22)

# now combine species and municipalities plots:
combi_plot <- cm_plot / cm_plot_gemeente+ 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 40, face = "bold"))
combi_plot
# save images:
ggsave("./output/figuren/verslag_250/cm_spec_munip.png", combi_plot, height = 25, width = 25, dpi = 320)
#############################################################################
# finally, plot sensitivity and specificity values per species and per municipality
# data:
cm_sens_spec
# add sample size
sample_size <- cm_sens_spec %>% group_by(Species) %>% summarize(n = n())
cm_sens_spec <- cm_sens_spec %>%
  left_join(sample_size) %>%
  mutate(n = paste0(Species,"\n", "n = ", n))
cm_TSS$n <- c('n = 5','n = 9','n = 7','n = 5','n = 5','n = 6','n = 9','n = 8','n = 7','n = 5','n = 5','n = 6','n = 9','n = 8','n = 7','n = 5','n = 9','n = 9','n = 8','n = 7','n = 5','n = 6','n = 9','n = 8','n = 7','n = 5','n = 6','n = 9','n = 8','n = 7','n = 5','n = 5','n = 6','n = 9','n = 1','n = 8','n = 7','n = 8','n = 6','n = 9','n = 8')
# make figure:
labs_species <- c('A. fragilis', 'C. austriaca', 'L. agilis', 'N. helvetica', 'P. muralis', 'V. berus', 'Z. vivipara')
sample_sizes_spec <- c('6', '5', '7', '8', '1', '5', '9')
x_labels <- paste0(labs_species, "\n(n = ", sample_sizes_spec, ")")
# need to put scores under each other to be able to plot:
long_data <- melt(cm_sens_spec, id.vars = "Species", 
                  measure.vars = c("Sensitivity", "Specificity"), 
                  variable.name = "Metric", value.name = "Score")
# figure:
sens_spec_species <- ggplot(data = long_data, aes(x = factor(Species, levels = unique(Species)), y = Score, fill = Metric)) + 
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +  # Dodge places boxes side by side
  theme_bw() +
  labs(x = "Species", y = "Score Value") +
  scale_fill_manual(values = c("Sensitivity" = "skyblue", "Specificity" = "orange")) + 
  # Apply custom x-axis labels with formatting
  scale_x_discrete(labels = c(
   "<i>C. austriaca</i><br><b>(n = 5)</b>",
   "<i>Z. vivipara</i><br><b>(n = 9)</b>",
    "<i>L. agilis</i><br><b>(n = 7)</b>",
    "<i>V. berus</i><br><b>(n = 5)</b>",
    "<i>A. fragilis</i><br><b>(n = 6)</b>",
    "<i>N. helvetica</i><br><b>(n = 8)</b>",
    "<i>P. muralis</i><br><b>(n = 1)</b>"
  )) +
  theme(axis.title = element_text(colour = 'grey2', size = 32, face = 'bold'), 
        # Format x-axis text using ggtext::element_markdown
        axis.text.x = element_markdown(size = 24, colour = 'grey2', face = 'bold'),
        axis.text.y = element_text(colour = 'grey2', size = 24, face = 'bold'),
        axis.ticks = element_blank()) +
  theme(axis.title.x = element_text(margin = margin(t = 20))) +
  theme(axis.title.y = element_text(margin = margin(t = 30))) +
  guides(fill = guide_legend(keysize = 6)) +
  theme(legend.text = element_text(size = 30, colour = 'grey2'),
        legend.title = element_text(size = 34, colour = 'grey2', face = 'bold'),legend.background = element_rect(fill = "white", colour = "black", linewidth = 1))
ggsave("./output/figuren/verslag_250/sens_spec_species.png", sens_spec_species, height = 16, width = 25, dpi = 320)
###########################################################################
# Now, also per municipality
# data:
cm_sens_spec
# add sample size
cm_sens_spec$n <- c('3','3','3','6','6','6','6','6','6','6','6', '6', '6','6','6','2','2','3', '3', '3','5','5', '5', '5','5','5','5', '5', '5','5','7','7', '7', '7','7','7', '7','1', '3', '3', '3')
sample_sizes <- c('3', '6', '6', '2', '3','5', '5', '7', '1', '3')


# Get the unique municipality names in the correct order
municipality_levels <- c(unique(cm_sens_spec$Municipality))

# Combine municipality names with sample sizes (horizontally aligned below)
municipality_labels <- paste0(municipality_levels, "\n(n = ", sample_sizes, ")")

labels <- c("Deurne \n (n = 3)", "Ede \n (n = 6)", "Epe \n (n = 6)", "Horst \n (n = 2)", "Huizen \n (n = 3)", "Ommen \n (n = 5)", "Roermond \n (n = 5)", "Westerveld \n (n = 7)", "Wormerland \n (n = 1)", "Zutphen \n (n = 3)")

# need to put scores under each other to be able to plot:
long_data <- melt(cm_sens_spec, id.vars = "Municipality", 
                  measure.vars = c("Sensitivity", "Specificity"), 
                  variable.name = "Metric", value.name = "Score")
long_data$Municipality <- as.factor(long_data$Municipality)
# figure:
sens_spec_munip <- ggplot(data = long_data, aes(x = factor(Municipality, levels = unique(Municipality)), y = Score, fill = Metric)) + 
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +  # Dodge places boxes side by side
  theme_bw() +
  labs(x = "Municipality", y = "Score Value") +
  scale_fill_manual(values = c("Sensitivity" = "skyblue", "Specificity" = "orange")) + 
  # Apply custom x-axis labels with formatting
  #scale_x_discrete(labels = labels) +
  scale_x_discrete(labels = c(
    "Deurne<br><b>(n = 3)</b>",
    "Ede<br><b>(n = 6)</b>",
    "Epe<br><b>(n = 6)</b>",
    "Horst<br><b>(n = 2)</b>",
    "Huizen<br><b>(n = 3)</b>",
    "Ommen<br><b>(n = 5)</b>",
    "Roermond<br><b>(n = 5)</b>",
    "Westerveld<br><b>(n = 7)</b>",
    "Wormerland<br><b>(n = 1)</b>",
    "Zutphen<br><b>(n = 3)</b>"
  )) +
  theme(axis.title = element_text(colour = 'grey2', size = 32, face = 'bold'), 
        # Format x-axis text using ggtext::element_markdown
        axis.text.x = element_markdown(size = 24, colour = 'grey2', face = 'bold'),
        axis.text.y = element_text(colour = 'grey2', size = 24, face = 'bold'),
        axis.ticks = element_blank()) +
  theme(axis.title.x = element_text(margin = margin(t = 25))) +
  theme(axis.title.y = element_text(margin = margin(t = 40))) +
  guides(fill = guide_legend(keysize = 6)) +
  theme(legend.text = element_text(size = 26, colour = 'grey2'),
        legend.title = element_text(size = 30, colour = 'grey2', face = 'bold'),legend.background = element_rect(fill = "white", colour = "black", linewidth = 1))
# save image
ggsave("./output/figuren/verslag_250/sens_spec_munip.png",sens_spec_munip, height = 16, width = 25, dpi = 320)
