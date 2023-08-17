library(dplyr)
library(xlsx)
library(ggpubr)
library(tidyverse)
#install compassR
#devtools::install_github("YosefLab/compassR")
library(compassR)
library(ggrepel)
library(forcats)
source("plot_pdf_png.r")

project_dir <- "project_path/"
path_output <- paste0(project_dir,"results/")
path_deposit <- paste0(project_dir,"dataDeposit/")


compass_settings <- CompassSettings$new(
    user_data_directory = paste0(project_dir,"input"),
    metabolic_model_directory = system.file("extdata", "RECON2", package = "compassR", mustWork = TRUE),
    cell_id_col_name = "cell_id",
    gene_id_col_name = "HGNC.symbol",
    
)
compass_data <- CompassData$new(compass_settings)
compass_analyzer <- CompassAnalyzer$new(compass_settings)

fib1_cell_ids <-
    compass_data$cell_metadata %>%
    filter(final_celltype_refined == "Fib_1") %>%
    pull(cell_id)
fib2_cell_ids <-
    compass_data$cell_metadata %>%
    filter(final_celltype_refined == "Fib_2") %>%
    pull(cell_id)

reaction_consistencies <- read.csv(paste0(project_dir, "results/reaction_consistencies.csv"), row.names = 1, header = F)
colnames(reaction_consistencies) <- reaction_consistencies[1,]
reaction_consistencies <- reaction_consistencies[-1,]
reaction_consistencies.matrix <- 
  sapply(reaction_consistencies, function(x) as.numeric(as.character(x)))
rownames(reaction_consistencies.matrix) <- rownames(reaction_consistencies)


wilcoxon_results <- compass_analyzer$conduct_wilcoxon_test(
    reaction_consistencies.matrix,
    fib2_cell_ids,
    fib1_cell_ids,
    for_metareactions = FALSE
)


cohens_d_by_subsystem2 <-
    wilcoxon_results %>%
    left_join(
        select(compass_data$reaction_partitions, "reaction_id", "reaction_no_direction"),
        by = "reaction_id"
    ) %>%
    left_join(
        compass_data$reaction_metadata,
        by = "reaction_no_direction"
    ) %>%
    # Keep only "confident reactions", as defined in our paper.
    #filter(!is.na(EC_number)) %>%
    #filter(confidence == "0" | confidence == "4") %>%
    # Keep only "interesting subsystems", as defined in our paper.
   #filter(!(subsystem == "Miscellaneous" | subsystem == "Unassigned")) %>%
  #  filter(!(startsWith(subsystem, "Transport") | startsWith(subsystem, "Exchange"))) %>%
    # Keep only subsystems of non-negligible size.
    group_by(subsystem) %>%
    filter(n() > 5) %>%
    ungroup() %>%
    # Order subsystems in a manner that will lend itself to a visually aesthetic plot.
    mutate(
        subsystem_priority = factor(subsystem) %>%
        fct_reorder2(
            cohens_d,
            adjusted_p_value,
            .fun = function(cohens_d, adjusted_p_value) {
                abs(median(cohens_d[adjusted_p_value < 0.1]))
            },
            .desc = FALSE
        )
    )

## forest plot
p2 <- ggplot(
    cohens_d_by_subsystem2,
    aes(
        x = subsystem_priority,
        y = cohens_d,
        color = if_else(cohens_d > 0, "Fib1", "Fib2"),
        alpha = if_else(adjusted_p_value < 0.01, "significant", "insignificant")
    )
) +
ggtitle("Fib1 vs Fib2 activated metabolic pathways ") +
xlab("") + ylab("Cohen's d") +
scale_color_manual(
    values = c(Fib1 = "#ca0020", Fib2 = "#0571b0"),
    guide = FALSE
) +
scale_alpha_manual(
    name = "",
    values = c(significant = 1, insignificant = 0.25),
    labels = c(significant = "BH-adjusted p-value < 0.01", insignificant = "insignificant")
) +
coord_flip() +
geom_point() +
geom_hline(yintercept = 0, linetype = "dashed") +
theme_bw() +
theme(legend.position = "bottom", legend.direction = "horizontal")
pdf_and_png(p1, paste0(path_output,"forest_plot_with_filter_fib1_vs_fib2"), width = 8, h_to_w_ratio = 1.2)


facets <- c( "Glycolysis", "TCA cycle", "Fatty acid oxidation", "Amino acid metabolism")

compass_scores_by_cell_type <-
    wilcoxon_results %>%
    left_join(
        select(compass_data$reaction_partitions, "reaction_id", "reaction_no_direction"),
        by = "reaction_id"
    ) %>%
    left_join(
        compass_data$reaction_metadata,
        by = "reaction_no_direction"
    ) %>%
    # Keep only "confident reactions", as defined in our paper.
    filter(!is.na(EC_number)) %>%
    filter(confidence == "0" | confidence == "4") %>%
    # Exclude non-mitochondrially localized reactions from TCA.
    mutate(subsystem = case_when(
        reaction_id == "SPMDOX_pos" ~ "Arginine and Proline Metabolism",
        subsystem == "Citric acid cycle" & !grepl("[m]", formula, fixed = TRUE) ~ "Other",
        TRUE ~ subsystem
    )) %>%
    # Assign reactions to the appropriate subsystem.
    mutate(
        subsystem_priority = factor(subsystem) %>%
        fct_recode(
            "Glycolysis" = "Glycolysis/gluconeogenesis",
            "TCA cycle" = "Citric acid cycle"
        ) %>%
        fct_collapse("Amino acid metabolism" = c(
            "Alanine and aspartate metabolism",
            "Arginine and Proline Metabolism",
            "beta-Alanine metabolism",
            "Cysteine Metabolism",
            "D-alanine metabolism",
            "Folate metabolism",
            "Glutamate metabolism",
            "Glycine, serine, alanine and threonine metabolism",
            "Histidine metabolism",
            "Lysine metabolism",
            "Methionine and cysteine metabolism",
            "Taurine and hypotaurine metabolism",
            "Tryptophan metabolism",
            "Tyrosine metabolism",
            "Urea cycle",
            "Valine, leucine, and isoleucine metabolism"
        )) %>%
           fct_collapse("Fatty acid oxidation" = c(
            "Fatty acid oxidation", "Fatty acid synthesis")) %>% 
        fct_other(keep = facets) %>%
        fct_relevel(facets)
    ) %>%
    # Keep only the subsystems for which we want to plot a facet.
    filter(subsystem_priority != "Other") %>%
    # Lower-bound the adjusted p-value.
    mutate(adjusted_p_value = if_else(
        subsystem_priority == "Amino acid metabolism" & adjusted_p_value <= 1e-12,
        1e-12,
        adjusted_p_value
    )) %>%
    # Assign descriptive labels to various reactions.
    mutate(label = case_when(
        reaction_id == "PGM_neg" ~ "phosphoglycerate mutase (PGAM)",
        reaction_id == "LDH_L_neg" ~ "lactate dehydrogenase",
        reaction_id == "PDHm_pos" ~ "pyruvate dehydrogenase (PDH)",
        reaction_id == "TPI_neg" ~ "triosephosphate isomerase (DHAP forming)",
        reaction_id == "FACOAL1821_neg" ~ "long-chain fatty-acid-CoA ligase",
        reaction_id == "r1257_pos" ~ "long-chain fatty-acid-CoA ligase",
        reaction_id == "FACOAL1831_neg" ~ "long-chain fatty-acid-CoA ligase",
        reaction_id == "CSNATr_neg" ~ "carnitine O-acetyltransferase",
        reaction_id == "C160CPT1_pos" ~ "carnitine O-palmitoyltransferase",
        reaction_id == "ACONTm_pos" ~ "aconitate hydratase",
        reaction_id == "SUCOASm_pos" ~ "succinate-CoA ligase",
        reaction_id == "AKGDm_pos" ~ "alpha-ketoglutarate dehydrogenase",
        reaction_id == "SUCD1m_pos" ~ "succinate dehydrogenase",
        reaction_id == "ICDHyrm_pos" ~ "isocitrate dehydrogenase",
        reaction_id == "CK_pos" ~ "creatine\nkinase",
        reaction_id == "PGCD_pos" ~ "phosphoglycerate dehydrogenase",
        reaction_id == "ARGSS_pos" ~ "arginosuccinate synthase",
        reaction_id == "r0281_neg" ~ "putrescine diamine oxidase",
        reaction_id == "SPMDOX_pos" ~ "spermidine dehydrogenase (spermidine -> GABA)",
        reaction_id == "ARGDCm_pos" ~ "arginine decarboxylase",
        reaction_id == "AGMTm_pos" ~ "agmatinase",
        reaction_id == "GHMT2r_pos" ~ "serine hydroxymethyltransferase",
        reaction_id == "AHC_pos" ~ "adenosylhomocysteinase",
        reaction_id == "METAT_pos" ~ "methionine adenosyltransferase",
        reaction_id == "METS_pos" ~ "methionine\nsynthase",
        reaction_id == "ARGN_pos" ~ "arginase",
        TRUE ~ ""
    ))


#write to excel with multiple sheets
library(xlsx)

for(i in facets){
  output <- compass_scores_by_cell_type %>% filter(subsystem_priority == i) %>% arrange(desc(cohens_d))
  
  write.xlsx(as.data.frame(output), file= paste0(path_output,"compass_score_4_category.xlsx"), sheetName=i, row.names=FALSE, append = T)
}

q <- ggplot(
    compass_scores_by_cell_type,
    aes(
        x = cohens_d,
        y = -log10(adjusted_p_value),
        color = subsystem_priority
    )
) +
ggtitle("Differential COMPASS Scores for Fib2 vs. Fib1") +
xlab("Cohen's d") + ylab("-log(BH-adjusted p-value)") +
xlim(-2.2, 2.2) +
facet_wrap(vars(subsystem_priority), scales = "free_y", ncol = 2) +
scale_color_manual(values = c(
    "Glycolysis" = "#662D8C",
    "TCA cycle" = "#B87013",
    "Fatty acid oxidation" = "#0B0D9D",
    "Amino acid metabolism" = "#B82130"
)) +
guides(color = FALSE) +
geom_point(size = 1, alpha = 0.5) +
geom_hline(yintercept = 1, linetype="dashed", color = "blue") +
geom_vline(xintercept = 0, linetype="dashed", color = "blue") +
geom_text_repel(
    aes(label = label),
    min.segment.length = 0.1,
    point.padding = 0.5,
    size = 2,
    seed = 7
) +
theme_bw()

pdf_and_png(q, paste0(path_output,"volcano_plot_with_random_labels"), width = 8, h_to_w_ratio = 0.8)