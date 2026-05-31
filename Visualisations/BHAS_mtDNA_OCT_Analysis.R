# =============================================================================
# BHAS mtDNA Haplogroup × Optic Nerve Head Morphology Analysis
# -----------------------------------------------------------------------------
# Companion code for the manuscript:
#   "Mitochondrial Haplogroup J Is Associated with Reduced Optic Nerve Head
#    Area in a Healthy Ageing Cohort"
#
# Pipeline:
#   1. Setup and shared helper functions (parent map, haplogroup classifier)
#   2. Sunburst visualisation of haplogroup hierarchy
#   3. Linear mixed-effects models for haplogroup × OCT phenotype associations
#      with Benjamini-Hochberg FDR correction
#
#
# Expected input dataframe `df` columns:
#   PatientID       - participant identifier (used for random intercept)
#   Haplogroup      - Haplogrep 3 terminal sub-haplogroup assignment
#   Age, Gender     - demographic covariates
#   OCT phenotypes  - BMO_Area, MRW_Global, RNFLMean_G, IOP, etc.
# =============================================================================

# 1. SETUP --------------------------------------------------------------------

library(tidyverse)
library(stringr)
library(plotly)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(writexl)


# 2. SHARED HELPERS -----------------------------------------------------------

# 2.1 Phylogenetic parent map (simplified PhyloTree-style hierarchy)
parent_map <- c(
  "L3" = "L1-6", "M" = "L3", "N" = "L3", "L0" = "L-Root",
  "L1" = "L1-6", "L2" = "L1-6", "L4" = "L1-6", "L5" = "L1-6", "L6" = "L1-6",
  "C" = "M", "Z" = "M", "D" = "M", "E" = "M", "G" = "M", "Q" = "M",
  "R" = "N", "I" = "N", "W" = "N", "X" = "N", "Y" = "N",
  "A" = "N", "S" = "N", "O" = "N",
  "R0" = "R", "HV" = "R0", "H" = "HV", "V" = "HV",
  "JT" = "R", "J" = "JT", "T" = "JT",
  "U" = "R", "K" = "U8", "U8" = "U",
  "B" = "R", "F" = "R", "P" = "R"
)

# 2.2 Strip Haplogrep modifiers (e.g. "H1a+5G", "H*" -> "H1a", "H")
clean_hg <- function(hgs) {
  str_split(hgs, "[\\+\\*\\@\\']", simplify = TRUE)[, 1]
}

# 2.3 Walk a terminal haplogroup up the tree, returning the full lineage
get_lineage_vector <- function(hg) {
  lineage <- character()
  curr <- hg
  visited <- c()
  while (nchar(curr) > 0) {
    if (curr %in% visited) break
    visited <- c(visited, curr)
    lineage <- c(lineage, curr)
    if (curr %in% names(parent_map)) {
      curr <- parent_map[[curr]]
    } else if (grepl("-", curr)) {
      break
    } else {
      if (nchar(curr) > 1) curr <- str_sub(curr, 1, nchar(curr) - 1) else curr <- ""
    }
  }
  return(unique(lineage))
}

# 2.4 Dynamic 1-vs-all classifier: returns `target_group` if the input haplogroup is a descendant, otherwise "Other".
get_dynamic_comparison <- function(hg, target_group) {
  if (is.na(hg) || hg == "") return("Other")
  curr <- str_split(hg, "[\\+\\*\\@\\']", simplify = TRUE)[, 1]
  if (is.na(curr)) return("Other")
  visited <- c()
  while (nchar(curr) > 0) {
    if (curr == target_group) return(curr)
    if (curr %in% visited) break
    visited <- c(visited, curr)
    if (curr %in% names(parent_map)) {
      curr <- parent_map[[curr]]
    } else if (grepl("-", curr)) {
      return("Other")
    } else {
      if (nchar(curr) > 1) curr <- str_sub(curr, 1, nchar(curr) - 1) else curr <- ""
    }
  }
  return("Other")
}


# 3. SUNBURST VISUALISATION ---------------------------------------------------

# 3.1 Per-participant lineage expansion and node counts
raw_haplos    <- df$Haplogroup
total_n       <- length(raw_haplos)
cleaned_input <- clean_hg(raw_haplos)
unique_hgs    <- unique(cleaned_input)
lineage_list  <- setNames(lapply(unique_hgs, get_lineage_vector), unique_hgs)

expanded_df <- data.frame(Raw_Input = cleaned_input) %>%
  left_join(
    data.frame(Raw_Input = names(lineage_list), Node = I(lineage_list)) %>%
      unnest(Node),
    by = "Raw_Input"
  )

final_counts <- expanded_df %>%
  group_by(Node) %>%
  summarise(Total_Participants = n()) %>%
  ungroup()

# 3.2 Assign each node to a major-clade colour family
major_clades <- c("L0", "L1-6", "M", "N", "R", "H", "V", "J", "T", "U", "K", "I", "W", "X", "A", "B", "C", "D", "G", "F", "Z")

get_color_group <- function(node) {
  curr <- node
  visited <- c()
  while (nchar(curr) > 0) {
    if (curr %in% visited) break
    visited <- c(visited, curr)
    if (curr %in% major_clades) return(curr)
    if (curr %in% names(parent_map)) {
      curr <- parent_map[[curr]]
    } else if (grepl("-", curr)) {
      return("Root")
    } else {
      if (nchar(curr) > 1) curr <- str_sub(curr, 1, nchar(curr) - 1) else curr <- ""
    }
  }
  return("Other")
}

plotly_df <- final_counts %>%
  rowwise() %>%
  mutate(
    Color_Group   = get_color_group(Node),
    Percentage    = (Total_Participants / total_n) * 100,
    Display_Label = paste0(Node, " (", round(Percentage, 1), "%)")
  ) %>%
  ungroup()

# 3.3 Build colour palette and parent assignments
unique_families <- unique(plotly_df$Color_Group)
my_palette      <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique_families))
names(my_palette) <- unique_families
plotly_df$Hex_Color <- my_palette[plotly_df$Color_Group]

get_immediate_parent <- function(node) {
  if (node %in% names(parent_map)) return(parent_map[[node]])
  if (grepl("-", node)) return("")
  if (nchar(node) > 1) return(str_sub(node, 1, nchar(node) - 1))
  return("")
}

plotly_df <- plotly_df %>%
  rowwise() %>%
  mutate(Parent = get_immediate_parent(Node)) %>%
  ungroup() %>%
  mutate(Parent = ifelse(Parent %in% Node, Parent, ""))

# 3.4 Render sunburst
plot_ly(
  data         = plotly_df,
  ids          = ~Node,
  labels       = ~Display_Label,
  parents      = ~Parent,
  values       = ~Total_Participants,
  type         = "sunburst",
  branchvalues = "total",
  marker       = list(colors = ~Hex_Color,
                      line   = list(color = "white", width = 0.5)),
  hovertemplate = "<b>%{label}</b><br>Count: %{value}<br>Percent: %{percentRoot:.1%}<extra></extra>"
) %>%
  layout(title = paste0("mtDNA Sunburst (N = ", total_n, ")"))


# 4. LINEAR MIXED-EFFECTS MODELS ----------------------------------------------

# 4.1 Target haplogroups (those meeting the ≥400-participant threshold or with
#     established LHON / glaucoma literature support) and OCT phenotypes.
target_haplos <- c("J", "J1c", "H", "U", "K", "T", "U5a", "U5b")

primary_phenotypes <- c("BMO_Area", "MRW_Global", "RNFLMean_G", "IOP")

sector_phenotypes  <- c("MRW_Tmp", "MRW_TS", "MRW_TI", "MRW_Nas", "MRW_NS", "MRW_NI",
                        "RNFLMean_T", "RNFLMean_TS", "RNFLMean_TI",
                        "RNFLMean_N", "RNFLMean_NS", "RNFLMean_NI")

# 4.2 Fit one lmer for every (haplogroup × phenotype) combination
run_lmm <- function(data, target_haplo, phenotype) {
  temp_data <- data %>%
    mutate(
      Comparison = sapply(Haplogroup, function(x) get_dynamic_comparison(x, target_haplo)),
      Comparison = relevel(as.factor(Comparison), ref = "Other")
    )
  
  form  <- as.formula(paste0("`", phenotype, "` ~ Comparison + Age + Gender + (1|PatientID)"))
  model <- lmer(form, data = temp_data)
  
  tidy(model, conf.int = TRUE) %>%
    filter(effect == "fixed", grepl("Comparison", term)) %>%
    mutate(Target_Haplogroup = target_haplo,
           Phenotype         = phenotype) %>%
    select(Target_Haplogroup, Phenotype, term, estimate, std.error,
           conf.low, conf.high, p.value)
}

# 4.3 Run the full grid (all haplogroups × all phenotypes)
run_grid <- function(data, haplos, phenos) {
  expand_grid(target = haplos, phenotype = phenos) %>%
    pmap_dfr(~ run_lmm(data, ..1, ..2))
}

primary_results <- run_grid(df, target_haplos, primary_phenotypes)
sector_results  <- run_grid(df, target_haplos, sector_phenotypes)

# 4.4 Benjamini-Hochberg FDR correction (applied separately within each
#     analysis family, as reported in the manuscript)
primary_results <- primary_results %>%
  mutate(p_adj_BH = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value)

sector_results <- sector_results %>%
  mutate(p_adj_BH = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value)

# 4.5 Export
write_xlsx(
  list(Primary = primary_results, Sector = sector_results),
  path = "LMM_Results.xlsx"
)
