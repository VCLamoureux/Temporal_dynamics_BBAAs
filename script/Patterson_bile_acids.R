# script use to generate counts of bile acids per organ regardless of ZT time
## set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth")

# Load required libraries
library(tidyverse)
library(duckplyr)
library(purrr)
library(ComplexUpset)
library(colorspace)
library(mgcv)

# import the table
table_simplified <- read_csv("Erin_patterson_lab/All_ikoo_area_1_2024_11_11_16_54_24 (1).csv")
BAs_to_remove <- read_csv("BA_to_remove_Patterson.csv")
cols_to_remove <- BAs_to_remove$BA

table_simplified_1 <- table_simplified |> 
  dplyr::select(-any_of(cols_to_remove)) |> 
  dplyr::rename(TissueType = `Tissue Type`) 

table_simplified_ext <- table_simplified_1 |> 
  dplyr::select(-c(`Injection order`, `Sample Name...4`, `Sample Name...117`)) |>
  dplyr::mutate(ZT = as.numeric(ZT)) #|> 
  #bind_rows(table_simplified |> 
      #dplyr::filter(ZT == 0) |> 
      #dplyr::mutate(ZT = 24))

# remove d4, d5, and some samples column and the ISF 
table_simplified_red <- table_simplified_ext |> 
  dplyr::select(-contains("-d4"), -contains("Sample"), -contains("-d5"), -"choloic acid [M-3H2O+H]+_9859")

# apply some filtering for a threshold; if the value is not at least 10 000 replace with 0. 
table_simplified_red_thres <- table_simplified_red |> 
  dplyr::mutate(across(.cols = -c(TissueType, ZT), .fns  = ~ ifelse(. < 1e4, 0, .)))

# target only the bile acid column for applying the threshold for number of replicates
feature_cols <- table_simplified_red_thres |>  
  dplyr::select(-TissueType, -ZT) |>  
  colnames()

# Now we will apply the threshold for number of replicates (3 out of 5 n)
df2_thresholded <- table_simplified_red_thres |> 
  group_by(ZT, TissueType) |> 
  dplyr::mutate(across(all_of(feature_cols), ~ {
    n_nonzero <- sum(. > 0, na.rm = TRUE)
    if (n_nonzero <= 2) {
      0L
    } else {
      as.integer(. > 0)
    }
  })) |> ungroup()


# Pivot longer so you have one row per (Tissue Type, Bile acid) and independent of ZT time
df2_long <- df2_thresholded |>
  dplyr::select(-ZT) |> 
  pivot_longer(cols = -TissueType, names_to  = "Bile_acid", values_to = "Present")

# report presence/absence of each bile acid per tissue type
df2_any <- df2_long |> 
  group_by(TissueType, Bile_acid) |> 
  summarise(Present = as.integer(any(Present == 1)), .groups = "drop")

# Create a summary of number of bile acids per Tissue type 
tissue_counts <- df2_any |> 
  group_by(TissueType) |> 
  summarise(Number_bile_acids = sum(Present), .groups = "drop")

# Filter for the presence, groupby tissue type and list the name of the bile acids
organ_bileacids_df <- df2_any |>  
  dplyr::filter(Present == 1) |>  
  group_by(TissueType) |>  
  summarise(Bile_acids = list(Bile_acid), Count = length(Bile_acid), .groups = "drop")

# list all of them, 1 per row
organ_bileacids_df_1 <-  organ_bileacids_df |> 
  tidyr::unnest(cols = Bile_acids)

# define groups of bile acids

free_bas <- c(
  "(4R)-4-((3R,5R,6S,7R,9S,10R,12S,13R,17R)-3,6,7,12-tetrahydroxy-10,13-dimethylhexadecahydro-1H-cyclopenta[a]phenanthren-17-yl)pentanoic acid_12021",
  "(4R)-4-((3R,5R,6S,7R,9S,10R,12S,13R,17R)-3,6,7,12-tetrahydroxy-10,13-dimethylhexadecahydro-1H-cyclopenta[a]phenanthren-17-yl)pentanoic acid_12027",
  "(4R)-4-((3R,5S,7S,9S,10S,12S,13R,14S,17R)-3,7,12-trihydroxy-10,13-dimethylhexadecahydro-1H-cyclopenta[a]phenanthren-17-yl)pentanoic acid_12154",
  "(R)-4-((3R,5S,7R,8R,9S,10S,12S,13R,14S,17R)-3,7,12-trihydroxy-10,13-dimethylhexadecahydro-1H-cyclopenta[a]phenanthren-17-yl)pentanoic acid_10974",
  "choloic acid [M-3H2O+H]+_9859",
  "deoxycholic acid_33153",
  "hyocholic acid_10971",
  "hyocholic acid_10972",
  "hyocholic acid_10973",
  "hyocholic acid_10975",
  "hyocholic acid_10976",
  "hyocholic acid_34140",
  "methyl (4R)-4-((3R,5S,7R,9S,10S,13R,15R,17R)-3,7,15-trihydroxy-10,13-dimethylhexadecahydro-1H-cyclopenta[a]phenanthren-17-yl)pentanoate_15545",
  "methyl (R)-4-((3R,5S,7R,8S,9S,10S,13R,14R,17R)-3,7,14-trihydroxy-10,13-dimethylhexadecahydro-1H-cyclopenta[a]phenanthren-17-yl)pentanoate_11881"
)

aa_conj_cols <- c(
  "Ala-bMCA_17791",
  "Ala-CA_17790",
  "Arg-bMCA_23369",
  "Arg-CA_23368",
  "Arg-DCA_22370",
  "Asn-CA_20646",
  "Gln-CA_21582",
  "His-CA_22165",
  "Ile-Leu-CA_20580",
  "Ile-Leu-gMCA_20581",
  "Lys-bMCA_21589",
  "Lys-CA_21586",
  "Lys-CA_21588",
  "Lys-gMCA_21590",
  "Met-CA_21765",
  "Met-gMCA_21769",
  "Ornithine-CA_20652",
  "Phe-bMCA_22813",
  "Phe-CA_22812",
  "Ser-CA_18843",
  "Thr-CA_19755",
  "Trp-bMCA_25122",
  "Trp-CA_25121",
  "Tyr-bMCA_23804",
  "Tyr-CA_23803",
  "Val-CA_19625"
)


taurine <- c(
  "taurocholic acid_16621",
  "taurocholic acid_16624",
  "taurocholic acid_17778",
  "taurocholic acid_18976",
  "taurocholic acid_18978",
  "taurocholic acid_18979",
  "taurocholic acid_18981",
  "taurocholic acid_18984",
  "taurocholic acid_20155",
  "taurocholic acid_20157",
  "taurodeoxycholic acid_19107",
  "taurodeoxycholic acid_19109",
  "taurohyocholic acid_18975",
  "taurohyocholic acid_18980",
  "taurohyocholic acid_20152",
  "taurohyocholic acid_20153",
  "taurohyocholic acid_20154",
  "taurohyocholic acid_20160",
  "taurohyodeoxycholic acid_16761",
  "taurohyodeoxycholic acid_16765",
  "taurolithocholic acid_18046",
  "tauroursodeoxycholic acid_17907",
  "tauroursodeoxycholic acid_17913"
)

glycine <- c(
  "glycochenodeoxycholic acid_15877",
  "glycochenodeoxycholic acid_15879",
  "glycochenodeoxycholic acid_15880",
  "glycocholic acid_16904",
  "glycohyocholic acid_15732"
)

# assign each group to each bile acids
organ_grouped <- organ_bileacids_df_1 |> 
  dplyr::mutate(Group = case_when(
    Bile_acids %in% glycine    ~ "glycine",
    Bile_acids %in% taurine   ~ "taurine",
    Bile_acids %in% free_bas ~ "free_bas",
    Bile_acids %in% aa_conj_cols              ~ "aa_conj", TRUE ~ "other"))

# counts how many per organs
group_counts <- organ_grouped |> 
  group_by(TissueType, Group) |> 
  summarise(N = n(), .groups = "drop")

# pivot wider
group_counts_wide <- group_counts |> 
  pivot_wider(names_from = Group, values_from = N, values_fill = 0) |> 
  dplyr::select(TissueType, glycine, taurine, free_bas, aa_conj)
#write_csv(group_counts_wide, "Patterson_bile_acids_counts_per_organs.csv")

# Plotting individual bile acids with tissue type combine 
## calculate the median per ZT per tissue type 
df_aa_summary <- table_simplified_red |>  
  group_by(TissueType, ZT) |>  
  summarise(across(all_of(aa_conj_cols),  ~ median(.x, na.rm = TRUE), .names = "{.col}"),.groups = "drop")

# pivot longer
df_aa_summary_long <- df_aa_summary |>  
  pivot_longer(
    cols      = all_of(aa_conj_cols),
    names_to  = "Metabolite",
    values_to = "Abundance")

# Normalize the abundance by the maximum value per tissue type and metabolite after median calculation
df_norm <- df_aa_summary_long |> 
  group_by(TissueType, Metabolite) |> 
  dplyr::mutate(RelAbundance = Abundance / max(Abundance, na.rm = TRUE) * 100) |> 
  ungroup()

#write_csv(df_norm, "Patterson_bile_acids_normalized_per_organs.csv")


# Separate contents (Contents, blood, urine, and feces) from GI organs and distal organs
df_norm_GI_organs <- df_aa_summary_long |> 
  dplyr::filter(TissueType %in% c("Cecum", "Colon", "Duodenum", "Ileum", "Jejunum")) |>
  group_by(TissueType, Metabolite) |> 
  dplyr::mutate(RelAbundance_norm = Abundance / max(Abundance, na.rm = TRUE) * 100, 
                RelAbundance_norm = if_else(is.nan(RelAbundance_norm), 0, RelAbundance_norm)) |> 
  ungroup()

df_norm_distal_organ <- df_aa_summary_long |> 
  dplyr::filter(TissueType %in% c("Blood", "Heart", "Eye", "Brain", "Kidney", "Liver", "Lung", "Skin", "Spleen", "Stomach")) |>
  group_by(TissueType, Metabolite) |> 
  dplyr::mutate(RelAbundance_norm = Abundance / max(Abundance, na.rm = TRUE) * 100, 
                RelAbundance_norm = if_else(is.nan(RelAbundance_norm), 0, RelAbundance_norm)) |> 
  ungroup()

df_norm_contents <- df_aa_summary_long |> 
  dplyr::filter(str_detect(TissueType, "Conts|Cont") | TissueType %in% c("Urine", "Feces")) |>
  group_by(TissueType, Metabolite) |> 
  dplyr::mutate(RelAbundance_norm = Abundance / max(Abundance, na.rm = TRUE) * 100,
                RelAbundance_norm = if_else(is.nan(RelAbundance_norm), 0, RelAbundance_norm)) |> 
  ungroup()

GI_cols <- c(
  Cecum     = "#56B4E9",  # sky blue
  Colon     = "#009E73",  # bluish green
  Duodenum   = "#0072B2",  # blue
  Ileum     = "#F0E442",
  Jejunum = "#000000")


distal_organ_cols <- c(
  Blood     = "#ee2a7b",  # orange
  Heart   = "#CC79A7",  # reddish purple
  Eye  = "#D55E00",  # vermilion
  Brain     = "#E69F00",  # orange
  Kidney = "#999999", 
  Liver = "#56B4E9", 
  Lung = "#009E73", 
  Skin = "#0072B2", 
  Spleen = "#F0E442", 
  Stomach = "#000000"
)

tissue_cols_contents <- c(
  Urine     = "#56B4E9",  # sky blue
  Feces     = "#009E73",  # bluish green
  StomachConts   = "#0072B2",  # blue
  DuodenumConts  = "#D55E00",  # vermilion
  JejunumConts   = "#CC79A7",  # reddish purple
  IlealCont     = "#F0E442",  # yellow
  ColonConts     = "#000000",  # black
  CecalConts     = "#999999"   # grey
)

# Function to apply spline interpolation
safe_spline_y <- function(x, y, xout){
  if (sum(!is.na(y))>=2) {
    spline(x = x, y = y, xout = xout)$y
  } else {
    rep(NA_real_, length(xout))
  }
}

# Create the spline interpolated data for GI organs
df_spline_norm_GI_organs <- df_norm_GI_organs |>
  group_by(TissueType, Metabolite) |>
  summarise(
    newZT = list(seq(min(ZT), max(ZT), length.out = 100)),
    newY  = list(safe_spline_y(ZT, RelAbundance_norm, newZT[[1]])),
    .groups = "drop"
  ) |>
  unnest(c(newZT, newY)) |>
  dplyr::rename(ZT = newZT, RelAbundance_norm = newY)

# Create the spline interpolated data for distal organs
df_spline_norm_distal_organ <- df_norm_distal_organ |>
  group_by(TissueType, Metabolite) |>
  summarise(
    newZT = list(seq(min(ZT), max(ZT), length.out = 100)),
    newY  = list(safe_spline_y(ZT, RelAbundance_norm, newZT[[1]])),
    .groups = "drop"
  ) |>
  unnest(c(newZT, newY)) |>
  dplyr::rename(ZT = newZT, RelAbundance_norm = newY)

# Create the spline interpolated data for contents
df_spline_norm_contents <- df_norm_contents |>
  group_by(TissueType, Metabolite) |>
  summarise(newZT = list(seq(min(ZT), max(ZT), length.out = 100)), newY = list(safe_spline_y(ZT, RelAbundance_norm, newZT[[1]])), .groups = "drop") |>
  unnest(c(newZT, newY)) |>
  dplyr::rename(ZT = newZT, RelAbundance_norm = newY)

# Create the plots for GI organs
df_norm_GI_organs_plot <- ggplot() +
  geom_line(data = df_spline_norm_GI_organs, aes(ZT, RelAbundance_norm, color = TissueType, group = TissueType), linewidth = 1, alpha = 0.8) +
  geom_point(data = df_norm_GI_organs, aes(ZT, RelAbundance_norm, color = TissueType), size   = 3, stroke = NA) +
  scale_color_manual(values = GI_cols) +
  scale_x_continuous(breaks = c(0,4,8,12,16,20), limits = c(0,20)) +
  facet_wrap(~Metabolite, scales="free_y", ncol=4) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.line  = element_line(color="black"), 
        legend.position="right", 
        strip.text=element_text(size=12), 
        axis.text=element_text(size=12),
        strip.background = element_blank(),
        )
df_norm_GI_organs_plot

# Create the plots for distal organs
df_norm_distal_organ_plot <- ggplot() +
  geom_line(data = df_spline_norm_distal_organ, aes(ZT, RelAbundance_norm, color = TissueType, group = TissueType), linewidth = 1, alpha = 0.8) +
  geom_point(data = df_norm_distal_organ, aes(ZT, RelAbundance_norm, color = TissueType), size = 3, stroke = NA) +
  scale_color_manual(values = distal_organ_cols) +
  scale_x_continuous(breaks = c(0,4,8,12,16,20), limits = c(0,20)) +
  facet_wrap(~Metabolite, scales="free_y", ncol=4) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.line  = element_line(color="black"), legend.position="right", strip.text=element_text(size=12), axis.text=element_text(size=12),strip.background = element_blank())
df_norm_distal_organ_plot

# Create the plots for contents
df_norm_contents_plot <- ggplot() +
  geom_line(data = df_spline_norm_contents, aes(ZT, RelAbundance_norm, color = TissueType, group = TissueType), linewidth = 1, alpha     = 0.8) +
  geom_point(data = df_norm_contents, aes(ZT, RelAbundance_norm, color = TissueType), size   = 3, stroke = NA) +
  scale_color_manual(values = tissue_cols_contents) +
  scale_x_continuous(breaks = c(0,4,8,12,16,20), limits = c(0,20)) +
  facet_wrap(~Metabolite, scales="free_y", ncol=4) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.line  = element_line(color="black"), legend.position="right", strip.text=element_text(size=12), axis.text=element_text(size=12),strip.background = element_blank())
df_norm_contents_plot

ggsave(
  filename = "df_norm_GI_organs_plot.pdf",
  plot     = df_norm_GI_organs_plot,
  device   = "pdf",
  width    = 15,
  height   = 10,
  units    = "in")

ggsave(
  filename = "df_norm_distal_organ_plot.pdf",
  plot     = df_norm_distal_organ_plot,
  device   = "pdf",
  width    = 15,
  height   = 10,
  units    = "in")

ggsave(
  filename = "df_norm_contents_plot.pdf",
  plot     = df_norm_contents_plot,
  device   = "pdf",
  width    = 15,
  height   = 10,
  units    = "in")




# separate and plot by organs for individual bile acids (aa conjugates)
organs <- unique(df_norm$TissueType)

plots_by_organ <- setNames(
  lapply(organs, function(org) {
    df_norm |>
      dplyr::filter(TissueType == org) |> 
      ggplot(aes(ZT, RelAbundance, color = Metabolite)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = FALSE) +
      facet_wrap(~ Metabolite, scales = "free_y", ncol = 4) +
      scale_x_continuous(breaks = unique(df_norm$ZT)) +
      labs(title = org, x = "ZT", y = "Abundance") +
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(size = 8),
        axis.text  = element_text(size = 7)
      )
  }),
  nm = organs
)

plots_by_organ

# export individual plots for aa conjugates and then combine them into a single PDF afterwards
for (org in names(plots_by_organ)) {
  ggsave(
    filename = paste0("AAconj_", org, ".pdf"),
    plot     = plots_by_organ[[org]],
    width    = 10,
    height   = 8,
    units    = "in"
  )
}


# Creating UpSet plot of aa conjugates
df_aa_presence <- df2_long |>
  dplyr::filter(!str_detect(TissueType, "Cont") &!TissueType %in% c("Blood", "Urine", "Feces")) |>
  group_by(TissueType, Bile_acid) |>
  summarise(Present = as.integer(any(Present == 1)), .groups = "drop")

df_aa_presence_contents <- df2_long |>
  dplyr::filter(str_detect(TissueType, "Cont") | TissueType %in% c("Blood", "Urine", "Feces")) |>
  group_by(TissueType, Bile_acid) |>
  summarise(Present = as.integer(any(Present == 1)), .groups = "drop")

df_upset <- df_aa_presence |>
  pivot_wider(
    names_from  = TissueType,
    values_from = Present,
    values_fill = 0) |> dplyr::filter(if_any(all_of(tissue_cols), ~ . > 0))

df_upset_contents <- df_aa_presence_contents |>
  pivot_wider(
    names_from  = TissueType,
    values_from = Present,
    values_fill = 0)

tissue_cols <- setdiff(names(df_upset), "Bile_acid")
tissue_cols_contents <- setdiff(names(df_upset_contents), "Bile_acid")

# Create the UpSet plot
upset_plot <- ComplexUpset::upset(
  df_upset,
  intersect        = tissue_cols,
  name             = "Tissue",
  width_ratio      = 0.2,
  min_size         = 1,
  base_annotations = list(
    `Intersection size` = intersection_size(
      counts = TRUE) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()))) +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    text        = element_text(size = 10))
upset_plot


upset_plot_contents <- ComplexUpset::upset(
  df_upset_contents,
  intersect        = tissue_cols_contents,
  name             = "Tissue",
  width_ratio      = 0.2,
  min_size         = 1,
  base_annotations = list(
    `Intersection size` = intersection_size(
      counts = TRUE) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()))) +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    text        = element_text(size = 10))

upset_plot_contents

# save plots
ggsave(
  filename = "upset_plot.pdf",
  plot     = upset_plot,
  device   = "pdf",
  width    = 11,    # inches
  height   = 6,
  units    = "in")

ggsave(
  filename = "upset_plot_contents.pdf",
  plot     = upset_plot_contents,
  device   = "pdf",
  width    = 11,
  height   = 6,
  units    = "in")



 #Explore intersections
# 1) rebuild your bile_intersections if needed:
bile_intersections <- df_upset |>
  pivot_longer(
    cols     = -Bile_acid,
    names_to = "Tissue",
    values_to= "Present"
  ) |>
  filter(Present == 1) |>
  group_by(Bile_acid) |>
  summarise(
   Intersection = paste(sort(Tissue), collapse = "|"),
    .groups     = "drop"
  ) |>
  group_by(Intersection) |>
  summarise(
    Bile_Acids  = list(Bile_acid),
    Count_Acids = n(),
    .groups     = "drop"
  )

bile_intersections_ordered <- bile_intersections |>
  arrange(
    desc(Count_Acids),
    Intersection
  ) |>
  mutate(
    bar_number = row_number()
  )

# 3) inspect the mapping of bar → intersection → acids
bile_intersections_ordered |> 
  select(bar_number, Intersection, Count_Acids, Bile_Acids) |>
  print(n = Inf)

