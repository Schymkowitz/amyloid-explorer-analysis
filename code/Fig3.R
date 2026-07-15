# ---------------------------------------------------------------------------
# 20260608_JS_AmylEx_Fig3.R — Residue-level grammar, segmental architecture,
#   and framework polymorphism (Fig 3 in the May 2026 manuscript v2.3 NL).
#
# Numbering history (manuscript was restructured):
#   - revision3 (March 2026): this was Fig 4, titled "Thermodynamic rules
#     of amyloid structural polymorphism".
#   - v2.3 NL (May 2026): renumbered to Fig 3, renamed to "Residue-level
#     grammar, segmental architecture, and framework polymorphism".  Same
#     panel layout (A–I), same underlying data.
#   - The PDF in the manuscript folder (Fig3.pdf) still has Illustrator
#     metadata "Title: Fig4" because the .ai source wasn't renamed — only
#     the manuscript text was updated.
#
# 9 panels A–I, rebuilt against the 775-structure energetics cohort:
#
#   A  Per–amino-acid energetics: violin/box per AA, ordered by Kyte-Doolittle
#      hydropathy + aligned sample-size bar
#   B  Volcano plot: correlation of mean per-AA energy vs ~547 physchem scales
#   C  Representative correlation: hydrophobicity (GOLD730101)
#   D  Representative correlation: side-chain size (GRAR740102, Grantham 1974)
#   E  Representative correlation: aggregation propensity (Aggrescan)
#   F  Stabilizing vs destabilizing run-length distributions (uses
#      Avg_SW_Energy_Stack_1).  Order swap (was Panel G, June 2026):
#      moved before Panel G to match manuscript narrative — the run-length
#      / segmental-architecture subsection precedes the framework-
#      polymorphism subsection in the body text.
#   G  Per-protein energy diversity: SD per protein, bars colored by
#      number-of-structures (was Panel F, June 2026)
#   H  Global quadrant map: per-residue mean energy vs SD across structures
#   I  Heatmap of AA log2-enrichment per quadrant
#
# Source for rewrite: code/20250925 JS AmylEx.R (the 1524-line legacy
#   script — historically labelled "Fig2" then "Fig4", now Fig 3).
#
# Output: paper_figures/output/Fig3.pdf
#
# WORK IN PROGRESS — Panel A only at this stage; B–I added incrementally
# after panel-A review (per Joost's per-file review cadence).
# ---------------------------------------------------------------------------

# ── Run-it-anywhere preamble (same as Fig1) ───────────────────────────────
.find_script_dir <- function() {
  if (sys.nframe() >= 1) {
    for (i in seq_len(sys.nframe())) {
      fr <- sys.frame(i)
      ofile <- tryCatch(fr$ofile, error = function(e) NULL)
      if (!is.null(ofile) && nzchar(ofile) && file.exists(ofile)) {
        return(dirname(normalizePath(ofile)))
      }
    }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) "")
    if (nzchar(p) && file.exists(p)) return(dirname(normalizePath(p)))
  }
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).+", args, perl = TRUE))
  if (length(m) > 0 && file.exists(m[1])) return(dirname(normalizePath(m[1])))
  getwd()
}
setwd(.find_script_dir())
cat(sprintf("Working directory: %s\n", getwd()))

# ── Libraries ─────────────────────────────────────────────────────────────
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(ggrepel)
library(stringr)
library(tibble)
library(patchwork)
library(readr)

# Shared palette + theme + loader
source("utils.R")
source("load_energies.R")

# ── Data ──────────────────────────────────────────────────────────────────
output_dir <- file.path("..", "output")
dir.create(output_dir, showWarnings = FALSE)

# Tidy energies × biology (775 PDBs, ~50k residues, drops AA_1L == "X").
en <- load_energies(keep_all_stacks = FALSE)

# Sanity print so this isn't silently wrong
cat(sprintf(
  "Loaded energies: %d residue rows across %d PDBs, %d proteins\n",
  nrow(en),
  dplyr::n_distinct(en$PDB),
  dplyr::n_distinct(en$Protein[!is.na(en$Protein)])
))

# Energy column used throughout this figure.  Panel F additionally uses
# Avg_SW_Energy_Stack_1 (the solvation-window variant).
energy_col <- "Average_Energy_Stack_1"

# ───────────────────────────────────────────────────────────────────────────
# Panel A — Per-amino-acid energetics
# ───────────────────────────────────────────────────────────────────────────
#
# Violin + inset boxplot of per-residue energy grouped by AA, ordered by
# Kyte-Doolittle hydropathy (hydrophobic → hydrophilic, left → right).
# Lower aligned bar plot shows n per AA.

# Kyte-Doolittle hydropathy scale (kept inline — small, well-known, no
# reason to require physprops.txt for this single panel).
kd <- c(
  I = 4.5, V = 4.2, L = 3.8, F = 2.8, C = 2.5, M = 1.9, A = 1.8,
  G = -0.4, T = -0.7, S = -0.8, W = -0.9, Y = -1.3, P = -1.6,
  H = -3.2, E = -3.5, Q = -3.5, D = -3.5, N = -3.5, K = -3.9, R = -4.5
)

plot_df_A <- en %>%
  mutate(AA_1L = toupper(AA_1L),
         KD = kd[AA_1L]) %>%
  filter(!is.na(.data[[energy_col]]), AA_1L %in% names(kd)) %>%
  mutate(AA_1L = fct_reorder(AA_1L, KD, .desc = TRUE))

n_df_A <- plot_df_A %>%
  count(AA_1L, name = "n")

# Dynamic y-limits from global box-whisker IQR (so whiskers never clip)
whisk_A <- plot_df_A %>%
  group_by(AA_1L) %>%
  summarise(
    q1 = quantile(.data[[energy_col]], 0.25, na.rm = TRUE),
    q3 = quantile(.data[[energy_col]], 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    iqr = q3 - q1,
    upper_w = q3 + 1.5 * iqr,
    lower_w = q1 - 1.5 * iqr
  )
y_upper_A <- max(whisk_A$upper_w, na.rm = TRUE)
y_lower_A <- min(whisk_A$lower_w, na.rm = TRUE)
pad_A     <- 0.05 * (y_upper_A - y_lower_A)
ylims_A   <- c(y_lower_A - pad_A, y_upper_A + pad_A)

# Hydropathy gradient — purposely a continuous 3-stop ramp (not one of the
# semantic categorical palettes in utils.R) because this encodes a
# continuous physical property, not a category.  Endpoints picked to be
# colour-blind safe and to mimic the Dark2 ramp used in the 2025 version.
pal_kd <- c("#1b9e77",   # green   — hydrophilic
            "#d95f02",   # orange  — neutral
            "#7570b3")   # purple  — hydrophobic

p_main_A <- ggplot(plot_df_A,
                   aes(x = AA_1L, y = .data[[energy_col]], fill = KD)) +
  geom_violin(width = 0.9, trim = FALSE, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.18, fill = "white", color = "black",
               linewidth = 0.4, outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", size = 1.3, color = "black") +
  scale_fill_gradientn(colours = pal_kd, name = "Hydropathy (KD)") +
  coord_cartesian(ylim = ylims_A, clip = "on") +
  labs(
    title    = "Per-residue energetics by amino acid",
    subtitle = paste(energy_col, "; AAs ordered by Kyte-Doolittle hydropathy"),
    x = NULL, y = expression(Delta*G~"(kcal/mol)")
  ) +
  theme_paper(base_size = 10) +
  theme(
    axis.text.x   = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  )

p_counts_A <- ggplot(n_df_A, aes(x = AA_1L, y = n)) +
  geom_col(fill = "grey35", width = 0.9) +
  geom_text(aes(label = n), vjust = -0.2, size = 2.6, color = "grey15") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = "Amino acid (hydrophobic → hydrophilic)", y = "Count (n)") +
  theme_paper(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 0, r = 6, b = 6, l = 6)
  )

pA <- p_main_A / p_counts_A + plot_layout(heights = c(3.2, 1))

# ───────────────────────────────────────────────────────────────────────────
# Shared data pipeline for panels B–E
# ───────────────────────────────────────────────────────────────────────────
#
# All four panels rest on the same computation: for each of the ~547
# physicochemical scales in physprops.txt, correlate the mean per-AA
# ΔG against the scale values.  Panel B = volcano (every scale plotted),
# C/D/E = scatter for the three highlighted scales.
#
# Three scales are highlighted in the manuscript figure:
#   - GOLD730101  : hydrophobicity (Goldsack & Chalifoux, 1973)
#   - GRAR740102  : side-chain size  (Grantham, 1974)
#   - Aggrescan   : aggregation propensity

# 1) Mean per-AA energy across all (PDB, residue) rows in the cohort.
aa_energy <- en %>%
  mutate(AA_1L = toupper(AA_1L)) %>%
  filter(!is.na(.data[[energy_col]]), AA_1L %in% names(aa_class)) %>%
  group_by(AA_1L) %>%
  summarise(mean_energy = mean(.data[[energy_col]], na.rm = TRUE),
            n_residues  = dplyr::n(),
            .groups     = "drop")

# 2) physprops.txt — 547 scales × 20 AAs.  Header row's first cell is
# literally "AA"; rename to "Scale" to match the old script convention.
props_raw <- readr::read_tsv(file.path("..", "data", "physprops.txt"),
                             show_col_types = FALSE, progress = FALSE)
names(props_raw)[1] <- "Scale"

aa_cols    <- setdiff(names(props_raw), "Scale")
aa_overlap <- intersect(aa_energy$AA_1L, aa_cols)
cat(sprintf(
  "Loaded physprops: %d scales x %d AAs (%d AAs overlap with cohort)\n",
  nrow(props_raw), length(aa_cols), length(aa_overlap)
))

# Align energy vector to the same AA order used in the scales
e_vec <- aa_energy %>%
  filter(AA_1L %in% aa_overlap) %>%
  arrange(match(AA_1L, aa_overlap)) %>%
  pull(mean_energy)

# 3) Compute Pearson r + p-value for every scale.
cor_df <- props_raw %>%
  rowwise() %>%
  mutate(
    Pearson  = {
      x <- as.numeric(c_across(all_of(aa_overlap)))
      if (sum(!is.na(x)) >= 3 && sd(x, na.rm = TRUE) > 0) {
        suppressWarnings(cor(x, e_vec, use = "pairwise.complete.obs"))
      } else NA_real_
    },
    p_val = {
      x <- as.numeric(c_across(all_of(aa_overlap)))
      if (sum(!is.na(x)) >= 3 && sd(x, na.rm = TRUE) > 0) {
        suppressWarnings(cor.test(x, e_vec)$p.value)
      } else NA_real_
    }
  ) %>%
  ungroup() %>%
  select(Scale, Pearson, p_val) %>%
  mutate(
    fdr           = p.adjust(p_val, method = "BH"),
    neglog10_fdr  = -log10(fdr)
  ) %>%
  filter(!is.na(Pearson))

cat(sprintf(
  "Computed correlations for %d scales (median |r| = %.3f, n_significant_at_FDR_0.05 = %d)\n",
  nrow(cor_df), median(abs(cor_df$Pearson)), sum(cor_df$fdr < 0.05, na.rm = TRUE)
))

# Scales to highlight in panel B (and to scatter in C/D/E).  These
# names must match physprops.txt exactly (case-sensitive).
highlighted_scales <- c("GOLD730101", "GRAR740102", "Aggrescan")

missing_h <- setdiff(highlighted_scales, cor_df$Scale)
if (length(missing_h)) {
  warning(sprintf("Highlighted scales not found in physprops.txt: %s",
                  paste(missing_h, collapse = ", ")))
}

# ───────────────────────────────────────────────────────────────────────────
# Panel B — Volcano: every scale vs mean per-AA energy
# ───────────────────────────────────────────────────────────────────────────

pB <- ggplot(cor_df, aes(x = Pearson, y = neglog10_fdr)) +
  # All scales as grey background points
  geom_point(data = filter(cor_df, !(Scale %in% highlighted_scales)),
             color = "grey75", size = 1.6, alpha = 0.7) +
  # Highlighted scales in vermillion, larger
  geom_point(data = filter(cor_df, Scale %in% highlighted_scales),
             color = unname(okabe_ito["vermillion"]), size = 3.5) +
  # Reference lines: r = 0 (vertical) and FDR = 0.05 (horizontal)
  geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = 3, color = "grey60", linewidth = 0.4) +
  # Label only the three highlighted scales — bigger, italic, off-axis
  ggrepel::geom_text_repel(
    data = filter(cor_df, Scale %in% highlighted_scales),
    aes(label = Scale),
    size = 4.2, fontface = "bold", color = "black",
    min.segment.length = 0, segment.size = 0.4, box.padding = 0.6
  ) +
  labs(
    title    = "Volcano: mean AA energy vs property scales",
    subtitle = sprintf("AAs used (n) = %d, scales tested = %d",
                       length(aa_overlap), nrow(cor_df)),
    x = "Pearson correlation (r)",
    y = expression(-log[10] * "(FDR)")
  ) +
  theme_paper(base_size = 12)

# ───────────────────────────────────────────────────────────────────────────
# Panels C, D, E — scatter for each highlighted scale
# ───────────────────────────────────────────────────────────────────────────
#
# One point per AA, coloured by biochemistry class (pal_aa_class from
# utils.R).  Each panel: title = "Correlation with <Scale>", subtitle =
# Pearson r + p + R².

# Long-format props for joining one scale at a time
props_long <- props_raw %>%
  select(Scale, all_of(aa_overlap)) %>%
  pivot_longer(-Scale, names_to = "AA_1L", values_to = "Property")

make_scale_scatter <- function(scale_name) {
  scatter_df <- aa_energy %>%
    filter(AA_1L %in% aa_overlap) %>%
    inner_join(
      props_long %>% filter(Scale == scale_name) %>% select(AA_1L, Property),
      by = "AA_1L"
    ) %>%
    mutate(Class = factor(unname(aa_class[AA_1L]),
                          levels = names(pal_aa_class)))

  r_row <- cor_df %>% filter(Scale == scale_name)
  r_val <- r_row$Pearson
  p_val <- r_row$p_val
  r2    <- r_val^2

  # Note: per reviewer request, the human-readable name for the AAindex
  # code (e.g. "Goldsack-Chalifoux hydrophobicity" for GOLD730101) should
  # be added in Illustrator at assembly time alongside the code shown
  # below.  Mapping for the record:
  #   GOLD730101 = Goldsack-Chalifoux hydrophobicity
  #   GRAR740102 = Grantham polarity
  #   Aggrescan  = Aggrescan aggregation propensity (already readable)

  ggplot(scatter_df,
         aes(x = Property, y = mean_energy, label = AA_1L, color = Class)) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                fill = "grey80", alpha = 0.25, linewidth = 0.6) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(min.segment.length = 0, size = 4.2,
                             color = "black", fontface = "bold",
                             box.padding = 0.4) +
    scale_color_manual(values = pal_aa_class, name = "AA class") +
    labs(
      x = sprintf("%s value", scale_name),
      y = "Mean per-residue ΔG (kcal/mol)"
    ) +
    theme_paper(base_size = 12) +
    theme(legend.position = "right")
}

pC <- make_scale_scatter("GOLD730101")
pD <- make_scale_scatter("GRAR740102")
pE <- make_scale_scatter("Aggrescan")

# Combined 2x2 panel for B, C, D, E.  Same typography settings as the
# individual panels above (base_size = 12), but patchwork ensures the
# four sub-panels share consistent text scaling without Illustrator
# tuning at assembly time.  Legends on C/D/E are identical (AA class)
# so patchwork's `guides = "collect"` pulls them into a single shared
# legend on the right.
#
# Layout matches the published figure (text dump shows "B C / D E").
pBCDE <- (pB + pC) / (pD + pE) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("B", "C", "D", "E"))) &
  theme(plot.tag = element_text(face = "bold", size = 14),
        legend.position = "right")

# ───────────────────────────────────────────────────────────────────────────
# Shared data pipeline for panels F, H, I (per-residue stats per protein)
# ───────────────────────────────────────────────────────────────────────────

# Per-residue stats across structures (within a protein).  Each row is a
# (Protein, residue position) and the SD is the within-protein, across-
# structures variability at that position.  We only keep residues seen in
# ≥ 2 structures (otherwise SD is undefined).
res_stats <- en %>%
  filter(!is.na(.data[[energy_col]]),
         AA_1L %in% names(aa_class),
         !is.na(Protein), Protein != "") %>%
  group_by(Protein, nA, AA_1L) %>%
  summarise(
    n_struct    = dplyr::n_distinct(PDB),
    mean_energy = mean(.data[[energy_col]], na.rm = TRUE),
    sd_energy   = ifelse(n_struct >= 2,
                         stats::sd(.data[[energy_col]], na.rm = TRUE),
                         NA_real_),
    .groups = "drop"
  ) %>%
  filter(!is.na(sd_energy))

cat(sprintf(
  "Per-residue stats: %d (Protein, residue) rows over %d proteins\n",
  nrow(res_stats), dplyr::n_distinct(res_stats$Protein)
))

# ───────────────────────────────────────────────────────────────────────────
# Panel F — Stabilizing vs destabilizing run-length distributions
# ───────────────────────────────────────────────────────────────────────────
#
# (Was Panel G prior to June 2026.  Promoted to Panel F so the figure-letter
# order matches the manuscript narrative: the segmental-architecture
# subsection — "Cooperative stabilising cores alternate with frustrated
# segments along the primary sequence" — cites this panel, and precedes the
# framework-polymorphism subsection that cites the protein-diversity bars
# now in Panel G.)
#
# For each PDB, take Avg_SW_Energy_Stack_1 (sliding-window energy) ordered
# along the sequence, run-length encode the SIGN, and harvest the lengths
# of all "negative" (stabilizing) and "positive" (destabilizing) runs.
# Then pool across all PDBs and show violin + boxplot per sign.
#
# Avg_SW_Energy_Stack_1 isn't loaded by default — load_energies() with
# keep_all_stacks = FALSE drops Stack_2..6 but keeps both Average_Energy
# and Avg_SW_Energy variants of Stack_1.

run_lengths_by_sign <- function(values) {
  v <- ifelse(values < 0, "neg",
              ifelse(values > 0, "pos", "zero"))
  r <- rle(v)
  data.frame(sign = r$values, len = r$lengths)
}

runs_tbl <- en %>%
  filter(!is.na(Avg_SW_Energy_Stack_1)) %>%
  arrange(PDB, nA) %>%
  group_by(PDB) %>%
  summarise(runs = list(run_lengths_by_sign(Avg_SW_Energy_Stack_1)),
            .groups = "drop") %>%
  unnest(runs) %>%
  filter(sign %in% c("neg", "pos")) %>%
  mutate(sign = factor(sign, levels = c("neg", "pos"),
                       labels = c("Stabilizing\n(ΔG < 0)",
                                  "Destabilizing\n(ΔG > 0)")))

pal_runs <- setNames(
  c(unname(okabe_ito["bluish_green"]),   # stabilizing — green
    unname(okabe_ito["vermillion"])),    # destabilizing — orange-red
  c("Stabilizing\n(ΔG < 0)", "Destabilizing\n(ΔG > 0)")
)

pF <- ggplot(runs_tbl, aes(x = sign, y = len, fill = sign)) +
  geom_violin(trim = FALSE, color = "black", linewidth = 0.4, alpha = 0.85,
              width = 0.9) +
  geom_boxplot(width = 0.15, fill = "white", color = "black",
               linewidth = 0.4, outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", color = "black", size = 1.3) +
  scale_fill_manual(values = pal_runs, guide = "none") +
  coord_cartesian(ylim = c(0, 40), clip = "on") +
  labs(
    x = NULL,
    y = "Run length (residues)"
  ) +
  theme_paper(base_size = 10)

# ───────────────────────────────────────────────────────────────────────────
# Panel G — Protein-level energy diversity
# ───────────────────────────────────────────────────────────────────────────
#
# (Was Panel F prior to June 2026; see Panel F header for the swap rationale.)
#
# One value per protein: the MEAN of the per-residue SDs.  Bars sorted by
# diversity, coloured by # of structures contributing (more structures =
# better-supported diversity estimate).  Annotation: n= structures at bar
# end.

protein_diversity <- res_stats %>%
  group_by(Protein) %>%
  summarise(
    mean_sd       = mean(sd_energy),
    residues_used = dplyr::n(),
    .groups       = "drop"
  )

# Count of structures per protein from the energies table (NOT from
# res_stats, which would undercount because of the ≥2-structures filter)
n_struct_tbl <- en %>%
  filter(!is.na(Protein), Protein != "") %>%
  summarise(n_struct = dplyr::n_distinct(PDB), .by = Protein)

plot_tbl_G <- protein_diversity %>%
  left_join(n_struct_tbl, by = "Protein") %>%
  mutate(Protein_f = fct_reorder(Protein, mean_sd))

pG <- ggplot(plot_tbl_G, aes(x = Protein_f, y = mean_sd, fill = n_struct)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.2) +
  geom_text(aes(label = paste0("n=", n_struct)),
            hjust = -0.15, size = 2.5, color = "black") +
  coord_flip(ylim = c(0, max(plot_tbl_G$mean_sd, na.rm = TRUE) * 1.15),
             clip = "off") +
  scale_fill_viridis_c(name = "Structures (n)", option = "C") +
  labs(
    x = NULL,
    y = "Mean per-residue SD across structures"
  ) +
  theme_paper(base_size = 10) +
  theme(
    axis.text.y     = element_text(size = 7),
    legend.position = "right",
    plot.margin     = margin(t = 6, r = 14, b = 6, l = 6)
  )

# ───────────────────────────────────────────────────────────────────────────
# Panels H & I — Per-residue mean vs SD quadrant map + AA enrichment
# ───────────────────────────────────────────────────────────────────────────
#
# Assign each (Protein, residue) to a quadrant by GLOBAL median cuts of
# mean_energy and sd_energy.  Panel H = the scatter with backgrounds; Panel
# I = log2-enrichment of AAs within each quadrant relative to the pooled
# background, with AAs hierarchically clustered by their enrichment
# profile.

thr_mean_global <- median(res_stats$mean_energy, na.rm = TRUE)
thr_sd_global   <- median(res_stats$sd_energy,   na.rm = TRUE)

res_class <- res_stats %>%
  mutate(
    quadrant = case_when(
      mean_energy <= thr_mean_global & sd_energy <= thr_sd_global ~ "core_stable",
      mean_energy <= thr_mean_global & sd_energy >  thr_sd_global ~ "flexible_stabilizing",
      mean_energy >  thr_mean_global & sd_energy <= thr_sd_global ~ "rigid_destabilizing",
      mean_energy >  thr_mean_global & sd_energy >  thr_sd_global ~ "unstable_destabilizing",
      TRUE                                                         ~ "other"
    ),
    quadrant = factor(quadrant,
                      levels = c("core_stable", "flexible_stabilizing",
                                 "rigid_destabilizing", "unstable_destabilizing"))
  )

cat(sprintf(
  "Quadrant cuts: thr_mean = %.3f, thr_sd = %.3f\n",
  thr_mean_global, thr_sd_global
))

# Rectangles for the four quadrant backgrounds
rects_H <- data.frame(
  xmin     = c(-Inf, -Inf, thr_mean_global, thr_mean_global),
  xmax     = c(thr_mean_global, thr_mean_global, Inf, Inf),
  ymin     = c(-Inf, thr_sd_global, -Inf, thr_sd_global),
  ymax     = c(thr_sd_global, Inf, thr_sd_global, Inf),
  quadrant = factor(c("core_stable", "flexible_stabilizing",
                      "rigid_destabilizing", "unstable_destabilizing"),
                    levels = names(pal_quadrant))
)

pH <- ggplot() +
  # pale quadrant backgrounds
  geom_rect(data = rects_H,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = quadrant),
            alpha = 0.15, color = NA) +
  # all residues coloured by quadrant
  geom_point(data = res_class,
             aes(x = mean_energy, y = sd_energy, color = quadrant),
             alpha = 0.55, size = 1.0, stroke = 0) +
  # median cut lines
  geom_vline(xintercept = thr_mean_global, linetype = 2, color = "grey40",
             linewidth = 0.4) +
  geom_hline(yintercept = thr_sd_global,   linetype = 2, color = "grey40",
             linewidth = 0.4) +
  scale_color_manual(values = pal_quadrant, name = "Category") +
  scale_fill_manual(values = pal_quadrant, guide = "none") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 3)) +
  labs(
    title    = paste("Per-residue", energy_col, ": mean vs SD across structures"),
    subtitle = "Dashed lines indicate median mean-energy and SD boundaries",
    x = "Mean energy across PDBs",
    y = "SD across PDBs"
  ) +
  theme_paper(base_size = 10) +
  theme(legend.position = "right")

# ── Panel I — AA log2-enrichment heatmap ──────────────────────────────────

# Per-quadrant AA counts → frequencies
comp_quad <- res_class %>%
  count(quadrant, AA_1L, name = "n") %>%
  group_by(quadrant) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

# Pooled background frequencies (all quadrants together)
comp_bg <- res_class %>%
  count(AA_1L, name = "n_bg") %>%
  mutate(freq_bg = n_bg / sum(n_bg))

eps <- 1e-9
enrich <- comp_quad %>%
  left_join(comp_bg, by = "AA_1L") %>%
  mutate(log2_enrich = log2((freq + eps) / (freq_bg + eps)))

# Hierarchical clustering of AAs by their log2-enrichment profile across
# the four quadrants — gives the "AAs that behave similarly cluster
# together" reading the manuscript points at.
enrich_mat <- enrich %>%
  select(AA_1L, quadrant, log2_enrich) %>%
  pivot_wider(names_from = quadrant, values_from = log2_enrich,
              values_fill = 0)
mat_num <- as.matrix(enrich_mat[, -1])
rownames(mat_num) <- enrich_mat$AA_1L
hc       <- hclust(dist(mat_num))
aa_order <- enrich_mat$AA_1L[hc$order]

enrich <- enrich %>%
  mutate(AA_1L    = factor(AA_1L, levels = aa_order),
         quadrant = factor(quadrant, levels = names(pal_quadrant)))

pI <- ggplot(enrich,
             aes(x = AA_1L, y = quadrant, fill = log2_enrich)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low      = "#2166ac",   # blue (depleted)
    mid      = "white",
    high     = "#b2182b",   # red (enriched)
    midpoint = 0,
    name     = expression(log[2] * " enrich")
  ) +
  labs(
    title    = "Quadrant-specific AA enrichment vs global background",
    subtitle = "AAs clustered by their enrichment profiles",
    x = "Amino acid (clustered)",
    y = NULL
  ) +
  theme_paper(base_size = 10) +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.ticks      = element_line(linewidth = 0.3),
    legend.position = "right",
    panel.grid      = element_blank()
  )

# ── Write outputs: assembled preview + individual panels ──────────────────
# See utils.R / save_panel() for the convention rationale.
#
# Two outputs per figure:
#   (1) `output/Fig3.pdf`               — assembled patchwork preview
#   (2) `output/Fig3_panels/Fig3_*.pdf` — per-panel PDFs for Illustrator

# (1) Assembled preview — all 9 panels.  Layout reads top-to-bottom in
#     manuscript-narrative order (A → B-E → F → G → H → I):
#       row 1   : A (per-AA violin, wide, takes 4 cols) | B (volcano)
#       row 2   : A continues                            | C (GOLD scatter)
#       row 3   : A continues                            | D (GRAR scatter)
#       row 4   : A continues                            | E (Aggrescan scatter)
#       row 5-6 : F (run-length violins) | H (quadrant scatter, 3x wider)
#       row 7-8 : G (protein-diversity bars, full width)
#       row 9   : I (AA enrichment heatmap, full width)
#     F sits beside H rather than on its own row because the run-length
#     plot is narrow (two groups) and pairs naturally with the wide
#     quadrant scatter; G then gets a full-width row for the ~25 protein
#     bars.  Patchwork's layout is only for the quick-look preview —
#     Illustrator does final assembly from the per-panel PDFs.
design <- "
AAAABB
AAAACC
AAAADD
AAAAEE
FHHHHH
FHHHHH
GGGGGG
GGGGGG
IIIIII
"
fig3_preview <- pA + pB + pC + pD + pE + pF + pG + pH + pI +
  plot_layout(design = design,
              heights = c(1, 1, 1, 1, 1, 1, 1.1, 1.1, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

output_pdf <- file.path(output_dir, "Fig3.pdf")
ggsave(output_pdf, fig3_preview, width = 16, height = 20, units = "in",
       device = cairo_pdf)
cat(sprintf("Wrote %s (assembled preview — all 9 panels A-I)\n",
            normalizePath(output_pdf)))

# (2) Individual panel PDFs for Illustrator.  Sizes tuned to each panel's
#     natural aspect ratio.
#
# IMPORTANT for B/C/D/E: prefer the combined `Fig3_BCDE.pdf` below over
# the individual panel PDFs.  When the four are rendered together by
# patchwork they share typography and one shared AA-class legend, so
# Illustrator assembly is just one drag-and-drop instead of four
# individually-tuned drops.  Individual B/C/D/E saves are kept for
# debugging / panel-by-panel tweaking only.
cat("Wrote per-panel PDFs:\n")
save_panel(pA, 3, "A", width = 8.0, height = 6.0)   # AA violins + sample-size bar
save_panel(pB, 3, "B", width = 5.5, height = 4.5)   # volcano (square-ish)         — superseded by Fig3_BCDE.pdf
save_panel(pC, 3, "C", width = 5.5, height = 4.0)   # GOLD730101 scatter           — superseded by Fig3_BCDE.pdf
save_panel(pD, 3, "D", width = 5.5, height = 4.0)   # GRAR740102 scatter           — superseded by Fig3_BCDE.pdf
save_panel(pE, 3, "E", width = 5.5, height = 4.0)   # Aggrescan scatter            — superseded by Fig3_BCDE.pdf
save_panel(pF, 3, "F", width = 4.0, height = 4.5)   # run-length violins (2 groups)
save_panel(pG, 3, "G", width = 7.0, height = 5.5)   # protein diversity bars (~25 proteins)
save_panel(pH, 3, "H", width = 6.5, height = 5.0)   # quadrant scatter (wide)
save_panel(pI, 3, "I", width = 8.0, height = 3.5)   # heatmap (wide, short)

# Combined B/C/D/E panel — the one to use for Illustrator assembly.
# Larger overall canvas so each of the four sub-panels gets ~6x5 inches,
# matching the area of an individual panel save, with patchwork managing
# text scaling and the shared legend.
save_panel(pBCDE, 3, "BCDE", width = 13.0, height = 11.0)
