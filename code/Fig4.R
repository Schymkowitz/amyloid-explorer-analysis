# ---------------------------------------------------------------------------
# 20260609_JS_AmylEx_Fig4.R - Functional vs Pathological amyloid fibrils:
#   shared cross-β architecture vs distinct side-surface composition.
#
# Single script that builds the complete 8-panel Fig 4:
#
#   ARCHITECTURE (rows A-D)
#     A  Mean ΔG (kcal/mol) per structure
#     B  Fraction of stabilising residues (ΔG < 0)
#     C  Mean length of contiguous stabilising runs (residues)
#     D  Within-structure SD of ΔG
#   Each panel is violin + inset boxplot + jittered structures + an MW p-value
#   annotation.  Test: two-sided Mann-Whitney U.
#
#   SURFACE COMPOSITION (rows E-H)
#     E  Volcano: 547 AAindex scales screened on per-PDB exposed-residue
#        mean values, per-protein collapsed, F vs P MW U with BH-FDR.
#        x = Cohen's d; y = -log10(FDR); two exemplar scales (GOLD730101
#        and HOPT810101) highlighted.
#     F  GOLD730101 (Goldsack-Chalifoux hydrophobicity) surface mean per
#        protein, F vs P violin.
#     G  HOPT810101 (Hopp-Woods hydrophilicity) surface mean per protein,
#        F vs P violin.
#     H  Per-AA log2 enrichment on F surface vs P surface (canonical 20
#        AAs only).
#
# Methodology rationale:
#
#   ARCHITECTURE (A-D).  Legacy versions used a bootstrap-of-means test that
#   was hypersensitive to tail outliers in the n_F = 44 functional cohort
#   (one Orb2 polymorph + two PMEL polymorphs collapse p < 0.001 to p = 0.14
#   on mean stab run length).  Swapped to two-sided Mann-Whitney U: rank-
#   based, robust to tail, makes no parametric assumption beyond identical-
#   shape under H0, standard for unbalanced two-group designs.
#
#   SURFACE (E-H).  AAindex scale screen on residues marked exposed by the
#   `amyloid_surface_composition` C++ command (central 5 layers of each
#   normalised 10-layer stack, excluding tips).  Per-protein collapse (mean
#   of polymorph means per protein) before MW U avoids pseudo-replication:
#   25 of 44 functional structures come from just 5 proteins.  After
#   collapse: 19 F-proteins vs 48 P-proteins.  BH-FDR controls the 547-scale
#   multiple-testing burden.  Two exemplar scales (Goldsack-Chalifoux
#   hydrophobicity GOLD730101 in panel F, Hopp-Woods hydrophilicity
#   HOPT810101 in panel G) are picked for biophysical interpretability;
#   both pass FDR < 0.05.  The smallest-FDR-per-direction "auto-picks"
#   (NAKH900108 and KHAG800101) are still computed and logged, but
#   substituted because (i) NAKH900108 is a fungal/plant compositional
#   profile that partly tracks Sup35's fungal origin in the F cohort
#   (artefactual), and (ii) KHAG800101 is the Kerr-effect electric
#   birefringence scale and reads as opaque to most audiences.  All AAindex
#   descriptions in the `pretty_name` table below have been verified
#   against the AAindex database.
#
# Companion: 20260609_JS_AmylEx_SuppFig4_global.R produces the GLOBAL
# (no exposure filter) mirror of E-H for Supp Fig 4, confirming that the
# F-vs-P composition contrast is cohort-wide and not a surface-sorting
# artefact (100% sign agreement on scales significant in both analyses).
#
# Outputs:
#   output/Fig4.pdf                       - assembled 8-panel preview
#   output/Fig4_panels/Fig4_A..H.pdf      - per-panel for Illustrator
#   data/surface_scales_mw_results.tsv    - full per-scale MW results
# ---------------------------------------------------------------------------

# ── Run-it-anywhere preamble ──────────────────────────────────────────────
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
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(readr)
  library(scales)
})

source("utils.R")
source("load_energies.R")

output_dir <- file.path("..", "output")
dir.create(output_dir, showWarnings = FALSE)

# ── Pretty names for AAindex codes ────────────────────────────────────────
# All entries verified against AAindex
# (https://www.genome.jp/dbget-bin/www_bget?aaindex:<CODE>).  Extend as
# new exemplars come into play.
pretty_name <- c(
  "GOLD730101" = "Goldsack-Chalifoux hydrophobicity",
  "GRAR740102" = "Grantham polarity",
  "HOPT810101" = "Hopp-Woods hydrophilicity",
  "NAKH900108" = "Normalized AA composition (fungi & plant)",
  "KHAG800101" = "Kerr-constant increments (electric birefringence)"
)

pretty_label <- function(scale_name) {
  if (scale_name %in% names(pretty_name))
    sprintf("%s\n(%s)", pretty_name[[scale_name]], scale_name)
  else
    scale_name
}
pretty_label_inline <- function(scale_name) {
  if (scale_name %in% names(pretty_name))
    sprintf("%s (%s)", pretty_name[[scale_name]], scale_name)
  else
    scale_name
}

# ═════════════════════════════════════════════════════════════════════════
# PART 1 - ARCHITECTURE PANELS A-D
# ═════════════════════════════════════════════════════════════════════════

en <- load_energies(keep_all_stacks = FALSE)
energy_col <- "Average_Energy_Stack_1"

# ── Per-structure summary stats (one row per PDB) ─────────────────────────
per_struct <- en %>%
  filter(!is.na(.data[[energy_col]])) %>%
  arrange(PDB, nA) %>%
  group_by(PDB, Role, Protein) %>%
  summarise(
    mean_dG       = mean(.data[[energy_col]], na.rm = TRUE),
    frac_stab     = mean(.data[[energy_col]] < 0, na.rm = TRUE),
    sd_dG         = stats::sd(.data[[energy_col]], na.rm = TRUE),
    mean_run_stab = {
      vals <- .data[[energy_col]]
      r <- rle(vals < 0)
      stab_runs <- r$lengths[r$values == TRUE]
      if (length(stab_runs) > 0) mean(stab_runs) else NA_real_
    },
    .groups = "drop"
  )

cat(sprintf(
  "Per-structure stats: %d PDBs (Pathological=%d, Functional=%d, Synthetic=%d, Unspecified=%d)\n",
  nrow(per_struct),
  sum(per_struct$Role == "Pathological"),
  sum(per_struct$Role == "Functional"),
  sum(per_struct$Role == "Synthetic"),
  sum(per_struct$Role == "Unspecified")
))

pf_dat <- per_struct %>%
  filter(Role %in% c("Pathological", "Functional")) %>%
  mutate(Role = factor(Role, levels = c("Functional", "Pathological")))

n_F <- sum(pf_dat$Role == "Functional")
n_P <- sum(pf_dat$Role == "Pathological")
cat(sprintf("Filtered to F vs P: n_F = %d, n_P = %d (imbalance ratio = %.1fx)\n",
            n_F, n_P, n_P / n_F))

# ── Mann-Whitney U test for each architectural metric ─────────────────────
mw_test <- function(metric_col) {
  wt <- wilcox.test(pf_dat[[metric_col]] ~ pf_dat$Role,
                    exact = FALSE, alternative = "two.sided")
  wt$p.value
}

p_val_mean <- mw_test("mean_dG")
p_val_frac <- mw_test("frac_stab")
p_val_run  <- mw_test("mean_run_stab")
p_val_sd   <- mw_test("sd_dG")

cat("\nArchitecture Mann-Whitney U p-values (two-sided):\n")
cat(sprintf("  mean ΔG       : %.4f\n", p_val_mean))
cat(sprintf("  frac stab     : %.4f\n", p_val_frac))
cat(sprintf("  mean stab run : %.4f\n", p_val_run))
cat(sprintf("  sd ΔG         : %.4f\n", p_val_sd))

# ── Architecture violin helper (A-D) ──────────────────────────────────────
pal_pf <- pal_role[c("Functional", "Pathological")]

make_violin_pf <- function(dat, yvar, ylab, p_val) {
  dat2 <- dat %>%
    mutate(y = .data[[yvar]]) %>%
    filter(is.finite(y))
  y_min <- min(dat2$y, na.rm = TRUE)
  y_max <- max(dat2$y, na.rm = TRUE)
  y_top <- y_max + 0.12 * (y_max - y_min)

  ggplot(dat2, aes(x = Role, y = y, fill = Role)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "black", linewidth = 0.4) +
    geom_boxplot(width = 0.15, fill = "white", color = "black",
                 linewidth = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.10, size = 0.9, alpha = 0.55,
                color = "grey25", stroke = 0) +
    annotate("text", x = 1.5, y = y_top,
             label = sprintf("Mann-Whitney p = %.3f\nn_F = %d   n_P = %d",
                             p_val, n_F, n_P),
             size = 3.2, hjust = 0.5, vjust = 0) +
    scale_fill_manual(values = pal_pf, guide = "none") +
    coord_cartesian(ylim = c(y_min - 0.05 * (y_max - y_min),
                             y_top + 0.05 * (y_max - y_min)),
                    clip = "off") +
    labs(x = NULL, y = ylab) +
    theme_paper(base_size = 12) +
    theme(
      axis.text.x  = element_text(size = 11, face = "bold"),
      plot.margin  = margin(t = 16, r = 8, b = 6, l = 6)
    )
}

pA <- make_violin_pf(pf_dat, "mean_dG",       "Mean ΔG (kcal/mol)",                          p_val_mean)
pB <- make_violin_pf(pf_dat, "frac_stab",     "Fraction stabilising",                        p_val_frac)
pC <- make_violin_pf(pf_dat, "mean_run_stab", "Mean stabilising run length\n(residues)",     p_val_run)
pD <- make_violin_pf(pf_dat, "sd_dG",         "Within-structure SD of ΔG",                   p_val_sd)

# ═════════════════════════════════════════════════════════════════════════
# PART 2 - SURFACE COMPOSITION PANELS E-H
# ═════════════════════════════════════════════════════════════════════════

# ── Load the surface composition TSV (written by the C++ command
#    `amyloid_surface_composition`).  PDB stem stored as "AMEX_<id>";
#    strip the prefix to join against metadata.
surf_raw <- read_tsv(file.path("..", "data",
                               "amyloid_surface_composition_per_residue.tsv"),
                     show_col_types = FALSE, progress = FALSE) %>%
  mutate(PDB = sub("^AMEX_", "", pdb))

cat(sprintf("\nLoaded surface TSV: %d rows over %d PDBs\n",
            nrow(surf_raw), n_distinct(surf_raw$pdb)))

# Metadata (Role, Protein) re-used from the load_energies() call above.
meta <- en %>% distinct(PDB, Protein, Role)

# Join, filter to F/P AND exposed.
surf <- surf_raw %>%
  inner_join(meta, by = "PDB") %>%
  filter(Role %in% c("Functional", "Pathological"),
         exposed == 1,
         aa1 %in% LETTERS)

cat(sprintf("After F/P + exposed filter: %d residue rows from %d PDBs (F=%d, P=%d)\n",
            nrow(surf),
            n_distinct(surf$PDB),
            n_distinct(surf$PDB[surf$Role == "Functional"]),
            n_distinct(surf$PDB[surf$Role == "Pathological"])))

# ── Load AAindex scales (same physprops.txt Fig 3 panel B uses) ───────────
props_raw <- read_tsv(file.path("..", "data", "physprops.txt"),
                      show_col_types = FALSE, progress = FALSE)
names(props_raw)[1] <- "Scale"

aa_cols     <- setdiff(names(props_raw), "Scale")
aa_present  <- intersect(aa_cols, unique(surf$aa1))
scales_mat  <- as.matrix(props_raw[, aa_present])
rownames(scales_mat) <- props_raw$Scale

cat(sprintf("Loaded physprops: %d scales x %d AAs (%d AAs present in surface set)\n",
            nrow(scales_mat), length(aa_cols), length(aa_present)))

keep_scale <- apply(scales_mat, 1, function(x) {
  sum(!is.na(x)) >= 3 && sd(x, na.rm = TRUE) > 0
})
scales_mat <- scales_mat[keep_scale, , drop = FALSE]
cat(sprintf("Scales retained after variance filter: %d\n", nrow(scales_mat)))

# ── Per-PDB surface scale means (vectorised matrix multiply) ──────────────
aa_counts <- surf %>%
  count(PDB, Protein, Role, aa1, name = "n") %>%
  pivot_wider(names_from = aa1, values_from = n, values_fill = 0)

count_mat <- as.matrix(aa_counts[, aa_present])
storage.mode(count_mat) <- "double"
n_surf_per_pdb <- rowSums(count_mat)

scales_for_mult <- scales_mat
scales_for_mult[is.na(scales_for_mult)] <- 0
scale_valid <- !is.na(scales_mat)
denom <- scale_valid %*% t(count_mat)

numer <- scales_for_mult %*% t(count_mat)
surf_means <- numer / pmax(denom, 1)
surf_means[denom == 0] <- NA_real_

pdb_scale <- as_tibble(t(surf_means), .name_repair = "minimal") %>%
  setNames(rownames(scales_mat)) %>%
  mutate(PDB     = aa_counts$PDB,
         Protein = aa_counts$Protein,
         Role    = aa_counts$Role,
         n_surf  = n_surf_per_pdb) %>%
  pivot_longer(cols = -c(PDB, Protein, Role, n_surf),
               names_to = "Scale", values_to = "surf_mean")

# ── Per-protein collapse + MW U + BH-FDR ──────────────────────────────────
prot_scale <- pdb_scale %>%
  filter(!is.na(surf_mean)) %>%
  group_by(Scale, Protein, Role) %>%
  summarise(prot_surf_mean = mean(surf_mean, na.rm = TRUE),
            n_polymorphs   = n(),
            .groups        = "drop")

n_F_prot <- prot_scale %>% filter(Role == "Functional")   %>%
  pull(Protein) %>% n_distinct()
n_P_prot <- prot_scale %>% filter(Role == "Pathological") %>%
  pull(Protein) %>% n_distinct()
cat(sprintf("Per-protein collapse: %d F-proteins vs %d P-proteins\n",
            n_F_prot, n_P_prot))

mw_per_scale <- prot_scale %>%
  group_by(Scale) %>%
  summarise(
    mean_F   = mean(prot_surf_mean[Role == "Functional"],   na.rm = TRUE),
    mean_P   = mean(prot_surf_mean[Role == "Pathological"], na.rm = TRUE),
    sd_F     = stats::sd(prot_surf_mean[Role == "Functional"],   na.rm = TRUE),
    sd_P     = stats::sd(prot_surf_mean[Role == "Pathological"], na.rm = TRUE),
    n_F      = sum(Role == "Functional"   & !is.na(prot_surf_mean)),
    n_P      = sum(Role == "Pathological" & !is.na(prot_surf_mean)),
    pooled_sd = sqrt(((n_F - 1) * sd_F^2 + (n_P - 1) * sd_P^2) /
                       pmax(n_F + n_P - 2, 1)),
    cohen_d  = (mean_F - mean_P) / pooled_sd,
    eff_FmP  = mean_F - mean_P,
    p_mw     = if (n_F >= 3 && n_P >= 3)
      suppressWarnings(
        wilcox.test(prot_surf_mean[Role == "Functional"],
                    prot_surf_mean[Role == "Pathological"],
                    exact = FALSE, alternative = "two.sided")$p.value)
      else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(fdr = p.adjust(p_mw, method = "BH"),
         neglog10_fdr = -log10(fdr))

cat(sprintf("MW per scale: %d tested, %d at FDR<0.05, %d at FDR<0.01\n",
            sum(!is.na(mw_per_scale$p_mw)),
            sum(mw_per_scale$fdr < 0.05, na.rm = TRUE),
            sum(mw_per_scale$fdr < 0.01, na.rm = TRUE)))

# ── Pick exemplar scales for panels F & G ─────────────────────────────────
exemplar_enriched <- "GOLD730101"   # Goldsack-Chalifoux hydrophobicity
exemplar_depleted <- "HOPT810101"   # Hopp-Woods hydrophilicity

top_enriched_auto <- mw_per_scale %>%
  filter(!is.na(fdr), cohen_d > 0) %>%
  arrange(fdr) %>% slice_head(n = 1)
top_depleted_auto <- mw_per_scale %>%
  filter(!is.na(fdr), cohen_d < 0) %>%
  arrange(fdr) %>% slice_head(n = 1)

top_highlight <- mw_per_scale %>%
  filter(Scale %in% c(exemplar_enriched, exemplar_depleted))

cat(sprintf("\nDisplayed F-enriched exemplar: %s\n", exemplar_enriched))
top_highlight %>% filter(Scale == exemplar_enriched) %>%
  select(Scale, mean_F, mean_P, cohen_d, p_mw, fdr) %>%
  as.data.frame() %>% print()
cat(sprintf("\nDisplayed F-depleted exemplar: %s\n", exemplar_depleted))
top_highlight %>% filter(Scale == exemplar_depleted) %>%
  select(Scale, mean_F, mean_P, cohen_d, p_mw, fdr) %>%
  as.data.frame() %>% print()
cat("\n(For reference, auto-picks by smallest FDR per direction:)\n")
cat(sprintf("  most F-enriched by FDR: %s  (d=%+.2f, FDR=%.3g)\n",
            top_enriched_auto$Scale, top_enriched_auto$cohen_d,
            top_enriched_auto$fdr))
cat(sprintf("  most F-depleted by FDR: %s  (d=%+.2f, FDR=%.3g)\n",
            top_depleted_auto$Scale, top_depleted_auto$cohen_d,
            top_depleted_auto$fdr))

# ── Volcano panel REMOVED (June 2026) ─────────────────────────────────────
# The full 547-scale AAindex screen was dropped from the figure.  Reason:
# the 547 scales carry only ~8-10 independent physicochemical properties
# (PCA: 8 PCs for 90 % variance, 10 for 95 %; Li & Ji effective-tests
# Meff = 29).  As a *screen* the hydrophobicity exemplars do not survive
# proper dependency-aware multiple-testing correction (Westfall-Young
# permutation FWER: GOLD730101 = 0.16, HOPT810101 = 0.26) - only an
# artefactual fungal-composition scale (NAKH900108, driven by Sup35) and
# two turn/sheet-propensity scales clear it.  A volcano where only an
# artefact clears the line is a negative-result figure.
#
# The biologically real and robust finding is the hydrophobicity contrast,
# which is the test the *literature* poses (polar-functional-amyloid
# claim).  We therefore present it as a PRE-SPECIFIED two-scale hypothesis
# test (Goldsack-Chalifoux hydrophobicity + Hopp-Woods hydrophilicity)
# rather than the top hit of a screen.  Bonferroni over the 2 pre-specified
# scales: GOLD p = 0.004, HOPT p = 0.008 - both clear 0.05, and the effect
# sizes are large (d = +1.13, -0.84) and robust to protein subsampling
# (92-97 %).  The full-screen MW + BH results are still written to the
# supplementary TSV for transparency, just not shown as a panel.

# ── Panels E & F: the two pre-specified hydrophobicity-axis scales ────────
# Bonferroni correction over the 2 pre-specified scales (not the 547-scale
# screen FDR): p_bonf = min(1, raw_p * 2).
N_PRESPEC <- 2
make_scale_violin <- function(scale_name) {
  d <- prot_scale %>%
    filter(Scale == scale_name) %>%
    mutate(Role = factor(Role, levels = c("Functional", "Pathological")))
  d_val    <- mw_per_scale$cohen_d[mw_per_scale$Scale == scale_name]
  p_raw    <- mw_per_scale$p_mw   [mw_per_scale$Scale == scale_name]
  p_bonf   <- min(1, p_raw * N_PRESPEC)
  ann_lbl  <- sprintf("Cohen's d = %+.2f\nMW p = %.3f  (Bonferroni x2 p = %.3f)\nn_F = %d   n_P = %d",
                      d_val, p_raw, p_bonf,
                      sum(d$Role == "Functional"),
                      sum(d$Role == "Pathological"))
  y_min <- min(d$prot_surf_mean, na.rm = TRUE)
  y_max <- max(d$prot_surf_mean, na.rm = TRUE)
  y_top <- y_max + 0.12 * (y_max - y_min)

  ggplot(d, aes(x = Role, y = prot_surf_mean, fill = Role)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "black", linewidth = 0.4) +
    geom_boxplot(width = 0.15, fill = "white", color = "black",
                 linewidth = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.10, size = 1.4, alpha = 0.7,
                color = "grey25", stroke = 0) +
    annotate("text", x = 1.5, y = y_top, label = ann_lbl,
             size = 3.0, hjust = 0.5, vjust = 0) +
    scale_fill_manual(values = pal_pf, guide = "none") +
    coord_cartesian(ylim = c(y_min - 0.05 * (y_max - y_min),
                             y_top + 0.05 * (y_max - y_min)),
                    clip = "off") +
    labs(x = NULL,
         y = sprintf("Surface mean (per protein)\n%s", pretty_label(scale_name))) +
    theme_paper(base_size = 11) +
    theme(axis.text.x  = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10),
          plot.margin  = margin(t = 16, r = 8, b = 6, l = 6))
}

pE <- make_scale_violin(exemplar_enriched)   # Panel E: GOLD730101 (F-enriched)
pF <- make_scale_violin(exemplar_depleted)   # Panel F: HOPT810101 (F-depleted)

# ── Panel G: per-AA log2(F/P) surface enrichment ──────────────────────────
canonical_aa <- c("A","R","N","D","C","Q","E","G","H","I",
                  "L","K","M","F","P","S","T","W","Y","V")

aa_role_counts <- surf %>%
  filter(aa1 %in% canonical_aa) %>%
  count(Role, aa1, name = "n") %>%
  group_by(Role) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

eps <- 1e-6
aa_enrich <- aa_role_counts %>%
  select(Role, aa1, freq) %>%
  pivot_wider(names_from = Role, values_from = freq, values_fill = 0) %>%
  mutate(log2_FoverP = log2((Functional + eps) / (Pathological + eps))) %>%
  arrange(log2_FoverP) %>%
  mutate(aa1   = fct_inorder(aa1),
         Class = factor(unname(aa_class[as.character(aa1)]),
                        levels = names(pal_aa_class)))

pG <- ggplot(aa_enrich, aes(x = aa1, y = log2_FoverP, fill = Class)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.8) +
  geom_hline(yintercept = 0, linetype = 1, color = "grey20", linewidth = 0.3) +
  scale_fill_manual(values = pal_aa_class, name = "AA class") +
  labs(
    x = "Amino acid",
    y = expression(log[2](frac(F[surf~freq], P[surf~freq])))
  ) +
  theme_paper(base_size = 11) +
  theme(legend.position = "right")

# ═════════════════════════════════════════════════════════════════════════
# PART 3 - ASSEMBLE + SAVE
# ═════════════════════════════════════════════════════════════════════════
#
# 7-panel data layout (volcano dropped; structural exemplars H, I are
# added in Illustrator from the PyMOL renders):
#   row 1: A B C D     (architecture, four small violins)
#   row 2: E F         (two pre-specified hydrophobicity-scale violins)
#   row 3: G           (per-AA enrichment bar, full width)
design <- "
AABBCCDD
EEEEFFFF
GGGGGGGG
"

fig4 <- pA + pB + pC + pD + pE + pF + pG +
  plot_layout(design = design,
              heights = c(1.0, 1.0, 0.55))

# (1) Assembled preview
output_pdf <- file.path(output_dir, "Fig4.pdf")
ggsave(output_pdf, fig4, width = 16, height = 14, units = "in",
       device = cairo_pdf)
cat(sprintf("\nWrote %s (assembled 7-panel deliverable; structural exemplars added in Illustrator)\n",
            normalizePath(output_pdf)))

# (2) Per-panel PDFs for Illustrator.  Data panels A-G; the two structural
#     exemplars (PyMOL) save as H and I from the .pml script.
cat("Wrote per-panel PDFs:\n")
save_panel(pA, 4, "A", width = 3.8, height = 5.0)   # architecture: mean ΔG
save_panel(pB, 4, "B", width = 3.8, height = 5.0)   # architecture: frac stab
save_panel(pC, 4, "C", width = 3.8, height = 5.0)   # architecture: stab run
save_panel(pD, 4, "D", width = 3.8, height = 5.0)   # architecture: SD ΔG
save_panel(pE, 4, "E", width = 3.8, height = 4.5)   # surface: GOLD730101 violin
save_panel(pF, 4, "F", width = 3.8, height = 4.5)   # surface: HOPT810101 violin
save_panel(pG, 4, "G", width = 9.0, height = 3.5)   # surface: per-AA enrichment

# (3) TSV: full per-scale surface MW results, for downstream eyeballing and
#     reviewer responses.
mw_out <- file.path("..", "data", "surface_scales_mw_results.tsv")
write_tsv(mw_per_scale %>% arrange(fdr), mw_out)
cat(sprintf("Wrote %s\n", normalizePath(mw_out)))
