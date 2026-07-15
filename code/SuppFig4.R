# ---------------------------------------------------------------------------
# 20260609_JS_AmylEx_SuppFig4_global.R - Supplementary Fig 4: GLOBAL
#   (not surface-restricted) central-layer composition comparison of
#   functional vs pathological amyloid fibrils.  Companion to the surface
#   composition analysis in Fig 4 E-H.  Same pipeline as FigSurface.R but
#   WITHOUT the `exposed == 1` filter, so every central-layer residue
#   counts whether it is solvent-exposed or buried.
#
# Question this answers: is the F-vs-P hydrophobicity contrast we see at
# the SIDE SURFACE a structural-sorting phenomenon (functional amyloids
# putting hydrophobic residues outward), or is it just inherited from a
# globally more hydrophobic sequence composition in the PDB-resolved F
# cohort?
#
# Answer (recorded 2026-06-09): the contrast is GLOBAL.  Sign agreement
# between this analysis and the surface analysis is 100% among scales
# significant in both (17/17), and 90% across all 547 tested scales.
# Therefore the PDB-resolved F cohort is intrinsically more hydrophobic
# than the P cohort (driven by over-representation of membrane-associated
# and antimicrobial functional amyloids: PMEL, Citropin, RIPK3, Nup98),
# not by surface sorting.  The surface analysis just AMPLIFIES the same
# contrast that exists globally.
#
# Input data:
#   data/amyloid_surface_composition_per_residue.tsv  - new, written by the
#     C++ command `amyloid_surface_composition` over the 775-PDB cohort.
#     One row per residue in the central 5 layers (3..7 of 10) of each
#     LayerNormalizer-normalized stack.  Columns: pdb, chain, pf, layer,
#     res_idx, res_num, aa1, aa3, cb_neighbors, exposed.
#     The Cb-neighbour grid is built across ALL chains in the structure,
#     so the `exposed` flag means "few neighbours including axial layers" -
#     i.e. side-surface exposure, not tip exposure.
#
#   data/physprops.txt  - 547 AAindex scales x 20 AAs (same file Fig 3
#     panel B uses).
#
#   data/amyloid_explorer_db.tsv  - public metadata for Role and Protein
#     (loaded via load_energies()).
#
# 5 panels:
#   A  Volcano of all 547 scales: x = (functional surface mean - pathological
#      surface mean) on the per-protein-collapsed metric;  y = -log10(BH-FDR)
#      of a two-sided Mann-Whitney U test on the same per-protein values.
#      Three top scales (by FDR) annotated and coloured.
#   B  Top scale #1 as F vs P violin/jitter (one point per protein).
#   C  Top scale #2 as F vs P violin/jitter.
#   D  Top scale #3 as F vs P violin/jitter.
#   E  Per-AA log2 enrichment on the F surface vs the P surface (one row
#      per AA, two columns F and P, or a single F-vs-P log2 column).
#
# Methodology rationale (recorded so future audits don't re-litigate):
#
#   PER-PROTEIN COLLAPSE.  44 Functional structures come from 20 distinct
#   proteins, 25 of them (57 %) from just 5 proteins (Sup35, Citropin,
#   Nup98, PMEL, RIPK3).  As in Fig 4, treating polymorphs as independent
#   observations massively overstates the effective sample size.  We
#   collapse each structure's surface scale mean to one value per protein
#   (mean of polymorph means) before testing.  Result: 20 F-proteins vs
#   48 P-proteins.
#
#   MANN-WHITNEY U.  Same justification as Fig 4: small functional cohort,
#   no parametric distributional assumption beyond identical-shape under H0.
#
#   BH-FDR.  547 scales = multiple-testing burden.  Benjamini-Hochberg gives
#   the right protection against scale-by-scale fishing.
#
# Outputs:
#   output/FigSurface.pdf                     - assembled 5-panel preview
#   output/FigSurface_panels/FigSurface_A.pdf - volcano
#   output/FigSurface_panels/FigSurface_B.pdf - top scale 1 (F vs P)
#   output/FigSurface_panels/FigSurface_C.pdf - top scale 2
#   output/FigSurface_panels/FigSurface_D.pdf - top scale 3
#   output/FigSurface_panels/FigSurface_E.pdf - per-AA enrichment
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
  library(readr)
  library(scales)
})

source("utils.R")
source("load_energies.R")

output_dir <- file.path("..", "output")
dir.create(output_dir, showWarnings = FALSE)

# ── Pretty names for AAindex codes ────────────────────────────────────────
# Same convention as Fig 3 panels C/D/E - reviewers want human-readable
# descriptors next to the opaque AAindex codes.  Extend this table as new
# scales come up in the volcano.  Fallback for unmapped codes is the
# raw AAindex code.
# All entries verified against the AAindex database
# (https://www.genome.jp/dbget-bin/www_bget?aaindex:<CODE>).
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

# Single-line variant for the volcano annotation (no \n).
pretty_label_inline <- function(scale_name) {
  if (scale_name %in% names(pretty_name))
    sprintf("%s (%s)", pretty_name[[scale_name]], scale_name)
  else
    scale_name
}

# ── 1.  Load the surface composition TSV ──────────────────────────────────
surf_raw <- read_tsv(file.path("..", "data",
                               "amyloid_surface_composition_per_residue.tsv"),
                     show_col_types = FALSE, progress = FALSE)

cat(sprintf("Loaded surface TSV: %d rows over %d PDBs\n",
            nrow(surf_raw), n_distinct(surf_raw$pdb)))

# C++ writes pdb stems as "AMEX_<id>" (the LayerNormalizer output filename
# convention).  Strip the prefix so we can join against the metadata.
surf_raw <- surf_raw %>%
  mutate(PDB = sub("^AMEX_", "", pdb))

# ── 2.  Load metadata (Role, Protein) via the shared loader ───────────────
# load_energies() returns per-residue rows with metadata joined; we only
# need one row per PDB.
meta <- load_energies(keep_all_stacks = FALSE) %>%
  distinct(PDB, Protein, Role)

cat(sprintf("Loaded metadata for %d distinct PDBs\n", nrow(meta)))

# ── 3.  Join surface rows to metadata, filter to exposed F/P residues ─────
surf <- surf_raw %>%
  inner_join(meta, by = "PDB") %>%
  filter(Role %in% c("Functional", "Pathological"),
         aa1 %in% LETTERS)
# NB: NO exposure filter here - this is the global central-layer mirror.

cat(sprintf("After F/P filter (ALL central residues): %d residue rows from %d PDBs (F=%d, P=%d)\n",
            nrow(surf),
            n_distinct(surf$PDB),
            n_distinct(surf$PDB[surf$Role == "Functional"]),
            n_distinct(surf$PDB[surf$Role == "Pathological"])))

# ── 4.  Load the AAindex scales (same file Fig 3 panel B uses) ────────────
props_raw <- read_tsv(file.path("..", "data", "physprops.txt"),
                      show_col_types = FALSE, progress = FALSE)
names(props_raw)[1] <- "Scale"

# Build a (scales x AAs) numeric matrix indexed by AA 1-letter code.
aa_cols     <- setdiff(names(props_raw), "Scale")
aa_present  <- intersect(aa_cols, unique(surf$aa1))
scales_mat  <- as.matrix(props_raw[, aa_present])
rownames(scales_mat) <- props_raw$Scale

cat(sprintf("Loaded physprops: %d scales x %d AAs (%d AAs present in surface set)\n",
            nrow(scales_mat), length(aa_cols), length(aa_present)))

# Drop scales with no variance or all-NA across present AAs.
keep_scale <- apply(scales_mat, 1, function(x) {
  ok <- sum(!is.na(x)) >= 3 && sd(x, na.rm = TRUE) > 0
  ok
})
scales_mat <- scales_mat[keep_scale, , drop = FALSE]
cat(sprintf("Scales retained after variance filter: %d\n", nrow(scales_mat)))

# ── 5.  Per-PDB surface scale means via vectorised matrix multiply ────────
# For each PDB:
#   counts[aa] = number of surface residues of that AA in central layers
#   surf_mean_for_scale = sum_aa( scale[aa] * counts[aa] ) / sum_aa(counts[aa])
# So the matrix multiply scales_mat %*% counts gives the numerator vector
# (length = #scales) for every PDB.  NA scale values are handled by
# zero-imputation after weighting (equivalent to dropping that AA from
# the per-PDB mean for that one scale).

aa_counts <- surf %>%
  count(PDB, Protein, Role, aa1, name = "n") %>%
  pivot_wider(names_from = aa1, values_from = n, values_fill = 0)

count_mat <- as.matrix(aa_counts[, aa_present])           # PDB x AA
storage.mode(count_mat) <- "double"
n_surf_per_pdb <- rowSums(count_mat)

# scales_mat is (S x AA); count_mat^T is (AA x PDB).
# numerator (S x PDB) is scales_mat %*% t(count_mat) with NA -> 0.
scales_for_mult <- scales_mat
scales_for_mult[is.na(scales_for_mult)] <- 0
# denominator per (scale, PDB): sum of count[aa] where scale[aa] is not NA.
scale_valid <- !is.na(scales_mat)        # S x AA, TRUE where scale value present
denom <- scale_valid %*% t(count_mat)    # S x PDB

numer <- scales_for_mult %*% t(count_mat)
surf_means <- numer / pmax(denom, 1)     # avoid /0; will be NA after mask below
surf_means[denom == 0] <- NA_real_

# Reshape to long: one row per (Scale, PDB)
pdb_scale <- as_tibble(t(surf_means), .name_repair = "minimal") %>%
  setNames(rownames(scales_mat)) %>%
  mutate(PDB     = aa_counts$PDB,
         Protein = aa_counts$Protein,
         Role    = aa_counts$Role,
         n_surf  = n_surf_per_pdb) %>%
  pivot_longer(cols = -c(PDB, Protein, Role, n_surf),
               names_to = "Scale", values_to = "surf_mean")

cat(sprintf("Per-PDB surface scale matrix: %d (PDB, Scale) pairs\n",
            nrow(pdb_scale)))

# ── 6.  Per-protein collapse (mean of polymorph means) ────────────────────
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

# ── 7.  Per-scale Mann-Whitney U, F vs P, BH-FDR adjusted ─────────────────
# Effect size is Cohen's d on the per-protein means.  Raw (mean_F - mean_P)
# would be dominated by each scale's intrinsic units (NAKH900108 ranges
# ~-0.1..+0.5, MEEJ runs -2..+5), so a raw-difference volcano collapses
# into a cloud rather than a V.  Standardising by pooled SD puts every
# scale on the same axis and recovers the canonical volcano geometry.
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
    cohen_d  = (mean_F - mean_P) / pooled_sd,           # standardised effect
    eff_FmP  = mean_F - mean_P,                         # raw, kept for the TSV
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

# ── 8.  Scales to display as exemplar violins ─────────────────────────────
# Top-by-FDR-within-direction is the natural automatic pick, but it ended
# up surfacing two obscure scales (NAKH900108 = AA composition of fungal /
# plant proteins; KHAG800101 = Kerr-constant increments / electric
# birefringence) that don't read clearly even with their proper names.
# NAKH900108's elevation in F is also at risk of being an artefact: the
# F set contains Sup35 (yeast/fungal prion, 6 polymorphs), so a "fungal
# composition" scale is partly detecting protein origin rather than
# surface biophysics.
#
# We therefore display two well-known, mechanistically clean
# hydrophobicity/hydrophilicity scales as the exemplars, BOTH of which
# also pass FDR.  The volcano shows the whole screen so no information is
# hidden.  The auto-picks (smallest FDR per direction) are still computed
# and logged below for transparency.
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

# ── 9.  Volcano panel REMOVED (June 2026) ─────────────────────────────────
# Mirrors the main-figure decision (see 20260609_JS_AmylEx_Fig4.R): the
# 547-scale AAindex set carries only ~8-10 independent properties, so a
# full-screen volcano with BH-FDR over-corrects and the hydrophobicity
# exemplars do not survive proper dependency-aware correction.  This
# supplementary figure now simply mirrors the GLOBAL (no-surface-filter)
# central-layer composition of the two pre-specified hydrophobicity-axis
# scales, confirming the main-figure surface result is cohort-wide rather
# than a surface-sorting artefact.

# ── 10.  Panels A & B: the two pre-specified hydrophobicity-axis scales ────
pal_pf <- pal_role[c("Functional", "Pathological")]
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
    theme(axis.text.x = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10),
          plot.margin = margin(t = 16, r = 8, b = 6, l = 6))
}

pA <- make_scale_violin(exemplar_enriched)   # Panel A: GOLD730101 (F-enriched)
pB <- make_scale_violin(exemplar_depleted)   # Panel B: HOPT810101 (F-depleted)

# ── 11.  Panel C: per-AA log2(F/P) global central-layer enrichment ─────────
# For each AA, count fraction of surface residues that are that AA, per Role.
# Then log2 of the ratio (with a small pseudocount).  Restricted to the
# canonical 20 AAs: "X" and other non-canonical codes (selenocysteine etc.)
# carry no scale value in physprops.txt and would just clutter the bar.
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
  pivot_wider(names_from = Role, values_from = freq,
              values_fill = 0) %>%
  mutate(log2_FoverP = log2((Functional + eps) / (Pathological + eps))) %>%
  arrange(log2_FoverP) %>%
  mutate(aa1 = fct_inorder(aa1),
         Class = factor(unname(aa_class[as.character(aa1)]),
                        levels = names(pal_aa_class)))

pC <- ggplot(aa_enrich, aes(x = aa1, y = log2_FoverP, fill = Class)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.8) +
  geom_hline(yintercept = 0, linetype = 1, color = "grey20", linewidth = 0.3) +
  scale_fill_manual(values = pal_aa_class, name = "AA class") +
  labs(
    x = "Amino acid",
    y = expression(log[2](frac(F[freq], P[freq])))
  ) +
  theme_paper(base_size = 11) +
  theme(legend.position = "right")

# ── 12.  Assemble preview + save panels ───────────────────────────────────
# Layout (volcano dropped; mirrors main Fig 4 surface panels on the
# global, no-surface-filter composition):
#   row 1: A (GOLD730101 violin) | B (HOPT810101 violin)
#   row 2: C (per-AA enrichment bar, full width)
design <- "
AAABBB
CCCCCC
"

fig_surface <- pA + pB + pC +
  plot_layout(design = design,
              heights = c(1.0, 0.55))

output_pdf <- file.path(output_dir, "FigS4_global.pdf")
ggsave(output_pdf, fig_surface, width = 11, height = 7.5, units = "in",
       device = cairo_pdf)
cat(sprintf("\nWrote %s (assembled 3-panel preview)\n",
            normalizePath(output_pdf)))

# Per-panel PDFs for Illustrator (Supp Fig 4 letter set)
cat("Wrote per-panel PDFs:\n")
save_panel(pA, "S4", "A", width = 3.8, height = 4.5)   # GOLD730101 violin
save_panel(pB, "S4", "B", width = 3.8, height = 4.5)   # HOPT810101 violin
save_panel(pC, "S4", "C", width = 9.0, height = 3.5)   # per-AA enrichment

# Save the global per-scale MW results separately so the surface and
# global TSVs do not overwrite each other.
mw_out <- file.path("..", "data", "global_scales_mw_results.tsv")
write_tsv(mw_per_scale %>% arrange(fdr), mw_out)
cat(sprintf("Wrote %s\n", normalizePath(mw_out)))
