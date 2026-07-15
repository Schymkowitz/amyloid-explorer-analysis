# ---------------------------------------------------------------------------
# 20260608_JS_AmylEx_SuppFig1.R — Robustness analysis (Supp Fig 1 in the
#   May 2026 manuscript v2.3 NL).  Manuscript caption: "Conserved
#   thermodynamic principles generalise beyond tau and α-synuclein".
#
# The on-disk reference (paper_figures/old/sFig2.pdf) is the same content
# — the "sFig2" filename predates the May 2026 supplementary renumbering.
# This script outputs under the new manuscript naming (FigS1.pdf).
#
# Four panels:
#   A  Per-amino-acid ΔG distributions for the FULL dataset.  Violin/box
#      per AA ordered by Kyte-Doolittle hydropathy + aligned sample-size
#      bar below.  (Identical recipe to Fig 3 Panel A but kept inline so
#      the supplementary stands on its own.)
#   B  Same distributions for the REDUCED dataset (tau + α-synuclein
#      excluded).  Same layout as A so eyeballing the comparison is easy.
#   C  Side-by-side violin/box of stabilising vs destabilising run lengths,
#      one violin per (run-type × dataset) pair.  Shows that the
#      stabilising/destabilising segmental architecture is essentially
#      unchanged when the two dominant proteins are removed.
#   D  Per-AA median ΔG correlation scatter: each point is one of the 20
#      AAs, x = full-dataset median, y = reduced-dataset median, with a
#      y = x reference line.  Pearson r + Spearman ρ annotated.  Makes
#      the robustness claim quantitative on the figure itself (rather
#      than leaving it implicit in the eyeball-the-distributions test
#      of panels A and B).
#
# Source for rewrite: legacy code/20260324_JS_new_analyses.R "ANALYSIS D"
#   block (lines 417-533).  Original output: fig_D_robustness_tau_asyn_
#   excluded.pdf in old_source_material/2026/.
#
# Outputs:
#   output/FigS1.pdf                       — assembled 4-panel preview
#   output/FigS1_panels/FigS1_A.pdf        — full-dataset AA violins (+ n bar)
#   output/FigS1_panels/FigS1_B.pdf        — tau+α-syn-excluded violins
#   output/FigS1_panels/FigS1_C.pdf        — run-length comparison
#   output/FigS1_panels/FigS1_D.pdf        — per-AA median ΔG correlation
#
# Naming history:
#   - revision3 (March 2026): "Thermodynamic principles are robust to the
#     exclusion of tau and α-synuclein"
#   - v2.3 NL  (May  2026): "Conserved thermodynamic principles generalise
#     beyond tau and α-synuclein" (renamed only; same content)
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
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(stringr)
library(tibble)
library(patchwork)
library(readr)
library(scales)

# Shared palette + theme + loader
source("utils.R")
source("load_energies.R")

# ── Data ──────────────────────────────────────────────────────────────────
output_dir <- file.path("..", "output")
dir.create(output_dir, showWarnings = FALSE)

# Energetics × biology (775 PDBs, ~50k residues, AA_1L == "X" already dropped)
en <- load_energies(keep_all_stacks = FALSE)
energy_col <- "Average_Energy_Stack_1"

# Identify the dominant proteins to exclude.  In the public DB the Greek
# α character is preserved in the Protein column.  Matching via grepl
# avoids the silent-fail trap where an in-source Greek α has an "unknown"
# Encoding() attribute while the DB-loaded string is marked "UTF-8" -
# byte-identical but identical() returns FALSE, so %in% drops the match.
# "synuclein" alone is a safe substring (no other Protein contains it).
# "^Tau$" anchored to avoid matching anything else.
is_dominant <- function(p) grepl("^Tau$", p) |
                            grepl("synuclein", p, ignore.case = TRUE)

n_dominant_pdb <- en %>% filter(is_dominant(Protein)) %>%
  dplyr::pull(PDB) %>% dplyr::n_distinct()
n_total_pdb    <- dplyr::n_distinct(en$PDB)
n_reduced_pdb  <- n_total_pdb - n_dominant_pdb
pct_dominant   <- round(100 * n_dominant_pdb / n_total_pdb)

# Guardrail: tau and α-syn are the two dominant proteins by a wide margin
# in any public-amyloid cohort.  If we matched < 300 PDBs the filter has
# silently failed (this exact bug bit us during the v2.3 rebuild).
stopifnot(
  "tau + α-syn filter matched < 300 PDBs - probable encoding regression" =
    n_dominant_pdb >= 300
)

cat(sprintf(
  "Full dataset: %d PDBs\nTau + α-syn = %d PDBs (%d%% of total)\nReduced dataset: %d PDBs\n",
  n_total_pdb, n_dominant_pdb, pct_dominant, n_reduced_pdb
))

en_reduced <- en %>% filter(!is_dominant(Protein))

# Kyte-Doolittle hydropathy (same scale used by main Fig 3 Panel A)
kd <- c(
  I = 4.5, V = 4.2, L = 3.8, F = 2.8, C = 2.5, M = 1.9, A = 1.8,
  G = -0.4, T = -0.7, S = -0.8, W = -0.9, Y = -1.3, P = -1.6,
  H = -3.2, E = -3.5, Q = -3.5, D = -3.5, N = -3.5, K = -3.9, R = -4.5
)

# Hydropathy gradient (continuous, NOT one of the categorical palettes
# in utils.R — encodes the physical property axis)
pal_kd <- c("#1b9e77",   # green   — hydrophilic
            "#d95f02",   # orange  — neutral
            "#7570b3")   # purple  — hydrophobic

# ── Helper: build a (violin + counts-bar) panel for one dataset ────────────
make_aa_violin_panel <- function(dat, title_str, y_limits = c(-4, 8)) {
  plot_df <- dat %>%
    mutate(AA_1L = toupper(AA_1L), KD = kd[AA_1L]) %>%
    filter(!is.na(.data[[energy_col]]), AA_1L %in% names(kd)) %>%
    mutate(AA_1L = fct_reorder(AA_1L, KD, .desc = TRUE))

  n_df <- plot_df %>% count(AA_1L, name = "n")

  p_main <- ggplot(plot_df,
                   aes(x = AA_1L, y = .data[[energy_col]], fill = KD)) +
    geom_violin(width = 0.9, trim = FALSE, color = "black", linewidth = 0.4) +
    geom_boxplot(width = 0.18, fill = "white", color = "black",
                 linewidth = 0.4, outlier.shape = NA) +
    stat_summary(fun = median, geom = "point", size = 1.3, color = "black") +
    scale_fill_gradientn(colours = pal_kd, name = "Hydropathy (KD)") +
    coord_cartesian(ylim = y_limits, clip = "on") +
    labs(
      x = NULL, y = expression(Delta*G~"(kcal/mol)")
    ) +
    theme_paper(base_size = 11) +
    theme(
      axis.text.x      = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
      legend.position  = "right"
    )

  p_counts <- ggplot(n_df, aes(x = AA_1L, y = n)) +
    geom_col(fill = "grey45", color = "black", linewidth = 0.2, width = 0.85) +
    geom_text(aes(label = comma(n)), vjust = -0.2, size = 2.4, color = "grey15") +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0.02, 0.15))) +
    labs(x = "Amino acid (hydrophobic → hydrophilic)", y = "n") +
    theme_paper(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
      plot.margin = margin(t = 0, r = 6, b = 6, l = 6)
    )

  p_main / p_counts + plot_layout(heights = c(3.2, 1))
}

# ───────────────────────────────────────────────────────────────────────────
# Panel A — Full dataset AA violins
# ───────────────────────────────────────────────────────────────────────────
pA <- make_aa_violin_panel(en,
        sprintf("All proteins (n = %d structures)", n_total_pdb))

# ───────────────────────────────────────────────────────────────────────────
# Panel B — Tau + α-synuclein excluded AA violins
# ───────────────────────────────────────────────────────────────────────────
pB <- make_aa_violin_panel(en_reduced,
        sprintf("Tau + α-syn excluded (n = %d structures)", n_reduced_pdb))

# ───────────────────────────────────────────────────────────────────────────
# Panel C — Run-length comparison: full vs reduced
# ───────────────────────────────────────────────────────────────────────────
#
# For each PDB, walk Avg_SW_Energy_Stack_1 along the sequence and run-
# length encode the sign.  Harvest stabilising (ΔG < 0) and destabilising
# (ΔG > 0) run lengths.  Compare distributions between full and reduced
# datasets side-by-side per run type.

compute_runs <- function(dat, label) {
  dat %>%
    filter(!is.na(Avg_SW_Energy_Stack_1)) %>%
    arrange(PDB, nA) %>%
    group_by(PDB) %>%
    summarise(
      r_obj = list(rle(Avg_SW_Energy_Stack_1 < 0)),
      .groups = "drop"
    ) %>%
    mutate(
      Stabilising   = purrr::map(r_obj, ~ .x$lengths[.x$values == TRUE]),
      Destabilising = purrr::map(r_obj, ~ .x$lengths[.x$values == FALSE])
    ) %>%
    select(PDB, Stabilising, Destabilising) %>%
    pivot_longer(c(Stabilising, Destabilising),
                 names_to = "type", values_to = "lengths") %>%
    unnest(lengths) %>%
    transmute(dataset = label, type, length = lengths)
}

if (!requireNamespace("purrr", quietly = TRUE)) {
  stop("Please install the 'purrr' package (`install.packages('purrr')`)")
}
library(purrr)

runs_full    <- compute_runs(en,         "Full dataset")
runs_reduced <- compute_runs(en_reduced, "Tau+α-syn excluded")
runs_both    <- bind_rows(runs_full, runs_reduced) %>%
  mutate(
    type    = factor(type,    levels = c("Stabilising", "Destabilising")),
    dataset = factor(dataset, levels = c("Full dataset", "Tau+α-syn excluded"))
  )

pal_dataset <- c(
  "Full dataset"       = unname(okabe_ito["blue"]),       # full = blue
  "Tau+α-syn excluded" = unname(okabe_ito["orange"])      # reduced = orange
)

pC <- ggplot(runs_both, aes(x = type, y = length, fill = dataset)) +
  geom_violin(position = position_dodge(0.8), trim = FALSE,
              color = "black", linewidth = 0.3, alpha = 0.7) +
  geom_boxplot(position = position_dodge(0.8), width = 0.12,
               fill = "white", color = "black", linewidth = 0.3,
               outlier.shape = NA) +
  scale_fill_manual(values = pal_dataset, name = "Dataset") +
  coord_cartesian(ylim = c(0, 30)) +
  labs(
    x = NULL, y = "Run length (residues)"
  ) +
  theme_paper(base_size = 11) +
  theme(legend.position = "bottom")

# ───────────────────────────────────────────────────────────────────────────
# Panel D — Per-AA median ΔG correlation: full vs reduced
# ───────────────────────────────────────────────────────────────────────────
#
# Each of the 20 AAs becomes one point.  X = median ΔG in the full
# dataset, Y = median in tau+α-syn-excluded subset.  A y = x reference
# line shows perfect agreement; deviation from it = the AA's median
# shifted when the dominant proteins were removed.
#
# Pearson r and Spearman ρ annotated in the top-left corner — this is
# what makes the robustness claim QUANTITATIVE on the figure.

med_per_aa <- function(dat) {
  dat %>%
    mutate(AA_1L = toupper(AA_1L)) %>%
    filter(!is.na(.data[[energy_col]]), AA_1L %in% names(kd)) %>%
    group_by(AA_1L) %>%
    summarise(median_dG = median(.data[[energy_col]], na.rm = TRUE),
              .groups   = "drop")
}

med_full    <- med_per_aa(en)         %>% rename(median_full    = median_dG)
med_reduced <- med_per_aa(en_reduced) %>% rename(median_reduced = median_dG)

med_compare <- inner_join(med_full, med_reduced, by = "AA_1L") %>%
  mutate(Class = factor(unname(aa_class[AA_1L]),
                        levels = names(pal_aa_class)))

r_pearson  <- cor(med_compare$median_full, med_compare$median_reduced,
                  method = "pearson")
r_spearman <- cor(med_compare$median_full, med_compare$median_reduced,
                  method = "spearman")
cat(sprintf("Per-AA median ΔG full vs reduced: Pearson r = %.4f, Spearman = %.4f\n",
            r_pearson, r_spearman))

# Symmetric axis range so the y = x line is visually meaningful
ax_range <- range(c(med_compare$median_full, med_compare$median_reduced),
                  na.rm = TRUE)
ax_pad   <- 0.08 * diff(ax_range)
ax_lim   <- ax_range + c(-ax_pad, +ax_pad)

pD <- ggplot(med_compare,
             aes(x = median_full, y = median_reduced,
                 label = AA_1L, color = Class)) +
  # y = x reference line (perfect agreement)
  geom_abline(intercept = 0, slope = 1,
              linetype = 2, color = "grey50", linewidth = 0.4) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(min.segment.length = 0, size = 4.2,
                           color = "black", fontface = "bold",
                           box.padding = 0.4) +
  scale_color_manual(values = pal_aa_class, name = "AA class") +
  coord_equal(xlim = ax_lim, ylim = ax_lim, clip = "off") +
  annotate("text",
           x = ax_lim[1] + 0.03 * diff(ax_lim),
           y = ax_lim[2] - 0.03 * diff(ax_lim),
           label = sprintf("Pearson r = %.4f\nSpearman ρ = %.4f",
                           r_pearson, r_spearman),
           hjust = 0, vjust = 1,
           size = 3.6, fontface = "italic") +
  labs(
    x = "Median ΔG, full dataset (kcal/mol)",
    y = "Median ΔG, reduced dataset (kcal/mol)"
  ) +
  theme_paper(base_size = 11) +
  theme(legend.position = "right")

# ───────────────────────────────────────────────────────────────────────────
# Assemble + save (two-track output convention)
# ───────────────────────────────────────────────────────────────────────────

# 2×2 layout: A above B on the left (each is itself a violin+counts
# patchwork → wrap_elements keeps patchwork from flattening their inner
# panels), C above D on the right.  C is the run-length comparison
# (eyeball test); D is the per-AA correlation scatter (quantitative
# robustness summary).
#
#   ┌─────────────────────┬──────────────────┐
#   │  A (full)           │  C (run lengths) │
#   │  violin + n bar     │                  │
#   ├─────────────────────┼──────────────────┤
#   │  B (excluded)       │  D (correlation) │
#   │  violin + n bar     │                  │
#   └─────────────────────┴──────────────────┘

left_col  <- wrap_elements(pA) / wrap_elements(pB)
right_col <- wrap_elements(pC) / wrap_elements(pD)

# No plot title / subtitle: titles and panel letters are added in
# Illustrator at assembly time; figure caption lives in the manuscript.
figS1 <- (left_col | right_col)

# (1) Assembled preview — Supp Fig 1
output_pdf <- file.path(output_dir, "FigS1.pdf")
ggsave(output_pdf, figS1, width = 16, height = 12, units = "in",
       device = cairo_pdf)
cat(sprintf("Wrote %s (assembled 4-panel preview)\n",
            normalizePath(output_pdf)))

# (2) Individual panel PDFs for Illustrator assembly.  A and B are each a
#     two-row patchwork (violin + counts bar), so save_panel() saves them
#     as compound units — Illustrator can treat each as one block.
cat("Wrote per-panel PDFs:\n")
save_panel(pA, "S1", "A", width = 7.0, height = 6.0)   # full-dataset violins + n
save_panel(pB, "S1", "B", width = 7.0, height = 6.0)   # excluded violins + n
save_panel(pC, "S1", "C", width = 6.0, height = 5.5)   # run-length comparison
save_panel(pD, "S1", "D", width = 6.0, height = 5.5)   # per-AA median correlation
