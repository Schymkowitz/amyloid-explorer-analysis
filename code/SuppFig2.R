# ---------------------------------------------------------------------------
# 20260609_JS_AmylEx_SuppFig2_resolution_energy.R — Supplementary Figure 2
#   for the Kyriazis et al. (v2.3 NL) manuscript.  Resolution and B-factor
#   relationships with FoldX per-residue ΔG, supporting the existing
#   resolution-vs-clash QC framework (Methods, "Structural quality control"
#   subsection of the manuscript).
#
# 2 panels (clean preliminary version):
#   A  Per-structure mean ΔG vs deposited resolution (cryo-EM Pathological
#      structures only).  Each point = one structure.  Linear fit + 95% CI
#      band, Spearman ρ + p annotated.  Shows that resolution and mean ΔG
#      are correlated at the BETWEEN-structure level.
#   B  Distribution of within-structure Spearman ρ between per-residue
#      B-factor (a per-residue local-resolution proxy) and per-residue ΔG.
#      One value per PDB (≥ 15 residues each).  Median annotated; %% of
#      PDBs with ρ > 0 annotated; Wilcoxon signed-rank test against ρ = 0
#      shows whether the within-structure relationship is consistent
#      across the cohort.
#
# Together: A shows the effect exists between structures; B shows it ALSO
# holds within each structure, controlling for between-structure
# composition.  The within-structure result is the stronger of the two
# because each PDB serves as its own control.
#
# Restricted to cryo-EM Pathological structures (solidNMR and X-ray
# B-factors mean different things and are not directly comparable).
#
# Data dependencies:
#   - paper_figures/data/Final_Structs_AE_energies.csv (energetics via load_energies)
#   - paper_figures/data/amyloid_explorer_db.tsv (metadata, Technique=cryoEM)
#   - paper_figures/data/bfactors_per_residue.tsv (preprocessed earlier)
#
# Output:
#   output/SuppFig2_resolution_energy.pdf
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
library(stringr)
library(tibble)
library(patchwork)
library(readr)
library(scales)

source("utils.R")
source("load_energies.R")

# ── Data ──────────────────────────────────────────────────────────────────
output_dir <- file.path("..", "output")
dir.create(output_dir, showWarnings = FALSE)

en <- load_energies(keep_all_stacks = FALSE)
energy_col <- "Average_Energy_Stack_1"

# Per-PDB summary
per_struct <- en %>%
  filter(!is.na(.data[[energy_col]])) %>%
  group_by(PDB, Protein, Role, Technique, Resolution) %>%
  summarise(mean_dG = mean(.data[[energy_col]], na.rm = TRUE),
            sd_dG   = stats::sd(.data[[energy_col]], na.rm = TRUE),
            n_res   = dplyr::n(),
            .groups = "drop")

# Restrict to all cryo-EM structures with a recorded resolution.  We do
# NOT restrict to a Role (Pathological / Functional / Synthetic): the
# question here is whether the force field responds to structural
# quality, and the physics is identical regardless of biological label.
# Solid-state and solution NMR are excluded because their reported
# resolution metadata is not directly comparable to FSC-based cryo-EM
# resolution.
ce <- per_struct %>%
  filter(Technique == "cryoEM",
         !is.na(Resolution),
         is.finite(Resolution),
         Resolution > 0, Resolution < 15,
         !is.na(mean_dG))

cat(sprintf("\nCohort for Panel A: %d cryo-EM structures (all Roles) with deposited resolution\n",
            nrow(ce)))
print(ce %>% count(Role) %>% as.data.frame())

# ── Panel A — Resolution vs per-structure mean ΔG ─────────────────────────
r_p   <- cor(ce$Resolution, ce$mean_dG, method = "pearson")
r_sp  <- cor(ce$Resolution, ce$mean_dG, method = "spearman")
ct_sp <- suppressWarnings(
  cor.test(ce$Resolution, ce$mean_dG, method = "spearman"))

annot_A <- sprintf("Spearman ρ = %+.3f\nPearson  r = %+.3f\np = %.2g\nn = %d structures",
                    r_sp, r_p, ct_sp$p.value, nrow(ce))

cat(sprintf("\nPanel A: Spearman ρ = %+.3f, Pearson r = %+.3f, p = %.2g\n",
            r_sp, r_p, ct_sp$p.value))

pA <- ggplot(ce, aes(x = Resolution, y = mean_dG)) +
  geom_smooth(method = "lm", color = "black", fill = "grey75",
              linewidth = 0.6, alpha = 0.45) +
  geom_point(shape = 21, size = 2.2, alpha = 0.7,
             fill = unname(okabe_ito["sky_blue"]),
             color = "grey20", stroke = 0.3) +
  annotate("label",
           x = max(ce$Resolution),
           y = min(ce$mean_dG) + 0.04 * diff(range(ce$mean_dG)),
           label = annot_A,
           hjust = 1, vjust = 0, size = 3.5,
           fontface = "italic", label.r = unit(0.1, "lines")) +
  labs(
    x = "Resolution (Å)",
    y = "Mean ΔG (kcal/mol)"
  ) +
  theme_paper(base_size = 12)

# ── Panel B — Within-structure ρ(B-factor, ΔG) distribution ──────────────
bf <- readr::read_tsv(file.path("..", "data", "bfactors_per_residue.tsv"),
                       show_col_types = FALSE) %>%
  mutate(PDB = toupper(PDB))

en_p <- en %>%
  filter(Technique == "cryoEM",
         !is.na(.data[[energy_col]]),
         !is.na(Protein), Protein != "")

joined <- en_p %>%
  inner_join(bf %>% select(PDB, nA, bfactor_mean), by = c("PDB", "nA")) %>%
  filter(bfactor_mean > 0, bfactor_mean < 500)

per_pdb_rho <- joined %>%
  group_by(PDB) %>%
  summarise(
    n_res  = dplyr::n(),
    rho_sp = if (dplyr::n() >= 15)
               suppressWarnings(cor(bfactor_mean,
                                      .data[[energy_col]],
                                      method = "spearman"))
            else NA_real_,
    .groups = "drop"
  ) %>%
  filter(!is.na(rho_sp))

med_rho <- median(per_pdb_rho$rho_sp)
pct_pos <- 100 * mean(per_pdb_rho$rho_sp > 0)
wt      <- wilcox.test(per_pdb_rho$rho_sp, mu = 0, alternative = "two.sided")

cat(sprintf("\nPanel B: median ρ = %+.3f, %% positive = %.1f%%, Wilcoxon p = %.2g, n PDBs = %d\n",
            med_rho, pct_pos, wt$p.value, nrow(per_pdb_rho)))

annot_B <- sprintf("Median ρ = %+.3f\n%.1f%% of PDBs ρ > 0\nWilcoxon p = %.2g\nn = %d PDBs",
                    med_rho, pct_pos, wt$p.value, nrow(per_pdb_rho))

pB <- ggplot(per_pdb_rho, aes(x = rho_sp)) +
  geom_histogram(binwidth = 0.04, fill = unname(okabe_ito["sky_blue"]),
                 color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40",
             linewidth = 0.5) +
  geom_vline(xintercept = med_rho, linetype = "solid",
             color = unname(okabe_ito["vermillion"]), linewidth = 0.9) +
  annotate("label",
           x = 0.85, y = Inf,
           label = annot_B,
           hjust = 1, vjust = 1.4, size = 3.5,
           fontface = "italic", label.r = unit(0.1, "lines")) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(
    x = "Spearman ρ within structure",
    y = "Number of structures"
  ) +
  theme_paper(base_size = 12)

# ── Assemble + save ───────────────────────────────────────────────────────
# No plot title / subtitle: titles and captions are added at Illustrator
# assembly time; figure caption lives in the manuscript document.
figS2 <- (pA | pB)

output_pdf <- file.path(output_dir, "SuppFig2_resolution_energy.pdf")
ggsave(output_pdf, figS2, width = 14, height = 6, units = "in",
       device = cairo_pdf)
cat(sprintf("\nWrote %s\n", normalizePath(output_pdf)))

# ── Per-panel saves for Illustrator ───────────────────────────────────────
save_panel(pA, "S2", "A", width = 6.5, height = 5.5)
save_panel(pB, "S2", "B", width = 6.5, height = 5.5)
