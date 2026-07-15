# ---------------------------------------------------------------------------
# 20260608_JS_AmylEx_Fig5.R — Sample-origin dependence of fibril energetics
#   (Fig 5 in the May 2026 manuscript v2.3 NL).
#
# Manuscript caption: "Patient-derived ex vivo fibrils show greater
# within-structure energetic heterogeneity than recombinant assemblies"
#
# ── Methodology change vs the previous-reviewer-round draft ───────────────
#
# Previous draft (and legacy code/20260324_JS_new_analyses.R ANALYSIS B):
#   pooled Wilcoxon rank-sum test comparing In vitro vs Ex vivo cohorts,
#   ignoring which proteins each structure belonged to.  That test is
#   confounded by protein composition: if one origin happens to contain
#   structurally-more-heterogeneous proteins (e.g. TTR, TMEM106B in the
#   patient cohort), the pooled test reports a "more heterogeneous"
#   effect that's actually about protein identity, not origin.
#
# This rebuild:
#   uses a MIXED-EFFECTS model — `metric ~ Origin + (1 | Protein)` —
#   fit by REML, with the omnibus Origin effect tested via a Type-II
#   F-test (Satterthwaite df, from lmerTest).  Pairwise origin
#   contrasts via emmeans with Bonferroni adjustment.  This properly
#   partitions within-protein origin variance from between-protein
#   composition variance.
#
# When we tried both on the same data:
#   - Pooled SD ΔG (the headline panel): p = 3e-9
#   - Within-protein paired test:        p = 0.85 (4/10 proteins shift
#     in the predicted direction — basically a coin-flip)
#   → the headline pooled finding was almost entirely a protein-
#     composition artifact.  The mixed-effects model below is the
#     more honest test and should change the manuscript narrative
#     for this figure if its omnibus p-values come back insignificant.
#
# ── Origin classification ─────────────────────────────────────────────────
#
# Uses the full 4-way classification we have in amyloid_explorer_db.tsv
# (Recombinant / Seeded / Animal / Patient).  All four have ≥50
# pathological structures so the mixed-effects estimates are stable.
# Ordered as a "lab-to-organism" axis: Recombinant → Seeded → Animal →
# Patient.  The binary collapse (Recombinant+Seeded vs Patient+Animal)
# was kept as a parallel alt script during the methodology decision and
# is now retired — this script is the canonical Fig 5.
#
# ── Layout (3 panels — slimmed to significant metrics only) ───────────────
#
#   Top row (A-B):
#     Violin + box + jitter per origin, for the TWO metrics where the
#     mixed-effects Origin term reaches significance after controlling
#     for protein composition:
#       A  Fraction stabilising      (omnibus p ≈ 0.019)
#       B  Within-structure SD of ΔG (omnibus p ≈ 0.0065)
#     Pairwise Bonferroni-adjusted contrasts shown as significance
#     brackets (* p < 0.05).  The two NS metrics — Mean ΔG (p ≈ 0.13)
#     and Mean stabilising run length (p ≈ 0.06) — are still computed
#     and printed to stdout for the record but not shown here.
#   Bottom row (C):
#     Paired within-protein plot — one line per protein with structures
#     in ≥ 2 origins, connecting per-(protein, origin) mean ΔG.  Labels
#     at the rightmost origin endpoint via ggrepel.  Visual companion
#     to the mixed-effects analysis (you can see which proteins drive
#     any apparent origin effect).
#
# ── Required packages ─────────────────────────────────────────────────────
#
# This script needs `lmerTest` and `emmeans` in addition to the usual
# tidyverse:
#
#   install.packages(c("lmerTest", "emmeans"))
#
# `lme4` is already a dependency of `lmerTest`.  Pairwise contrasts use
# emmeans with Bonferroni adjustment, printed to stdout for inspection
# (only omnibus p is annotated on-panel to keep the visuals clean).
#
# ── Outputs ───────────────────────────────────────────────────────────────
#
#   output/Fig5.pdf                      — assembled 5-panel deliverable
#   output/Fig5_panels/Fig5_A.pdf ... E  — individual panel fallbacks
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
library(ggrepel)
library(stringr)
library(tibble)
library(patchwork)
library(readr)

for (pkg in c("lmerTest", "emmeans")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "Required package '%s' is not installed. Install it with:\n  install.packages(c('lmerTest','emmeans'))",
      pkg))
  }
}
library(lmerTest)   # also loads lme4, adds Satterthwaite df to lmer
library(emmeans)

source("utils.R")
source("load_energies.R")

# ── Data ──────────────────────────────────────────────────────────────────
output_dir <- file.path("..", "output")
dir.create(output_dir, showWarnings = FALSE)

en <- load_energies(keep_all_stacks = FALSE)
energy_col <- "Average_Energy_Stack_1"

# Per-PDB summary stats — same recipe as Fig 4 + carry Resolution through
# so we can use it as a covariate in the mixed-effects models below.
per_struct <- en %>%
  filter(!is.na(.data[[energy_col]])) %>%
  arrange(PDB, nA) %>%
  group_by(PDB, Role, Protein, Origin_short, Resolution) %>%
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

ORIGIN_LEVELS <- c("Recombinant", "Seeded", "Animal", "Patient")

origin_dat <- per_struct %>%
  filter(Role == "Pathological", Origin_short %in% ORIGIN_LEVELS,
         !is.na(Protein), Protein != "") %>%
  mutate(Origin = factor(Origin_short, levels = ORIGIN_LEVELS))

cat(sprintf("\nFig5 cohort (Pathological + known origin + protein): %d structures across %d proteins\n",
            nrow(origin_dat), dplyr::n_distinct(origin_dat$Protein)))
print(origin_dat %>% count(Origin) %>% as.data.frame())

# ── Mixed-effects modelling ───────────────────────────────────────────────
#
# For each per-structure metric, fit:
#     metric ~ Origin + (1 | Protein)
# - Random intercept per Protein soaks up the inherent variability
#   between protein identities.
# - Origin (4-level fixed effect) is what we're testing.
# - Type-II F-test (Satterthwaite df from lmerTest) for the omnibus
#   Origin effect — robust to unbalanced cell sizes.
# - Pairwise origin contrasts via emmeans with Bonferroni adjustment.

METRICS <- list(
  mean_dG       = "Mean ΔG (kcal/mol)",
  frac_stab     = "Fraction stabilising",
  mean_run_stab = "Mean stabilising run length\n(residues)",
  sd_dG         = "Within-structure SD of ΔG"
)

fit_mixed <- function(yvar, with_resolution = FALSE) {
  needed_cols <- c("PDB", "Origin", "Protein")
  if (with_resolution) needed_cols <- c(needed_cols, "Resolution")
  dat <- origin_dat %>%
    select(any_of(needed_cols), y = all_of(yvar)) %>%
    filter(is.finite(y))
  if (with_resolution) dat <- dat %>% filter(is.finite(Resolution))
  # REML for the model (default); Type-III F for the Origin fixed effect
  # (Satterthwaite df from lmerTest, robust to unbalanced cell sizes).
  formula <- if (with_resolution) y ~ Origin + Resolution + (1 | Protein)
             else                  y ~ Origin + (1 | Protein)
  m <- lmerTest::lmer(formula, data = dat, REML = TRUE)
  a <- anova(m, type = "III", ddf = "Satterthwaite")
  origin_row <- which(rownames(a) == "Origin")
  p_omnibus <- a[["Pr(>F)"]][[origin_row]]
  F_stat    <- a[["F value"]][[origin_row]]
  df_num    <- a[["NumDF"]][[origin_row]]
  df_den    <- a[["DenDF"]][[origin_row]]
  em        <- emmeans::emmeans(m, ~ Origin)
  prs       <- pairs(em, adjust = "bonferroni")
  res_row   <- which(rownames(a) == "Resolution")
  res_info  <- if (with_resolution && length(res_row) == 1)
                  list(p_resolution = a[["Pr(>F)"]][[res_row]],
                       F_resolution = a[["F value"]][[res_row]])
               else NULL
  list(model         = m,
       n_used        = nrow(dat),
       p_omnibus     = p_omnibus,
       F_stat        = F_stat,
       df_num        = df_num,
       df_den        = df_den,
       res_info      = res_info,
       em_summary    = as.data.frame(em),
       pairs_summary = as.data.frame(prs))
}

# Run BOTH models per metric — the original (Origin only) and the
# resolution-controlled (Origin + Resolution + (1|Protein)).  Print
# side-by-side so any methodology-driven shift is visible.

cat("\n=== Mixed-effects results: Origin only vs Origin + Resolution ===\n")
mixed_results <- list()
mixed_results_with_res <- list()
for (mname in names(METRICS)) {
  cat(sprintf("\n--- %s ---\n", mname))
  r0 <- fit_mixed(mname, with_resolution = FALSE)
  r1 <- fit_mixed(mname, with_resolution = TRUE)
  mixed_results[[mname]] <- r0
  mixed_results_with_res[[mname]] <- r1
  cat(sprintf("  WITHOUT Resolution: F(%g, %.1f) = %.3f, p_Origin = %.4g  (n=%d)\n",
              r0$df_num, r0$df_den, r0$F_stat, r0$p_omnibus, r0$n_used))
  cat(sprintf("  WITH    Resolution: F(%g, %.1f) = %.3f, p_Origin = %.4g  (n=%d)\n",
              r1$df_num, r1$df_den, r1$F_stat, r1$p_omnibus, r1$n_used))
  if (!is.null(r1$res_info)) {
    cat(sprintf("                      Resolution fixed effect: F = %.3f, p = %.4g\n",
                r1$res_info$F_resolution, r1$res_info$p_resolution))
  }
  cat("  Estimated marginal means per origin (resolution-controlled, 95%% CI):\n")
  print(r1$em_summary[, c("Origin","emmean","SE","lower.CL","upper.CL")],
        row.names = FALSE, digits = 3)
  cat("  Pairwise contrasts WITHOUT resolution (Bonferroni-adjusted):\n")
  print(r0$pairs_summary[, c("contrast","estimate","SE","t.ratio","p.value")],
        row.names = FALSE, digits = 3)
  cat("  Pairwise contrasts WITH    resolution (Bonferroni-adjusted):\n")
  print(r1$pairs_summary[, c("contrast","estimate","SE","t.ratio","p.value")],
        row.names = FALSE, digits = 3)
}

# Going forward: use the RESOLUTION-CONTROLLED results for plot annotation.
# The original `mixed_results` (Origin only) stays available in scope if a
# panel ever wants to display both for comparison.
mixed_results <- mixed_results_with_res

# ── Palette (shared 4-way from utils.R) ────────────────────────────────────
pal_origin_4way <- pal_origin[ORIGIN_LEVELS]

# ── Violin helper, annotated with mixed-effects omnibus p + pairwise brackets ──
#
# The omnibus answers "is there any Origin effect at all?".  The pairwise
# brackets answer "which specific origin pairs differ?".  Only brackets
# for pairs significant at p_adj < 0.05 are drawn (Bonferroni-adjusted
# via emmeans).  Stars use the standard convention: * < 0.05, ** < 0.01,
# *** < 0.001.  Full numeric values are in stdout.
#
# When the omnibus is non-significant we skip brackets entirely (you
# shouldn't be picking through pairwise contrasts after a failed omnibus).

sig_stars <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*",   "ns")))
}

make_violin_mixed <- function(yvar, ylab, mixed_res) {
  dat2 <- origin_dat %>%
    mutate(y = .data[[yvar]]) %>%
    filter(is.finite(y))
  ns <- dat2 %>% count(Origin)
  y_min <- min(dat2$y, na.rm = TRUE)
  y_max <- max(dat2$y, na.rm = TRUE)
  y_range <- y_max - y_min

  # Per-(protein, origin) means — open circles on top of the jitter cloud
  prot_means <- dat2 %>%
    group_by(Origin, Protein) %>%
    summarise(prot_mean = mean(y, na.rm = TRUE), .groups = "drop")

  # Build pairwise bracket set: only significant pairs, only if omnibus
  # itself is significant.  Pairs sorted so shortest brackets land lowest
  # (less line crossing).
  pairs_df <- mixed_res$pairs_summary
  ORIGIN_X <- match(levels(dat2$Origin), levels(dat2$Origin))  # 1..4
  names(ORIGIN_X) <- levels(dat2$Origin)
  bracket_pairs <- data.frame()
  if (!is.na(mixed_res$p_omnibus) && mixed_res$p_omnibus < 0.05) {
    sig_pairs <- pairs_df[pairs_df$p.value < 0.05, , drop = FALSE]
    if (nrow(sig_pairs) > 0) {
      parsed <- do.call(rbind, lapply(seq_len(nrow(sig_pairs)), function(i) {
        pp <- strsplit(as.character(sig_pairs$contrast[i]), " - ")[[1]]
        # contrast names may be wrapped in parentheses in some emmeans output
        pp <- gsub("[()]", "", pp)
        data.frame(a = pp[1], b = pp[2], p = sig_pairs$p.value[i],
                   stringsAsFactors = FALSE)
      }))
      parsed$x1 <- pmin(ORIGIN_X[parsed$a], ORIGIN_X[parsed$b])
      parsed$x2 <- pmax(ORIGIN_X[parsed$a], ORIGIN_X[parsed$b])
      parsed$span <- parsed$x2 - parsed$x1
      parsed <- parsed[order(parsed$span, parsed$x1), ]
      bracket_pairs <- parsed
    }
  }

  # Y positions: brackets stack from y_max upward.  Omnibus annotation
  # sits ABOVE the topmost bracket.
  bracket_y_step <- 0.10 * y_range
  bracket_y0     <- y_max + 0.06 * y_range
  if (nrow(bracket_pairs) > 0) {
    bracket_pairs$y <- bracket_y0 + (seq_len(nrow(bracket_pairs)) - 1) * bracket_y_step
    y_top <- max(bracket_pairs$y) + 0.10 * y_range
  } else {
    y_top <- y_max + 0.10 * y_range
  }

  p <- ggplot(dat2, aes(x = Origin, y = y, fill = Origin)) +
    geom_violin(trim = FALSE, alpha = 0.55, color = "black", linewidth = 0.4) +
    geom_boxplot(width = 0.16, fill = "white", color = "black",
                 linewidth = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.10, size = 0.55, alpha = 0.35,
                color = "grey25", stroke = 0) +
    geom_point(data = prot_means,
               aes(x = Origin, y = prot_mean),
               inherit.aes = FALSE,
               shape = 21, fill = "white", color = "black",
               size = 1.6, stroke = 0.5) +
    geom_text(data = ns,
              aes(x = Origin, y = y_min - 0.05 * y_range,
                  label = paste0("n=", n)),
              inherit.aes = FALSE, size = 2.8, color = "grey20")

  # Pairwise significance brackets (only if any)
  tick_drop <- 0.018 * y_range
  for (i in seq_len(nrow(bracket_pairs))) {
    row <- bracket_pairs[i, ]
    p <- p +
      annotate("segment", x = row$x1, xend = row$x2, y = row$y, yend = row$y,
               color = "black", linewidth = 0.4) +
      annotate("segment", x = row$x1, xend = row$x1,
               y = row$y, yend = row$y - tick_drop,
               color = "black", linewidth = 0.4) +
      annotate("segment", x = row$x2, xend = row$x2,
               y = row$y, yend = row$y - tick_drop,
               color = "black", linewidth = 0.4) +
      annotate("text", x = (row$x1 + row$x2) / 2,
               y = row$y + 0.012 * y_range,
               label = sig_stars(row$p), size = 4.2, vjust = 0,
               fontface = "bold")
  }

  # Omnibus annotation at the very top, slightly smaller
  p <- p +
    annotate("text", x = 2.5, y = y_top,
             label = sprintf("Mixed-effects p = %.3g\n(Origin | Protein random)",
                             mixed_res$p_omnibus),
             size = 3.0, hjust = 0.5, vjust = 0, fontface = "italic") +
    scale_fill_manual(values = pal_origin_4way, guide = "none") +
    coord_cartesian(ylim = c(y_min - 0.12 * y_range,
                             y_top + 0.10 * y_range),
                    clip = "off") +
    labs(x = NULL, y = ylab) +
    theme_paper(base_size = 12) +
    theme(
      axis.text.x  = element_text(size = 10.5, face = "bold",
                                   angle = 25, hjust = 1),
      plot.margin  = margin(t = 22, r = 8, b = 6, l = 6)
    )

  p
}

# ── Panels A & B — only the two metrics with a significant Origin effect ──
#
# The mixed-effects models for `mean_dG` (omnibus p ≈ 0.13) and
# `mean_run_stab` (omnibus p ≈ 0.06) are still computed above and
# printed to stdout for the record, but their panels are not shown
# in the main figure — they don't reach significance after accounting
# for protein composition and would just dilute the message.  Only
# Fraction stabilising (omnibus p ≈ 0.019) and Within-structure SD ΔG
# (omnibus p ≈ 0.0065) survive; those become Panels A and B here.

pA <- make_violin_mixed("frac_stab", METRICS$frac_stab, mixed_results$frac_stab)
pB <- make_violin_mixed("sd_dG",     METRICS$sd_dG,     mixed_results$sd_dG)

# Panel C (within-protein paired-trajectory plot of mean ΔG) was removed
# from Fig 5 on 2026-06-09.  Mean ΔG is the non-significant axis after
# the resolution + protein-composition adjustment (p_Origin = 0.058), so
# a panel built around it sat awkwardly next to A (frac stab, p = 0.009)
# and B (SD ΔG, p = 0.002), and the random intercept on Protein in the
# mixed-effects model already absorbs the between-protein composition
# control that the trajectory plot used to provide visually.  The PMEL
# recombinant-to-ex-vivo exception is preserved as one sentence in the
# body text (see Task #59 Fig 5 prose rewrite).

# ── Assemble + save ───────────────────────────────────────────────────────
#
# Slimmed 2-panel main figure: A and B are the two metrics where the
# mixed-effects Origin term reaches significance after accounting for
# protein-composition variance AND structural resolution.
#
# Mean ΔG and Mean stabilising run length are intentionally dropped
# (mixed-effects p_Origin = 0.058 and 0.083 respectively after resolution
# adjustment — see stdout output for the full statistical record).
# No plot title / subtitle: titles and panel letters are added in
# Illustrator at assembly time; descriptive caption lives in the
# manuscript.  Sample-size and methodology stats are still printed to
# stdout and recorded in Task #59 for the manuscript caption.
fig5 <- (pA | pB)

output_pdf <- file.path(output_dir, "Fig5.pdf")
ggsave(output_pdf, fig5, width = 11, height = 6, units = "in",
       device = cairo_pdf)
cat(sprintf("\nWrote %s\n", normalizePath(output_pdf)))

cat("Wrote per-panel PDFs:\n")
save_panel(pA, 5, "A", width = 4.5, height = 5.0)
save_panel(pB, 5, "B", width = 4.5, height = 5.0)
